using EMS, SHA, Printf

const SOURCE = "/Users/auroradysis/Workspace/EMS/artifacts/echo/t4"
const DEST = joinpath(@__DIR__, "fixtures-echo")
fmt(x) = @sprintf("%.17g", Float64(x))

function table(dir, name, header, rows)
    open(joinpath(dir, name), "w") do io
        println(io, join(header, '\t'))
        for row in rows
            println(io, join(fmt.(row), '\t'))
        end
    end
end

for member in ("B", "E")
    dir = joinpath(DEST, member)
    mkpath(dir)
    source = joinpath(SOURCE, "$member.trumpet")
    cp(source, joinpath(dir, "reference.trumpet"); force=true)
    p = read_trumpet(source)
    Rh = parse(Float64, p.meta["R_h"])
    radii = (1e-5, 1e-4, Rh / 2, Rh, 2Rh, 0.01, 0.05, 0.16, 0.5, 1.0, 5.0, 20.0)
    ss = (0.0, 1e-6, 0.01, 0.1, 0.5, 0.9, 0.9999, 1.0)
    etas = (0.0, 0.20273255, -0.20273255)
    table(dir, "jets.tsv", vcat(["s"], ["$(name)_d$j" for name in ("Y", "D", "S") for j in 0:3]),
        ([s; [EMS.trumpet_jet(p, k, s, j) for k in 1:3 for j in 0:3]] for s in ss))
    table(dir, "radial.tsv", ["R", "s", "r", "delta", "phi", "X", "alphaK", "betaR",
        "k", "Pi", "ER", "PR", "alphaR", "betaR_deriv"],
        ([R; (v = trumpet_sample(p, R); [v.s, v.r, v.δ, v.φ, v.X, v.α, v.β,
            v.k, v.Π, v.ER, v.PR, v.αR, v.βR])] for R in radii))
    function objectrow(L, eta, R)
        x, y, z = (L * R) .* (0.82, 0.50, 0.28)
        o = trumpet_object(p, x, y, z; mass=L, rapidity=eta)
        [L, eta, R, x, y, z, vec(o.gamma)..., vec(o.K)...,
            o.shift..., o.phi, o.Pi, o.E..., o.B..., o.lapse, o.chi, o.J, o.w]
    end
    table(dir, "objects.tsv", vcat(["L", "eta", "R", "x", "y", "z"],
        ["gamma_$i$j" for j in 1:3 for i in 1:3],
        ["K_$i$j" for j in 1:3 for i in 1:3],
        ["shift_$i" for i in 1:3], ["phi", "Pi"],
        ["E_$i" for i in 1:3], ["B_$i" for i in 1:3],
        ["lapse", "chi", "J", "w"]),
        (objectrow(L, eta, R) for L in (1.0, 2.0), eta in etas, R in radii))
    keyscc = keys(trumpet_ccz4(p, 1.0, 0.3))
    table(dir, "single_ccz4.tsv", vcat(["L", "eta", "R", "x", "y"], String.(keyscc)),
        ([L, eta, R, L * R * 0.82, L * R * 0.57,
            Tuple(trumpet_ccz4(p, L * R * 0.82, L * R * 0.57; mass=L, rapidity=eta))...]
         for L in (1.0, 2.0), eta in etas, R in radii))
    open(joinpath(dir, "manifest.txt"), "w") do io
        println(io, "schema=EMSTRUMPET/1")
        println(io, "generator=Tests/EMSTrumpet/echo_fixtures.jl")
        println(io, "source=$source")
        println(io, "source.sha256=$(bytes2hex(sha256(read(source))))")
        for name in ("src/trumpet.jl", "scripts/emstrumpet_fixtures.jl")
            path = joinpath(dirname(SOURCE), "..", "..", name)
            println(io, "julia.$name.sha256=$(bytes2hex(sha256(read(path))))")
        end
        println(io, "julia_version=$(VERSION)")
        println(io, "radii=$(join(fmt.(radii), ','))")
        println(io, "rapidities=$(join(fmt.(etas), ','))")
        for name in ("reference.trumpet", "jets.tsv", "radial.tsv", "objects.tsv", "single_ccz4.tsv")
            path = joinpath(dir, name)
            println(io, "fixture.$name.bytes=$(filesize(path))")
            println(io, "fixture.$name.sha256=$(bytes2hex(sha256(read(path))))")
        end
    end
    println("$member: ell=$(p.meta["ell"]) Rh=$Rh rows=$(length(radii))")
end
