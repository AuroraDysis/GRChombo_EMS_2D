# Run from EMS.jl: julia --project=test /path/to/this/file /path/to/output
# All output is in the fork; EMS.jl and its production tables are read-only.
using EMS, DoubleFloats, SHA

function fixtures(out::String)::Nothing
    mkpath(out)
    T = Double64
    rn = reissner_nordstrom(T; q=T("0.700784945779"), Q=T("0.700784945779"), degree=32, interior_degree=24)
    # phi=0 is the exact RN branch of this same quadratic coupling.
    c = EMSCoupling(4T(pi), zero(T), zero(T), -T(20))
    sol = Solution(rn.exterior, rn.interior, rn.e, c, rn.rh, rn.rexc)
    p, checks = trumpet_profile(sol)
    write_trumpet(joinpath(out, "rn.trumpet"), p)
    println("RN profile complete ", p.meta)
    open(joinpath(out, "manifest.sha256"), "w") do io
        println(io, bytes2hex(sha256(read(joinpath(out, "rn.trumpet")))), "  rn.trumpet")
    end
    nothing
end
fixtures(abspath(ARGS[1]))
