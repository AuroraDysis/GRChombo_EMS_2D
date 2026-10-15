using EMS, LinearAlgebra, ForwardDiff, DelimitedFiles
BLAS.set_num_threads(1)
const PROFILE=ARGS[1]
const OUT=ARGS[2]
const DATA=ARGS[3]
function packed(c,x;panel=0)
    a=EMS.ctt_pulse_fields(c,x;panel)
    vcat(vec(a.gamma),vec(a.K),a.phi,a.Pi,a.E,a.B,a.shift)
end
function main()
    points=Tuple{Vector{Float64},Int}[]
    c=EMS.read_ctt_pulse(joinpath(DATA,"n24-final.ctt"),PROFILE)
    for (j,p) in enumerate(c.panels)
        r=p.kind==:log ? sqrt(p.lo*p.hi) : 2/p.hi
        for n in ([.37,.48,sqrt(1-.37^2-.48^2)],[sqrt(1-1e-12),1e-6,0.],[-1.,0.,0.])
            push!(points,(r*n,0))
        end
        if j<length(c.panels)
            r=p.hi
            for panel in (0,j,j+1)
                push!(points,(r*[.6,.48,.64],panel))
            end
        end
    end
    for r in (1e-18,1e-10,.0071,.02,1.,16.72433994232569,17.,20.,30.,39.00182260142435,45.,55.,61.,62.06455126251464,128.,256.,1e4)
        push!(points,(r*[.6,.8,0.],0))
    end
    open(joinpath(OUT,"points.txt"),"w") do io
        for (x,p) in points;println(io,join(x,' '), ' ',p);end
    end
    for name in ("zero","half","n12","n18","n24")
        c=EMS.read_ctt_pulse(joinpath(DATA,name*"-final.ctt"),PROFILE)
        # n12/n18 panels have the same geometry; forced index remains meaningful.
        values=reduce(vcat,transpose(packed(c,x;panel=p)) for (x,p) in points)
        writedlm(joinpath(OUT,name*"-julia.txt"),values)
    end
    # Exact AD jets at off-grid points and every radial interface side.
    # Local Cartesian coordinate u=(x-x0)/rho0; regularization is fixed at x0.
    c=EMS.read_ctt_pulse(joinpath(DATA,"n24-final.ctt"),PROFILE)
    jets=Tuple{Vector{Float64},Int,Int}[]
    for j in eachindex(c.panels)
        p=c.panels[j]
        r=p.kind==:log ? sqrt(p.lo*p.hi) : 2/p.hi
        for axis in 1:6;push!(jets,(r*[.6,.48,.64],0,axis));end
        if j<length(c.panels)
            for panel in (0,j,j+1), axis in 1:6
                push!(jets,(p.hi*[.6,.48,.64],panel,axis))
            end
        end
    end
    for r in (1e-10,.0071,20.,39.,55.)
        for axis in 1:6;push!(jets,(r*[sqrt(1-1e-12),1e-6,0.],0,axis));end
    end
    input=open(joinpath(OUT,"stencil-points.txt"),"w")
    reference=open(joinpath(OUT,"stencil-julia.txt"),"w")
    scales=open(joinpath(OUT,"jet-scales.txt"),"w")
    for (x,p,axis) in jets
        r=norm(x);X=EMS.trumpet_sample(c.seed.profile,r).X
        scale=vcat(fill(X^2,18),[1.,1.],fill(X,6),ones(3))
        println(scales,join(scale,' '))
        dir=zeros(3)
        if axis<=3
            dir[axis]=r
        else
            i,j=((1,2),(1,3),(2,3))[axis-3]
            dir[i]=r/sqrt(2.);dir[j]=r/sqrt(2.)
        end
        for h in (.002,.001), k in (-2,-1,0,1,2)
            y=x+(k*h)*dir
            println(input,join(y,' '),' ',p)
            println(reference,join(packed(c,y;panel=p),' '))
        end
    end
    close(input);close(reference);close(scales)
    println("reference complete: points=$(length(points)) jets=$(length(jets)); derivative steps .002/.001 in local Cartesian units")
end
main()
