using EMS, LinearAlgebra
BLAS.set_num_threads(1)
function main()
    profile,data,out=ARGS
    for name in ("zero","n24")
        c=EMS.read_ctt_pulse(joinpath(data,name*"-final.ctt"),profile)
        function fields(x)
            a=EMS.ctt_pulse_fields(c,x)
            g=a.gamma;K=sum(inv(g).*a.K);chi=det(g)^(-1/3)
            h=chi*g;A=chi*(a.K-K*g/3)
            vcat(chi,h[1,1],h[1,2],h[2,2],h[3,3],K,A[1,1],A[1,2],A[2,2],A[3,3],a.phi,a.Pi,a.E,a.B,sqrt(chi),a.shift[1:2])
        end
        points=open(joinpath(out,"ccz4-"*name*"-points.txt"),"w")
        values=open(joinpath(out,"ccz4-"*name*"-julia.txt"),"w")
        scales=open(joinpath(out,"ccz4-"*name*"-scales.txt"),"w")
        radii=[1e-18,1e-10,.0071,.02,1.,20.,39.,55.,128.,256.]
        append!(radii,[p.hi for p in c.panels[1:end-1]])
        for r in radii, direction in ([1.,0.,0.],[0.,1.,0.],[inv(sqrt(2)),inv(sqrt(2)),0.])
            x=r*[sqrt(1-1e-12),1e-6,0.]
            X=EMS.trumpet_sample(c.seed.profile,r).X
            scale=ones(21);scale[1]=1/X^2;scale[13:18].=X;scale[19]=1/X
            println(scales,join(scale,' '))
            for k in -2:2
                y=x+.001*k*r*direction
                println(points,join(y,' ')," 0")
                println(values,join(fields(y),' '))
            end
        end
        close(points);close(values);close(scales)
    end
    println("CCZ4 value and consumed derivative reference complete")
end
main()
