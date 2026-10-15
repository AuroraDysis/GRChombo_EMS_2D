# Reproduce the archived exact-AD diagnostic using unchanged delivered fields.
# BigFloat is needed to preserve exp(-1/t) near taper faces; double-double
# shares Float64's exponent range and cannot represent these intermediates.
using EMS, ForwardDiff, LinearAlgebra, DelimitedFiles
BLAS.set_num_threads(1)
function main()
    profile,companion,points,output=ARGS
    c=EMS.read_ctt_pulse(companion,profile)
    packed(x,p)=begin
        a=EMS.ctt_pulse_fields(c,x;panel=p)
        vcat(vec(a.gamma),vec(a.K),a.phi,a.Pi,a.E,a.B,a.shift)
    end
    rows=readdlm(points)
    open(output,"w") do io
        for row in eachrow(rows)
            x=row[1:3];p=Int(row[4]);axis=Int(row[5]);r=norm(x)
            X=EMS.trumpet_sample(c.seed.profile,r).X
            scale=vcat(fill(X^2,18),[1.,1.],fill(X,6),ones(3))
            dir=zeros(3);dir[axis]=r
            v,d,dd=setprecision(512) do
                xb=BigFloat.(x);db=BigFloat.(dir);sb=BigFloat.(scale)
                f(u)=packed(xb+u[1]*db,p).*sb
                u=BigFloat[0]
                (Float64.(f(u)),Float64.(vec(ForwardDiff.jacobian(f,u))),
                    Float64.(vec(ForwardDiff.jacobian(u->vec(ForwardDiff.jacobian(f,u)),u))))
            end
            println(io,join(vcat(scale,v,d,dd),' '))
        end
    end
    println("AD reference: 512-bit arithmetic, $(size(rows,1)) directional samples")
end
main()
