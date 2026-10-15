# Optional diagnostic of the unchanged reference taper near |z|=6.
using EMS, ForwardDiff
function main()
    for t in (Float64,BigFloat)
        setprecision(512) do
            z=t(-6.000000000000002)
            f(x)=last(EMS.pulse_gaussian(x))
            println(t," value=",f(z)," D1=",ForwardDiff.derivative(f,z),
                " D2=",ForwardDiff.derivative(x->ForwardDiff.derivative(f,x),z))
        end
    end
end
main()
