using EMS, ForwardDiff, Printf

# Independent static-chart curvature and electrostatic first integral, using
# the supplied spectral metric (not the Cartesian Gauss/Codazzi reconstruction).
function reference_point(p::EMS.TrumpetProfile, R::Float64)
    s=EMS.trumpet_s(p,R)
    C=EMS.trumpet_meta(p,"C"); nu=EMS.trumpet_meta(p,"nu")
    q=EMS.trumpet_meta(p,"q_native")
    radius(z::Real)=EMS.trumpet_jet(p,1,z)/(1-z)
    delta(z::Real)=(1-z)^2*EMS.trumpet_jet(p,2,z)
    phi(z::Real)=EMS.trumpet_meta(p,"phi_inf")+(1-z)*EMS.trumpet_jet(p,3,z)
    function A(z::Real)
        Y=EMS.trumpet_jet(p,1,z); Y1=EMS.trumpet_jet(p,1,z,1)
        alpha=exp(-delta(z))*nu*z*(Y+(1-z)*Y1)/(Y*(1+(nu-1)*z))
        alpha^2-C^2/radius(z)^4
    end
    N(z::Real)=exp(2delta(z))*A(z)
    derivative(f,z::Real)=ForwardDiff.derivative(f,z)
    r=radius(s); rs=derivative(radius,s); rss=derivative(z->derivative(radius,z),s)
    function jet(f)
        fs=derivative(f,s); fss=derivative(z->derivative(f,z),s)
        (f(s),fs/rs,(fss*rs-fs*rss)/rs^3)
    end
    n,nr,nrr=jet(N); d,dr,drr=jet(delta)
    coupling=EMS.f(p.coupling,phi(s))
    Vr=q*exp(-d)/(coupling*r*r)
    psi2=(nrr/2-n*drr-3nr*dr/2+n*dr^2-nr/r+n*dr/r-(1-n)/r^2)/6
    phi1=exp(d)*Vr/2
    (;R,r,N=n,delta=d,Vr,psi2,phi1)
end

function main()
    length(ARGS)>=3 || error("usage: static_oracle.jl profile output R...")
    p=read_trumpet(ARGS[1])
    println("route=np-static-chart-oracle profile=",ARGS[1]," coefficient_type=",eltype(p.coeff))
    open(ARGS[2],"w") do io
        println(io,"r,R,N,delta,Vr,Psi2,Phi1")
        for text in ARGS[3:end]
            v=reference_point(p,parse(Float64,text))
            v.N>0 || error("oracle sample must be outside the static horizon")
            @printf(io,"%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",v.R,v.r,v.N,v.delta,v.Vr,v.psi2,v.phi1)
            if p.meta["Q_s"]=="0" || iszero(parse(Float64,p.meta["Q_s"]))
                Q=parse(Float64,p.meta["Q"])
                expected=-1/v.r^3+Q^2/v.r^4
                @assert abs(v.psi2-expected)<1e-10*max(abs(expected),1/v.r^3)
            end
        end
    end
end
main()
