# Run from EMS.jl with --project=test; no writes in EMS.jl.
using EMS, Printf

function postprocess(input::String, output::String)::Nothing
    lines=readlines(input); names=split(first(lines),',')
    pins=Dict("rn"=>"b98d6dc915b60656e848a6ca45d38bfb2f7d5ba2a5b6273395ab3ee3fe412c05",
              "scalarized-positive-n0"=>"dd4013ca19a746b5bab847fcaaf7e9ece4d5ddf7f4ba5b29d84cb86386b23a07")
    open(output,"w") do io
        println(io,first(lines),",R_A,M_irr,e,M_eq,delta_M,delta_M_measurement,delta_M_table,delta_A,delta_Q,model_phi_H,model_Qs_over_R_A,expansion_residual,mass_error,charge_relative_error,area_relative_error,phi_difference,mass_gate,charge_gate,budget_covers,area_bias_covers,RN_relation_error,table_hash,equilibrium_status,M_settled,spread_M_eq")
        for line in lines[2:end]
            row=Dict(zip(names,split(line,',')))
            val(k)=parse(Float64,row[k]); branch=row["branch"]
            coupling=EMSCoupling(val("ems_alpha"),val("f0"),val("f1"),val("f2"))
            stem=branch=="rn" ? "rn" : "alpha20-positive"
            # Bind every run's physical parameters, not metadata inferred from the table.
            f=read_horizon_family(joinpath(pwd(),"artifacts/horizon-family/production",stem*".hfamily");
                coupling,phi_inf=val("phi_inf"),branch,expected_sha256=pins[branch])
            A=val("A"); Q=val("Q"); n=parse(Int,row["N_theta"])
            b=(pi/(2n))/sin(pi/(2n))-1
            dA=abs(A)*b; dQ=abs(Q)*b
            m=horizon_mass(f,A,Q;delta_A=dA,delta_Q=dQ); s=horizon_sample(f,m.e)
            expected_Q=0.700784945779
            expected_A=branch=="rn" ? 4pi*(1+sqrt(1-expected_Q^2))^2 : 38.510457444941309
            error=abs(m.M_eq-1); qerror=abs(Q/expected_Q-1); aerror=abs(A/expected_A-1)
            rn=m.R_A*(1+m.e^2)/2
            duplicate_masses=[horizon_mass(f,parse.(Float64,split(pair,':'))...).M_eq
                              for pair in split(row["duplicate_A_Q"],';')]
            values=(m.R_A,m.R_A/2,m.e,m.M_eq,m.delta_M,m.delta_M_measurement,m.delta_M_table,
                    dA,dQ,s.phi_H,s.Q_s_over_R_A,m.R_A*sqrt(val("expansion_squared")),error,qerror,
                    aerror,val("phi_mean")-s.phi_H,error <= (n==48 ? 1e-3 : 3e-4),qerror<=5e-4,
                    error<=m.delta_M,aerror<=b,branch=="rn" ? abs(m.M_eq-rn) : NaN,f.sha256,
                    "BRANCH_UNVERIFIED_STATIC_PROJECTION","unset",maximum(duplicate_masses)-minimum(duplicate_masses))
            println(io,line,',',join(values,','))
        end
    end
    println("Julia horizon_mass post-processing complete: ",output)
    nothing
end
postprocess(abspath(ARGS[1]),abspath(ARGS[2]))
