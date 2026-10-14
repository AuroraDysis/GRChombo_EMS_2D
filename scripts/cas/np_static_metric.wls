(* Exact witness over Q(r,A,N,A',N',A'',N''). Static chart r>0,
   A>0,N>0,0<theta<Pi. The final N,delta formula extends analytically.
   CLOSE is Cancel[Together[residual]]===0; no general simplifier.
   No radical is used in the curvature witness. *)
ClearAll["Global`*"];
Check[
coords={t,r,th,ph}; metric=DiagonalMatrix[{-A[r],1/NN[r],r^2,r^2 Sin[th]^2}];
inverse=DiagonalMatrix[{-1/A[r],NN[r],1/r^2,1/(r^2 Sin[th]^2)}];
normal[x_]:=Cancel[Together[x]];
assert[claim_,label_]:=If[!TrueQ[claim],Print["NOT-CLOSED ",label];Exit[1]];
gamma=Table[normal[Sum[inverse[[a,e]] (D[metric[[e,b]],coords[[c]]]+
    D[metric[[e,c]],coords[[b]]]-D[metric[[b,c]],coords[[e]]])/2,{e,4}]],{a,4},{b,4},{c,4}];
riemann[a_,b_,c_,d_]:=normal[Sum[metric[[a,e]] (D[gamma[[e,d,b]],coords[[c]]]-
    D[gamma[[e,c,b]],coords[[d]]]+Sum[gamma[[e,c,f]] gamma[[f,d,b]]-
    gamma[[e,d,f]] gamma[[f,c,b]],{f,4}]),{e,4}]];
components={riemann[1,2,1,2] NN[r]/A[r],riemann[1,3,1,3]/(A[r] r^2),
    riemann[2,3,2,3] NN[r]/r^2,riemann[3,4,3,4]/(r^4 Sin[th]^2)};
expected={NN[r] (A''[r]/(2 A[r])-A'[r]^2/(4 A[r]^2))+NN'[r] A'[r]/(4 A[r]),
    NN[r] A'[r]/(2 A[r] r),-NN'[r]/(2 r),(1-NN[r])/r^2};
assert[And@@(TrueQ[#===0]& /@ (normal /@ (components-expected))),"four metric Riemann components"];
(* In an orthonormal frame x=R0101,y=R0202,z=R1212,w=R2323,
   Ric00=x+2y,Ric11=-x+2z,scalar=-2x-4y+4z+2w. *)
psiFromRicci=(x+(-x+2z)/2-(x+2y)/2-(-2x-4y+4z+2w)/6)/2;
assert[normal[psiFromRicci-(x-y+z-w)/6]===0,"Ricci subtraction"];
psiA=normal[(components[[1]]-components[[2]]+components[[3]]-components[[4]])/6];
psiND=(NN''[r]/2-NN[r] delta''[r]-3 NN'[r] delta'[r]/2+NN[r] delta'[r]^2-
    NN'[r]/r+NN[r] delta'[r]/r-(1-NN[r])/r^2)/6;
rules={A[r]->NN[r] Exp[-2 delta[r]],A'[r]->D[NN[r] Exp[-2 delta[r]],r],
    A''[r]->D[NN[r] Exp[-2 delta[r]],{r,2}]};
assert[normal[(psiA/.rules)-psiND]===0,"A to N delta"];
rn=1-2 mass/r+charge^2/r^2;
rnRules={NN[r]->rn,NN'[r]->D[rn,r],NN''[r]->D[rn,{r,2}],delta'[r]->0,delta''[r]->0};
assert[normal[(psiND/.rnRules)-(-mass/r^3+charge^2/r^4)]===0,"RN"];
assert[normal[((psiND/.rnRules)/.charge->0)+mass/r^3]===0,"Schwarzschild"];
(* Native Maxwell first integral: V'=Q exp(-delta)/(sqrt(8 Pi) F r^2).
   E_u=exp(delta) V', Phi1=E_u/2, positive radial E for Q>0. *)
phi1=charge/(2 Sqrt[8 Pi] coupling r^2);
assert[normal[Exp[dd] (charge Exp[-dd]/(Sqrt[8 Pi] coupling r^2))/2-phi1]===0,"native Phi1"];
Print["PROVED four static-metric Riemann projections, Ricci subtraction, Psi2(N,delta), Schwarzschild, RN and native Phi1"];
Print["Psi2 = ",ToString[psiND,InputForm]];
Print["Phi1 = ",ToString[phi1,InputForm]];
Print["Domain r>0,A>0,N>0,0<theta<Pi,F>0; poles excluded r=0,A=0,N=0,sin(theta)=0,F=0"];
Print["CAS scope: algebraic identity — production physics not certified"];
Exit[0],Print["NOT-CLOSED Wolfram message"];Exit[1]];
