(* Q[a,p0,p1,p2,p3,p4,f0,f1,f2], characteristic zero. No poles,
   branches or data-dependent assumptions. Witness: coefficient of a^1.
   The binomial laws follow by multilinearity of symmetric rank-4 Weyl
   and rank-2 Maxwell spinors under iota -> iota+a o (and the reverse).
   Every Expand is a deterministic polynomial CLOSE, no simplifier. *)
ClearAll["Global`*"];
Check[
  p={p0,p1,p2,p3,p4}; f={f0,f1,f2};
  xo[ps_,fs_]:=2 fs[[2]] ps[[4]]-3 ps[[3]] fs[[3]];
  xi[ps_,fs_]:=2 fs[[2]] ps[[2]]-3 ps[[3]] fs[[1]];
  rotL[v_]:=Table[Sum[Binomial[j,k] a^k v[[j-k+1]],{k,0,j}],{j,0,Length[v]-1}];
  rotK[v_]:=Reverse[rotL[Reverse[v]]];
  dl=Expand /@ {xo[rotL[p],rotL[f]]-xo[p,f],xi[rotL[p],rotL[f]]-xi[p,f]};
  dk=Expand /@ {xo[rotK[p],rotK[f]]-xo[p,f],xi[rotK[p],rotK[f]]-xi[p,f]};
  actual=Coefficient[#,a,1]& /@ Join[dl,dk];
  expected={2 f0 p3-6 f2 p1,2 f1 p0-4 f0 p1,
            2 f1 p4-4 f2 p3,2 f2 p1-6 f0 p3};
  If[!And@@Thread[(Expand /@ (actual-expected))==0],Print["NOT-CLOSED coefficient identity"];Exit[1]];
  Print["via: CLI-fallback; version: ",$Version];
  Print["Domain: polynomial ring over Q; a real rotation suffices to refute; no excluded loci"];
  Print["Linear coefficients (l-fixed out,in; k-fixed out,in): ",ToString[actual,InputForm]];
  Print["PROVED: the four exact coefficient identities"];
  background={p0->0,p1->0,p3->0,p4->0,f0->0,f2->0};
  If[!And@@Thread[(actual/.background)==0],Exit[1]];
  Print["PROVED: all four coefficients vanish on the aligned Coulomb background"];
  (* Real, arbitrarily small non-Coulomb perturbation; not checkpoint data. *)
  counter={p0->0,p1->0,p2->-1,p3->eps,p4->0,f0->eps,f1->1,f2->0};
  witness=Expand[actual[[1]]/.counter];
  If[Expand[witness-2 eps^2]=!=0,Exit[1]];
  Print["Exact counterexample, X_out about l: coefficient = ",ToString[witness,InputForm]," (eps != 0)"];
  Print["REFUTED: universal O(a^2) response at fixed perturbed fields"];
  Print["CAS scope: algebraic identity — production physics not certified"];
  Exit[0], Print["NOT-CLOSED: Wolfram message"];Exit[1]
];
