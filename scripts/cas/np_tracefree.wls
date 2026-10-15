(* Q[e11,e12,e13,e22,e23,b11,b12,b13,b22,b23,a], extension I,sqrt(2).
   Oriented orthonormal (-+++), positive NP contractions, real null rotation a.
   Witnesses: every algebraic Weyl component identity and NP binomial coefficient.
   No denominators depend on parameters; only 2,3,sqrt(2) occur. *)
ClearAll["Global`*"];
Check[
assert[q_,label_]:=If[!TrueQ[q],Print["NOT-CLOSED ",label];Exit[1]];
z[q_]:=Expand[q]===0;
e={{e11,e12,e13},{e12,e22,e23},{e13,e23,-e11-e22}};
b={{b11,b12,b13},{b12,b22,b23},{b13,b23,-b11-b22}};
eps[i_,j_,k_]:=(i-j)(j-k)(k-i)/2;
eta=DiagonalMatrix[{-1,1,1,1}];
c=ConstantArray[0,{4,4,4,4}];
Do[
 c[[1,i+1,1,j+1]]=c[[i+1,1,j+1,1]]=e[[i,j]];
 c[[i+1,1,1,j+1]]=c[[1,i+1,j+1,1]]=-e[[i,j]];
 Do[
  mixed=Sum[eps[j,k,l] b[[i,l]],{l,3}];
  c[[1,i+1,j+1,k+1]]=c[[j+1,k+1,1,i+1]]=mixed;
  c[[i+1,1,j+1,k+1]]=c[[j+1,k+1,i+1,1]]=-mixed;
  Do[c[[i+1,j+1,k+1,l+1]]=KroneckerDelta[i,k] e[[j,l]]-
   KroneckerDelta[i,l] e[[j,k]]-KroneckerDelta[j,k] e[[i,l]]+
   KroneckerDelta[j,l] e[[i,k]],{l,3}],{k,3}],{i,3},{j,3}];
Do[
 assert[z[c[[i,j,k,l]]+c[[j,i,k,l]]] && z[c[[i,j,k,l]]+c[[i,j,l,k]]] &&
  z[c[[i,j,k,l]]-c[[k,l,i,j]]] &&
  z[c[[i,j,k,l]]+c[[i,k,l,j]]+c[[i,l,j,k]]],"algebraic Weyl symmetries"],
 {i,4},{j,4},{k,4},{l,4}];
Do[assert[z[Sum[eta[[i,i]] c[[i,j,i,l]],{i,4}]],"trace-free"],{j,4},{l,4}];
q4[u_,v_,w_,x_]:=Expand[Sum[c[[i,j,k,l]] u[[i]] v[[j]] w[[k]] x[[l]],{i,4},{j,4},{k,4},{l,4}]];
ff=Table[If[i<j,f[i,j],If[i==j,0,-f[j,i]]],{i,4},{j,4}];
q2[u_,v_]:=Expand[u.ff.v];
ps[t_]:=With[{l=t[[1]],k=t[[2]],m=t[[3]],mb=t[[4]]},
 {q4[l,m,l,m],q4[l,k,l,m],q4[l,m,mb,k],q4[l,k,mb,k],q4[k,mb,k,mb]}];
ph[t_]:=With[{l=t[[1]],k=t[[2]],m=t[[3]],mb=t[[4]]},
 {q2[l,m],(q2[l,k]+q2[mb,m])/2,q2[mb,k]}];
t={{1,1,0,0},{1,-1,0,0},{0,0,1,I},{0,0,1,-I}}/Sqrt[2];
p=ps[t];fbase=ph[t];
rot[v_]:=Table[Sum[Binomial[j,k] a^k v[[j-k+1]],{k,0,j}],{j,0,Length[v]-1}];
Do[
 fixed=If[class==1,1,2];other=3-fixed;
 tr=t;tr[[other]]=t[[other]]+a(t[[3]]+t[[4]])+a^2 t[[fixed]];
 tr[[3]]=t[[3]]+a t[[fixed]];tr[[4]]=t[[4]]+a t[[fixed]];
 rp=If[class==1,rot[p],Reverse[rot[Reverse[p]]]];
 rf=If[class==1,rot[fbase],Reverse[rot[Reverse[fbase]]]];
 assert[And@@(z /@ (ps[tr]-rp)),"Weyl null rotation"];
 assert[And@@(z /@ (ph[tr]-rf)),"Maxwell null rotation"],{class,2}];
(* Generic Gauss spatial block has symmetric contracted Q, normal block
   Q-ricSpatial. The scalar trace mismatch and momentum antisymmetry are free. *)
qq=Table[q[Min[i,j],Max[i,j]],{i,3},{j,3}];
rr=Table[r[Min[i,j],Max[i,j]],{i,4},{j,4}];
raw=ConstantArray[0,{4,4,4,4}];
Do[
 raw[[1,i+1,1,j+1]]=raw[[i+1,1,j+1,1]]=qq[[i,j]]-rr[[i+1,j+1]];
 raw[[i+1,1,1,j+1]]=raw[[1,i+1,j+1,1]]=rr[[i+1,j+1]]-qq[[i,j]];
 Do[
  (* Arbitrary trace-free magnetic dual: Codazzi Bianchi, not momentum. *)
  mb=Sum[eps[j,k,l] (b[[i,l]]+Sum[eps[i,l,h] v[h],{h,3}]),{l,3}];
  raw[[1,i+1,j+1,k+1]]=raw[[j+1,k+1,1,i+1]]=mb;
  raw[[i+1,1,j+1,k+1]]=raw[[j+1,k+1,i+1,1]]=-mb;
  Do[raw[[i+1,j+1,k+1,l+1]]=KroneckerDelta[i,k] qq[[j,l]]-
   KroneckerDelta[i,l] qq[[j,k]]-KroneckerDelta[j,k] qq[[i,l]]+
   KroneckerDelta[j,l] qq[[i,k]]-Tr[qq]/2 (KroneckerDelta[i,k] KroneckerDelta[j,l]-
     KroneckerDelta[i,l] KroneckerDelta[j,k]),{l,3}],{k,3}],{i,3},{j,3}];
sr=Tr[eta.rr];
old=Table[Expand[raw[[i,j,k,l]]-(eta[[i,k]] rr[[l,j]]-eta[[i,l]] rr[[k,j]]-
 eta[[j,k]] rr[[l,i]]+eta[[j,l]] rr[[k,i]])/2+
 sr(eta[[i,k]] eta[[l,j]]-eta[[i,l]] eta[[k,j]])/6],{i,4},{j,4},{k,4},{l,4}];
eraw=Table[old[[1,i+1,1,j+1]],{i,3},{j,3}];
braw=Table[Sum[eps[j,k,l] old[[1,i+1,k+1,l+1]]/2,{k,3},{l,3}],{i,3},{j,3}];
assert[z[Tr[eraw]-(Tr[qq]-Tr[rr[[2;;4,2;;4]]]-rr[[1,1]])],"electric trace Hamiltonian"];
Do[assert[z[Sum[eps[i,j,k] braw[[j,k]],{j,3},{k,3}]-
 Sum[eta[[h,h]] raw[[h,1,h,i+1]],{h,4}]+rr[[1,i+1]]],"magnetic axial momentum"],{i,3}];
etf=(eraw+Transpose[eraw])/2-Tr[eraw] IdentityMatrix[3]/3;
btf=(braw+Transpose[braw])/2-Tr[braw] IdentityMatrix[3]/3;
(* Replace only the ten independent STF symbols, avoiding duplicate rules. *)
new=c/.{e11->etf[[1,1]],e12->etf[[1,2]],e13->etf[[1,3]],e22->etf[[2,2]],e23->etf[[2,3]],
 b11->btf[[1,1]],b12->btf[[1,2]],b13->btf[[1,3]],b22->btf[[2,2]],b23->btf[[2,3]]};
oldp=Block[{c=old},ps[t]];newp=Block[{c=new},ps[t]];
assert[z[oldp[[1]]-newp[[1]]] && z[oldp[[5]]-newp[[5]]],"radiative preservation"];
Print["PROVED: STF E/B reconstruction has all algebraic Weyl symmetries and zero trace"];
Print["PROVED: both null-rotation binomial laws from actual tensor and Faraday contractions"];
Print["PROVED: Psi0 and Psi4 preserved; removed electric trace=H, magnetic axial=-M (geometric constraints)"];
Print["Domain: oriented orthonormal (-+++), real arbitrary tensor components and a; no parameter poles"];
Exit[0],Print["NOT-CLOSED Wolfram message"];Exit[1]];
