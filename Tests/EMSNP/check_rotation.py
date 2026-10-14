"""Direct tensor contractions of a frozen-checkpoint jet; standard library only.

This is an acceptance probe, not the requested production NP extractor.
Run: python3 check_rotation.py results/matter/np_jet.txt
"""
import itertools
import json
import math
from pathlib import Path
import sys

INDICES = range(4)
PAIRS = list(itertools.product(INDICES, repeat=2))
QUADS = list(itertools.product(INDICES, repeat=4))


def load(path):
    data = {name: {} for name in ("metric", "ricci", "faraday", "frame", "riemann")}
    for line in Path(path).read_text().splitlines():
        name, *values = line.split()
        if name in data:
            data[name][tuple(map(int, values[:-1]))] = float(values[-1])
        else:
            data[name] = list(map(float, values))
    for name in ("metric", "ricci", "faraday", "frame"):
        assert set(data[name]) == set(PAIRS), name
    assert set(data["riemann"]) == set(QUADS)
    assert all(math.isfinite(x) for name in ("metric", "ricci", "faraday", "frame", "riemann")
               for x in data[name].values())
    return data


def add(*vectors):
    return [sum(v[i] for v in vectors) for i in INDICES]


def scale(a, v):
    return [a*x for x in v]


def contract(tensor, *vectors):
    total = 0j
    for indices, value in tensor.items():
        for vector, i in zip(vectors, indices):
            value *= vector[i]
        total += value
    return total


def absolute_contraction(tensor, *vectors):
    return sum(abs(value)*math.prod(abs(v[i]) for v,i in zip(vectors,indices))
               for indices,value in tensor.items())


def change_roundoff(C,F,t):
    l,k,m,mb=t
    cp=[absolute_contraction(C,*v) for v in ((l,m,l,m),(l,k,l,m),(l,m,mb,k),
                                            (l,k,mb,k),(k,mb,k,mb))]
    fp=[absolute_contraction(F,l,m),
        (absolute_contraction(F,l,k)+absolute_contraction(F,mb,m))/2,
        absolute_contraction(F,mb,k)]
    # 4096 eps covers the <=256-term contractions, products, and two evaluations.
    return [4096*sys.float_info.epsilon*x for x in
            (2*fp[1]*cp[3]+3*cp[2]*fp[2],2*fp[1]*cp[1]+3*cp[2]*fp[0],fp[2])]


def scalars(C, F, tetrad):
    l, k, m, mb = tetrad
    p = [contract(C, *v) for v in ((l,m,l,m), (l,k,l,m), (l,m,mb,k),
                                  (l,k,mb,k), (k,mb,k,mb))]
    f = [contract(F,l,m), (contract(F,l,k)+contract(F,mb,m))/2, contract(F,mb,k)]
    return p, f


def observables(C, F, t):
    p, f = scalars(C,F,t)
    return [2*f[1]*p[3]-3*p[2]*f[2], 2*f[1]*p[1]-3*p[2]*f[0], f[2]]


def rotate(t, a, fixed):
    l, k, m, mb = t
    if fixed == "l":
        return [l,add(k,scale(a,m),scale(a,mb),scale(a*a,l)),
                add(m,scale(a,l)),add(mb,scale(a,l))]
    return [add(l,scale(a,m),scale(a,mb),scale(a*a,k)),k,
            add(m,scale(a,k)),add(mb,scale(a,k))]


def encode(z):
    return [z.real,z.imag]


def self_check():
    # Independent electric/magnetic decomposition makes a genuine trace-free Weyl
    # tensor, so the binomial laws are checked against tensor rotations, not themselves.
    def eps(i,j,k):
        return 0 if len({i,j,k})<3 else (1 if (i,j,k) in ((0,1,2),(1,2,0),(2,0,1)) else -1)
    E=[[.7,.02,.03],[.02,-.4,.01],[.03,.01,-.3]]
    H=[[.2,.01,.02],[.01,.1,.03],[.02,.03,-.3]]
    C=dict.fromkeys(QUADS,0.)
    def set_component(a,b,c,d,x):
        for aa,bb,s in ((a,b,1),(b,a,-1)):
            for cc,dd,q in ((c,d,1),(d,c,-1)):
                C[aa,bb,cc,dd]=C[cc,dd,aa,bb]=s*q*x
    for i in range(3):
        for j in range(3):
            set_component(0,i+1,0,j+1,E[i][j])
            for k in range(3):
                set_component(0,i+1,j+1,k+1,sum(eps(j,k,h)*H[i][h] for h in range(3)))
                for l in range(3):
                    set_component(i+1,j+1,k+1,l+1,-sum(eps(i,j,h)*eps(k,l,q)*E[h][q]
                                                                    for h in range(3) for q in range(3)))
    F=dict.fromkeys(PAIRS,0.)
    for i in range(3):
        F[i+1,0]=(.3,.01,.02)[i]; F[0,i+1]=-F[i+1,0]
        for j in range(3): F[i+1,j+1]=sum(eps(i,j,k)*(.1,.02,.01)[k] for k in range(3))
    t=[[1/math.sqrt(2),1/math.sqrt(2),0,0],[1/math.sqrt(2),-1/math.sqrt(2),0,0],
       [0,0,1/math.sqrt(2),1j/math.sqrt(2)],[0,0,1/math.sqrt(2),-1j/math.sqrt(2)]]
    base=scalars(C,F,t)
    assert max(abs(sum((-1 if a==0 else 1)*C[a,b,a,d] for a in INDICES))
               for b,d in PAIRS) < 64*sys.float_info.epsilon
    worst=0
    for fixed in ("l","k"):
        for a in (.01,-.03):
            direct=scalars(C,F,rotate(t,a,fixed))
            for values,got in zip(base,direct):
                original=values if fixed=="l" else values[::-1]
                expected=[sum(math.comb(j,k)*a**k*original[j-k] for k in range(j+1)) for j in range(len(values))]
                if fixed=="k": expected=expected[::-1]
                worst=max(worst,max(abs(x-y) for x,y in zip(got,expected)))
    assert worst<4096*sys.float_info.epsilon, worst
    return worst


def check(path):
    d = load(path)
    g, ric = d["metric"], d["ricci"]
    n,u,v,w = [[d["frame"][i,j] for j in INDICES] for i in INDICES]
    t = [scale(1/math.sqrt(2),add(n,u)),scale(1/math.sqrt(2),add(n,scale(-1,u))),
         scale(1/math.sqrt(2),add(v,scale(1j,w))),scale(1/math.sqrt(2),add(v,scale(-1j,w)))]
    # Trace in the exported orthonormal frame avoids a separate inverse implementation.
    trace = -contract(ric,n,n)+sum(contract(ric,x,x) for x in (u,v,w))
    C = {(a,b,c,e): r-(g[a,c]*ric[e,b]-g[a,e]*ric[c,b]-g[b,c]*ric[e,a]+g[b,e]*ric[c,a])/2
         +trace*(g[a,c]*g[e,b]-g[a,e]*g[c,b])/6
         for (a,b,c,e),r in d["riemann"].items()}
    F = d["faraday"]
    p,f = scalars(C,F,t)
    baseline = observables(C,F,t)
    target = [[0,-1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,1,0]]
    norm_error = max(abs(contract(g,t[i],t[j])-target[i][j]) for i,j in PAIRS)
    assert norm_error < 64*sys.float_info.epsilon, norm_error
    # Magnitude-derived roundoff budget for the legacy normalization identity.
    cnorm = max(abs(x) for x in C.values())
    frame_scale = max(sum(abs(x) for x in q) for q in t)**4
    tolerance = 4096*sys.float_info.epsilon*cnorm*frame_scale
    legacy_error = abs(p[4]-d["legacy"][0])
    assert legacy_error <= tolerance, (legacy_error,tolerance)
    # Residual trace: the field-equation Ricci and finite-difference Gauss pieces
    # need not satisfy the numerical constraints exactly. Report, never project it away.
    trace_error = max(abs(sum((-1 if a==0 else 1)*contract(C,q,x,q,y)
                             for a,q in enumerate((n,u,v,w)))) for x in (n,u,v,w) for y in (n,u,v,w))
    output = {"path": str(path), "point_x_y_dx": d["point"], "psi": list(map(encode,p)),
              "phi": list(map(encode,f)), "X_out_X_in_Phi2":list(map(encode,baseline)),
              "tetrad_error":norm_error,"weyl_trace_residual":trace_error,
              "legacy_difference":legacy_error,"legacy_roundoff_tolerance":tolerance,"rotations":{}}
    for fixed in ("l","k"):
        rows=[]
        for a in (1e-2,1e-3,1e-4,1e-5,1e-6):
            tp,tm = rotate(t,a,fixed),rotate(t,-a,fixed)
            for tr in (tp,tm):
                assert max(abs(contract(g,tr[i],tr[j])-target[i][j]) for i,j in PAIRS) < 128*sys.float_info.epsilon
            plus,minus = observables(C,F,tp),observables(C,F,tm)
            delta = [x-y for x,y in zip(plus,baseline)]
            derivative = [(x-y)/(2*a) for x,y in zip(plus,minus)]
            rows.append({"a":a,"absolute_change":list(map(abs,delta)),
                         "roundoff_budget":[x+y for x,y in zip(change_roundoff(C,F,t),change_roundoff(C,F,tp))],
                         "central_derivative":list(map(encode,derivative))})
        output["rotations"][fixed]=rows
    output["decade_slopes"]={fixed:[
        [math.log(rows[i]["absolute_change"][j]/rows[i+1]["absolute_change"][j],10)
         for j in range(3)] for i in range(len(rows)-1)]
        for fixed,rows in output["rotations"].items()}
    # Expected coefficients for a trace-free Weyl tensor, for diagnostic comparison.
    output["ideal_linear_coefficients"]={
        "l":list(map(encode,(2*f[0]*p[3]-6*f[2]*p[1],2*f[1]*p[0]-4*f[0]*p[1],2*f[1]))),
        "k":list(map(encode,(2*f[1]*p[4]-4*f[2]*p[3],2*f[2]*p[1]-6*f[0]*p[3],0j)))}
    return output


if __name__ == "__main__":
    print(json.dumps({"self_check_residual":self_check(),
                      "jets":[check(p) for p in sys.argv[1:]]}, indent=2, allow_nan=False))
