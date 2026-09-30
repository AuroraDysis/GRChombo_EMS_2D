# T9 characteristic and sampling evidence

The claim is conditional: it concerns the frozen longitudinal shift block,
not the complete CCZ4 characteristic system. Set `p = partial_r beta_n`,
`mu > 0`, `g = n_i h^ij n_j > 0`, `c = sqrt(4 mu g / 3) > 0`, and
`b = beta_n` from the current front. For `(p,Gamma_n,B_driver_n)`, the code gives

```
P = [[b, mu, -1], [4g/3, b, 0], [0, 0, 0]],   U_t = P U_r + forcing.
```

The first row follows from ExperimentalGauge.hpp: beta_t = beta advection +
mu Gamma - eta beta - B_driver. The second follows from the normal projection
of CCZ4Cartoon.impl.hpp: `h^jk beta^i_,jk + h^ij beta^k_,jk/3`.
The third follows from `B_driver_t = -0.1 B_driver`, with no advection.
Metric/lapse/K/Theta principal perturbations, transverse derivatives, radial
cartoon terms, damping, matter and KO are excluded from this block. In particular
`-4 alpha h^ij K_,j/3 + 2 alpha h^ij Theta_,j` remain forcing, not zero terms in
the actual evolution. This is an approximate projection of saved states; no
gauge equation is altered. Its near-puncture accuracy is not asserted.

The algebra is `QQ[b,c,mu]`, with `g = 3c^2/(4mu)`. The exact witness is
`l P - lambda l = 0` for each left eigenvector, after rational cancellation.
The domain excludes `mu=0`, `c=0`, and `b=+/-c`; there are no hidden radical
branch changes. The `c>0` branch is declared, and current positive metrics
are checked. At `lambda=b-c`, the outgoing coordinate speed is `c-b`:

```
W_out = Gamma_n - (c/mu) p - c/[mu(c-b)] B_driver_n.
W_in  = Gamma_n + (c/mu) p - c/[mu(c+b)] B_driver_n.
```

[t9-characteristic-verify.py](t9-characteristic-verify.py) first probes a
nondegenerate special value, then closes both rational numerator witnesses
with exact assertions. The result is **PROVED within this conditional block**.
An independent mpmath matrix/vector computation at 32 points, seed 9009,
60 digits, c in [0.6,1.4], mu in [0.5,1], b in [-0.2,0.2], gives worst
residual 3.111508e-61 (**CORROBORATED numerically**). This is not a proof that
the saved front is a pure mode of the full CCZ4 system.

The separate finite sampling witness is `sum_i L_i(x) i^k = x^k` in QQ[x],
for n=6,8, every 0<=k<n, and derivatives 0..3: 56 exact polynomial
coefficient checks, **PROVED**. The handwritten runtime product-rule kernel
is checked against independently evaluated polynomial derivatives and the
existing T5 weight function. Worst scaled Float64 residual is 6.063298e-13.
No symbolic result is emitted as a generated arithmetic kernel.

Eight further exact QQ[nx,ny,r] derivative witnesses check radial cubic and
quintic vector projection onto a ray. On nx^2+ny^2=1 these yield beta_n=r^3
and Gamma_n=r^5. An end-to-end runtime test samples these vectors and the
linear driver on all three rays, with parity and directional derivatives
0..3, at both stencil sizes. Its worst scaled residual is 2.185266e-13.

The native Gamma-metric implementation reuses the sealed exp-0020 reducer's
fourth-order metric derivative helper and retains the separate 2D and ww
terms. Against the reducer's independently assembled Z, `C_Gamma=2Z/chi`
agrees on all three rung test blocks to scaled residual 3.469447e-18.
This verifies transcription, not spatial accuracy at an unresolved front.

For jump integration the dense P6/P8 values are joined with a C2 cubic spline.
Its piecewise linear second derivative is integrated by exact signed trapezoid
pieces split at each zero. The endpoint first-derivative difference is an
independent check on the sum. The cubic witness x^3 over [-1,1] gives signed
lobes +3/-3; every recorded data integral checks the corresponding endpoint
identity. Half-density reconstruction and the direct local-polynomial endpoint
derivative are recorded separately. These are numerical reconstruction
sensitivities, not certified continuum error bounds.
