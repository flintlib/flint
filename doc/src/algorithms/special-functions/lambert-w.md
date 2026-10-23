---
title: The Lambert W function
summary: Evaluation of all branches of the Lambert W function in real and complex ball arithmetic — starting values from asymptotic expansions, the series at the branch point and double-precision approximations, Halley iteration, and an a posteriori error bound from a bound for W' that also certifies the branch.
complexity: O(log p) Halley steps with doubling precision, each dominated by one exponential, so the cost is a small multiple of that of exp at precision p.
order: 95
---

The Lambert $W$ function solves $W(z)e^{W(z)} = z$. It has infinitely
many branches $W_k$, $k \in \mathbb{Z}$; the principal branch $W_0$ is
real on $[-1/e, \infty)$ and $W_{-1}$ is real on $[-1/e, 0)$. The
standard branch cuts are those of Corless et al. [@CGHJK1996]: $W_0$
has the cut $(-\infty, -1/e]$, the other branches have the cut
$(-\infty, 0]$, and $W_{-1}$ has in addition a cut along the real
segment $(-1/e, 0)$, where it meets $W_1$. `acb_lambertw` also offers
two nonstandard conventions through its flags: `ACB_LAMBERTW_LEFT`
rotates the cut on $(-\infty, 0)$ to $(0, +\infty)$, and
`ACB_LAMBERTW_MIDDLE` joins $W_{-1}$ in the upper half-plane and $W_1$
in the lower half-plane through $(-1/e, 0)$, which extends the real
function $W_{-1}(x)$ analytically. FLINT implements the function in
[`arb`](arb) (`arb_lambertw`, branches 0 and $-1$ for real arguments)
and [`acb`](acb) (`acb_lambertw`, any branch $k$ given as an `fmpz`),
following
Johansson [@Joh2017b]. Unlike the functions of the rest of this group,
$W$ is not given by a convenient series; it is computed by solving the
equation.

## Starting values {#initial}

- **Large $|z|$ or $k \ne 0$.** The asymptotic expansion
  $$
  W_k(z) \sim L_1 - L_2 + \frac{L_2}{L_1} + \frac{L_2(L_2 - 2)}{2L_1^2}
  + \cdots, \qquad L_1 = \log z + 2\pi i k, \quad L_2 = \log L_1,
  $$
  gives a starting value (in double precision when the exponents fit).
  `acb_lambertw_asymp` evaluates the full bivariate series in $1/L_1$
  and $L_2/L_1$ with $L$ and $M$ terms and a rigorous error bound,
  which is finite when $|\log z|$ is large; for huge $|z|$ or $|k|$
  this gives the final result directly.
- **Near the branch point $z = -1/e$.** With $t = 2(ez + 1)$, $W(z) =
  -1 + \sqrt{t} - t/3 + \frac{11}{72}t^{3/2} - \cdots$; the series is
  evaluated with a bound on its coefficients, and when $z$ is extremely
  close to $-1/e$ it gives the result directly. The quantity $ez + 1$
  is computed once at high precision, since it suffers from
  cancellation exactly where it matters.
- **Moderate $z$.** Double-precision approximations: `d_lambertw` for
  the principal branch, and rational approximations followed by a
  Halley step (`d_lambertw_branch1`) for $W_{-1}$ on $(-1/e, 0)$.
- **Small $|z|$ on the principal branch.** The Taylor series $W_0(z) =
  z - z^2 + \frac32 z^3 - \cdots$, whose coefficients are bounded by
  $4^k/16$ for $k \ge 2$.

## Halley iteration {#halley}

From an approximation $w$, Halley's method for $f(w) = we^w - z$ gives
$$
w' = w - \frac{we^w - z}{e^w(w+1) - \frac{(w+2)(we^w - z)}{2w+2}},
$$
which converges cubically. The iteration runs on midpoints (floating
point, no ball arithmetic) with a precision that roughly triples at
each step, so that the total cost is a small multiple of the last step,
which is dominated by one exponential at the target precision.

## Certification {#certify}

The result is certified a posteriori: for the final approximation $w$,
the residual $t = we^w$ is computed in ball arithmetic, and the error is
bounded by
$$
|W(z) - w| = |W(z) - W(t)| \le |z - t| \cdot \max_{u \in [t, z]} |W'(u)|,
$$
where the maximum is over the straight segment from $t$ to $z$, which
must not cross a branch cut of $W_k$. For real arguments $W'$ is
bounded with $W_0'(x) \le 1/(x+1)$ for $x \ge 0$, $|W_0'(x)| <
2/\sqrt{1 + ex}$ for $x > -1/e$ and $|W_{-1}'(x)| \le 2/\sqrt{1+ex} +
2/|x|$; in the complex case, `acb_lambertw_bound_deriv` uses bounds of
the form $|W'(z)| \le 2.25/\sqrt{|t|(1+|t|)}$ with $t = ez + 1$ near the
branch point and analogous ones elsewhere. A separate check
(`acb_lambertw_check_branch`) verifies that $w$ lies in the region of
the requested branch $k$, which the residual alone cannot detect.

For an inexact input ball, the function is evaluated at the midpoint
and the propagated error is bounded with the same derivative bounds.
If the input ball straddles a branch cut, the result is the union of
the values on both sides, computed separately.

## Use in FLINT {#use}

- `arb_lambertw`, `acb_lambertw` and the generic `gr_lambertw` and
  `gr_lambertw_fmpz` for ball contexts.
- The power series versions `arb_poly_lambertw_series` and
  `acb_poly_lambertw_series` use the scalar function for the constant
  term and Newton iteration on power series for the other coefficients.
- `d_lambertw` of [`double_extras`](double_extras) and the internal
  `d_lambertw_branch1` provide the double-precision starting values;
  other functions use them to estimate numbers of terms, for example
  the error function, the cosine integral and the gamma function of
  rational arguments in [`arb_hypgeom`](arb_hypgeom).
