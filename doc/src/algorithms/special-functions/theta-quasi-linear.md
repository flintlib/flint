---
title: Quasi-linear evaluation of theta functions
summary: The algorithm of Elkies and Kieffer for Riemann theta functions in any genus: duplication formulas as a generalised AGM, square roots chosen from low-precision summations, an auxiliary real shift for robustness, dimension lowering for skewed matrices, and derivatives by finite differences.
complexity: $O(M(p) \log p)$ operations at precision $p$ for theta values at a reduced $(z, \tau)$, for fixed genus $g$ (the constant grows exponentially with $g$); derivatives up to order $k$ cost about $(k+1)^g$ times as much.
order: 66
---

Summing the theta series over an ellipsoid costs $O(p^{g/2})$
multiplications at precision $p$ in genus $g$
([](algorithms/special-functions/theta-functions#summation)). For
$g = 1$, the theta constants can instead be computed with the
arithmetic-geometric mean in $O(\log p)$ steps, and Dupont [@Dup2006],
and Labrande and Thomé [@LT2016] for $g \ge 2$, generalised this to
quasi-linear algorithms based on the duplication formulas and Newton's
method, which however needed hypotheses on the input (good choices of
square roots) and were not uniform in $\tau$. Elkies and Kieffer
[@EK2025] give an algorithm that is uniform in $(z, \tau)$ for every
fixed $g$ and needs no Newton iteration. It is implemented in
`acb_theta_ql_exact` and used by every evaluation in
[`acb_theta`](acb_theta) above low precision.

Throughout, $(z, \tau)$ is reduced
([](algorithms/special-functions/theta-functions#reduction)), and the
characteristics $a$ run over $\{0,1\}^g$ with addition modulo 2.

## Duplication formulas {#duplication}

::: theorem Duplication formula {#thm-duplication}
For $z, z' \in \mathbb{C}^g$, $\tau \in \mathcal{H}_g$ and $a \in \{0,1\}^g$,

$$
\theta_{a,0}(z,\tau)\,\theta_{a,0}(z',\tau) = \sum_{a' \in \{0,1\}^g} \theta_{a',0}(z+z',2\tau)\,\theta_{a+a',0}(z-z',2\tau),
$$

and for any $b$,
$\theta_{a,b}(z,\tau)\theta_{a,b}(z',\tau) = \sum_{a'} (-1)^{a'^T b}\, \theta_{a',0}(z+z',2\tau)\,\theta_{a+a',0}(z-z',2\tau)$.
:::

::: proof
The product is the sum over $m, n \in \Z^g + \frac a2$ of
$\exp(\pi i (m^T\tau m + n^T \tau n) + 2\pi i (m^T z + n^T z'))$. Put
$P = (m+n)/2$ and $Q = (m-n)/2$. Then
$m^T\tau m + n^T\tau n = P^T (2\tau) P + Q^T (2\tau) Q$ and
$m^T z + n^T z' = P^T(z + z') + Q^T(z - z')$. As $(m, n)$ runs over
$(\Z^g + \frac a2)^2$, $(P, Q)$ runs exactly once over the pairs with
$P \in \Z^g + \frac{a'}2$ and $Q \in \Z^g + \frac{a + a'}{2}$ for some
$a'$: indeed $P - Q = n \in \Z^g + \frac a2$, and conversely
$m = P + Q$, $n = P - Q$ lie in $\Z^g + \frac a2$. The sum factors as
stated. For the general case, replace $z$ and $z'$ by $z + \frac b2$ and
$z' + \frac b2$: then $z - z'$ is unchanged and $z + z'$ is shifted by
$b$, which multiplies the term $P$ by
$e^{2\pi i P^T b} = (-1)^{a'^T b}$.
:::

The formula also holds for the normalised functions
$\widetilde\theta_{a,0}(z, \tau) = e^{-\pi y^T Y^{-1} y}\theta_{a,0}(z, \tau)$,
because the normalising factors match on both sides. Three special cases
are used:

$$
\begin{aligned}
\theta_{a,0}(z,\tau)^2 &= \sum_{a'} \theta_{a',0}(2z,2\tau)\, \theta_{a+a',0}(0,2\tau), \\
\theta_{a,0}(0,\tau)\,\theta_{a,0}(z,\tau) &= \sum_{a'} \theta_{a',0}(z,2\tau)\, \theta_{a+a',0}(z,2\tau), \\
\theta_{a,0}(0,\tau)^2 &= \sum_{a'} \theta_{a',0}(0,2\tau)\, \theta_{a+a',0}(0,2\tau).
\end{aligned}
$$

For $g = 1$, the third formula and its version with $b = 1$ give
$\theta_{0,0}(0,\tau)^2 \pm \theta_{0,1}(0,\tau)^2 = 2\theta_{a,0}(0,2\tau)^2$
(with $a = 0$ for the sign $+$ and $a = 1$ for $-$), and the second
formula at $z = \frac12$, where $\theta_{0,0}(\frac12, \tau) = \theta_{0,1}(0, \tau)$
and the odd $\theta_{1,1}(0, 2\tau)$ vanishes, gives
$\theta_{0,1}(0, 2\tau)^2 = \theta_{0,0}(0,\tau)\theta_{0,1}(0,\tau)$.
Hence $\alpha_n = \theta_{0,0}(0, 2^n\tau)^2$ and
$\beta_n = \theta_{0,1}(0, 2^n\tau)^2$ satisfy
$\alpha_{n+1} = (\alpha_n + \beta_n)/2$ and $\beta_{n+1} = \sqrt{\alpha_n \beta_n}$:
the AGM, read from $2^m\tau$ down to $\tau$ by taking square roots. In general, each formula is a sum of $2^g$ products,
computed for all $a$ at once by two Hadamard transforms of length $2^g$
and $2^g$ pointwise products (`acb_theta_agm_mul`, following [@LT2016]).

## The algorithm for theta constants {#constants}

To compute $\theta_{a,0}(0, \tau)$ for all $a$, choose $m$ with
$2^m \approx p$ (more precisely, so that the ellipsoid for $2^m\tau$ has
$O(1)$ points: $2^m \gamma_j^2 \approx p$ for the diagonal entries
$\gamma_j$ of the Cholesky factor of $\pi Y$).

1. Compute $\theta_{a,0}(0, 2^m\tau)$ for all $a$ by summation. Since
   $2^m Y$ is large, only $O_g(1)$ lattice points contribute at
   precision $p$, and the value is about $e^{-2^m d_a^2}$ where $d_a$ is
   the distance from 0 to $\Z^g + \frac a2$.
2. For $k = m - 1, \ldots, 0$: apply the third duplication formula to
   get $\theta_{a,0}(0, 2^k\tau)^2$ from the values at $2^{k+1}\tau$, and
   take square roots, choosing the signs with low-precision
   approximations of $\theta_{a,0}(0, 2^k\tau)$ computed by summation.

The cost is $m = O(\log p)$ steps of $O(2^g)$ multiplications and square
roots at precision $O(p)$, that is $O(M(p) \log p)$, provided the
precision loss per step is bounded.

**Precision.** The values $\theta_{a,0}(0, 2^k\tau)$ for $a \ne 0$ are
exponentially small. FLINT works with *shifted absolute precision*: a
value of expected magnitude $e^{-d^2}$ is computed with absolute error
$e^{-d^2} 2^{-p}$. The duplication formula is compatible with this,
because by the parallelogram identity the products
$\theta_{a',0}\theta_{a+a',0}$ at $2\tau$ have magnitude at most that of
$\theta_{a,0}^2$ at $\tau$ (see `acb_theta_agm_mul_tight`, which computes
the products on midpoints and adds error bounds of the right shape).
Assuming that every $|\theta_{a,0}(0, 2^k\tau)|$ is of its expected
order of magnitude, the loss is $O_g(1)$ bits per step, and the square
roots are determined by approximations to $O_g(1)$ bits.

## Robustness: an auxiliary shift {#shift}

The assumption can fail: some $\theta_{a,0}(0, 2^k\tau)$ may be much
smaller than expected (for instance, close to a zero of a theta
constant), and then a square root cannot be chosen, or precision is
lost. The algorithm then introduces a random real vector $t$. A real
shift does not change the expected magnitudes, and for almost every $t$
none of the values $\theta_{a,0}(2^k t, 2^k\tau)$ and
$\theta_{a,0}(2^{k+1} t, 2^k\tau)$ is abnormally small. At each step,
the first duplication formula gives $\theta_{a,0}(2^k t, 2^k\tau)$ and
$\theta_{a,0}(2^{k+1}t, 2^k\tau)$ by square roots (from the values at
$2^{k+1}\tau$ of the vectors $0, 2^{k+1}t, 2^{k+2}t$), and the second
formula with $z = 2^k t$ gives the theta constants
$\theta_{a,0}(0, 2^k\tau)$ by divisions, without square roots of small
numbers. `acb_theta_ql_setup` tries a few random vectors $t$ at low
precision, doubling the precision if none works, and determines how
many initial steps are *easy* (need no shift). This makes the algorithm
probabilistic, with a negligible probability of failure, in which case
the result is indeterminate rather than wrong.

**Nonzero $z$.** For theta values at $z \ne 0$, the first formula
computes $\theta_{a,0}(2^k z, 2^k\tau)$ from the values at
$(2^{k+1}z, 2^{k+1}\tau)$ and the theta constants, so the algorithm
carries the vectors $2^k v$ for $v \in \{0, t, 2t, z, z+t, z+2t\}$, of
which $z$ itself is only needed at the last step. Theta values with
nonzero $b$ are obtained at the last step from the general duplication
formula with signs $(-1)^{a'^T b}$.

## Dimension lowering {#lower-dim}

If the eigenvalues of $Y$ have very different sizes, doubling $\tau$
makes the ellipsoids thin in some directions long before they are thin
in others, and the number of duplication steps needed for the largest
direction would make the precision loss in the others large. The
algorithm then splits the variables: writing
$\tau = \begin{pmatrix} \tau_0 & x \\ x^T & \tau_1 \end{pmatrix}$ with
$\tau_0$ of size $s$ and $n = (n_0, n_1)$, the theta series is

$$
\theta_{a,0}(z, \tau) = \sum_{n_1 \in \Z^{g-s} + a_1/2} e^{\pi i (n_1^T \tau_1 n_1 + 2 n_1^T z_1)} \, \theta_{a_0,0}(z_0 + x n_1, \tau_0),
$$

and only the $n_1$ in a small ellipsoid of dimension $g - s$ (with the
lower right block of the Cholesky factor) contribute at the required
precision. Each term is a theta value in dimension $s$ at a shifted
vector, which is evaluated recursively by the same algorithm for many
vectors at once. `acb_theta_ql_nb_steps` decides, from the diagonal of
the Cholesky factor and the precision, how many duplication steps to
apply before lowering the dimension at each level (a vector called the
*pattern*). `acb_theta_ql_lower_dim` and `acb_theta_ql_recombine` perform the split
and the recombination, and `acb_theta_ql_exact` runs the whole recursion.
The effect is that every summation and every duplication step works at
precision $O(p)$ on ellipsoids with $O(1)$ points, which keeps the
algorithm quasi-linear uniformly in $\tau$.

## The algorithm {#algorithm}

::: algorithm Theta values at a reduced point (`acb_theta_ql_exact`) {#alg-ql}
**Input:** exact reduced $(z_i, \tau)$, $i = 1, \ldots, \mathit{nb}$, with
$z_1 = 0$; a pattern from `acb_theta_ql_nb_steps`; precision $p$.

**Output:** $\widetilde\theta_{a,0}(z_i, \tau)$ for all $a$ (or all
$\widetilde\theta_{a,b}$).

1. Let $m$ be the number of duplication steps of the pattern. If
   $m = 0$, sum the series and return.
2. With `acb_theta_ql_setup`, compute low-precision approximations of all
   the values whose square roots will be taken, choose the auxiliary
   vector $t$ if needed, and estimate the precision loss $\ell$.
3. Compute the theta values at $2^m\tau$ for the vectors $2^m v$ at
   precision $p + \ell$: by summation, or, if the pattern says so, by
   dimension lowering and a recursive call.
4. For $k = m - 1, \ldots, 0$, apply the duplication formulas
   (`acb_theta_agm_mul_tight`), taking square roots with
   `acb_theta_agm_sqrt` guided by the approximations, or divisions by the
   shifted values.
5. At $k = 0$, compute the requested values (with the general formula
   for $b \ne 0$).
:::

The main function `acb_theta_jet` reduces $(z, \tau)$, uses
$\theta_{0,0}(z, \tau) = \sum_a \theta_{a,0}(2z, 4\tau)$ when only one
value is needed, and
$\theta_{a,b}(z,\tau)^2 = \sum_{a'}(-1)^{a'^T b}\theta_{a',0}(0, 2\tau)\theta_{a+a',0}(2z, 2\tau)$
when the squares of all values are requested, so that the quasi-linear
algorithm only ever computes the $2^g$ functions $\theta_{a,0}$.

## Derivatives {#derivatives}

Partial derivatives with respect to $z$ up to total order $k$ are
computed from values (`acb_theta_ql_jet_fd`). Let $\zeta$ be a primitive
$(k+1)$-th root of unity, $\varepsilon > 0$, and
$h_n = (\varepsilon\zeta^{n_0}, \ldots, \varepsilon\zeta^{n_{g-1}})$ for
$n \in \{0, \ldots, k\}^g$. If $\theta(z + h) = \sum_{j \ge 0} c_j h^j$
(multi-index notation), a discrete Fourier transform in each coordinate
gives

$$
\sum_{n} \zeta^{-q^T n} \theta(z + h_n) = (k+1)^g \sum_{j \equiv q \bmod (k+1)} c_j \varepsilon^{|j|},
$$

which is $c_q \varepsilon^{|q|}$ plus terms of order
$\varepsilon^{|q| + k + 1}$. If $|\theta| \le c$ on the polydisc of
radius $\rho$ around $z$, Cauchy's estimate bounds those terms, and the
error of $c_q$ is at most $2cg\,\varepsilon^{k+1}/\rho^{k+1+|q|}$.
`acb_theta_ql_local_bound` gives $c$ and $\rho$ from the Cholesky factor
and $y$; $\varepsilon$ is chosen so that this error is $2^{-p}$, and the
precision of the $(k+1)^g$ evaluations is raised to compensate for the
division by $\varepsilon^{|q|}$. The evaluations are made at the
midpoints of $(z, \tau)$, and the effect of the input radii is added
with `acb_theta_ql_jet_error`, which bounds the variation of the
derivatives over the input balls with derivatives of order $k + 2$
(computed cheaply at low precision) and the heat equation. The total cost
is quasi-linear in $p (k+1)^g$. At low precision, direct summation of the
differentiated series is faster, and `acb_theta_ql_jet` chooses between
the two.

## Use in FLINT {#use}

- `acb_theta_agm_mul`, `acb_theta_agm_mul_tight`, `acb_theta_agm_sqrt`:
  duplication steps.
- `acb_theta_ql_nb_steps`, `acb_theta_ql_setup`, `acb_theta_ql_lower_dim`,
  `acb_theta_ql_recombine`, `acb_theta_ql_exact`: the algorithm for
  values.
- `acb_theta_ql_local_bound`, `acb_theta_ql_jet_error`,
  `acb_theta_ql_jet_fd`, `acb_theta_ql_jet`: derivatives.
- `acb_theta_jet_notransform`, `acb_theta_jet`, `acb_theta_one`,
  `acb_theta_all`, and the genus 2 modular forms
  (`acb_theta_g2_sextic_chi5`) built on them.
