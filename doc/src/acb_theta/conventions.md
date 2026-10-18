---
title: Characteristics, derivatives and precision
summary: How theta characteristics are numbered, how derivatives are normalised and ordered, how the symplectic group acts, and what the precision argument means for theta functions.
order: 1
---

This guide explains the conventions shared by all functions of
[`acb_theta`](acb_theta). The definitions follow [@Mum1983] and
[@Igu1972]; the algorithms are in
[](algorithms/special-functions/theta-functions).

## Theta functions with characteristics {#definition}

For $\tau \in \mathcal{H}_g$ (a symmetric complex $g \times g$ matrix
with positive definite imaginary part $Y$), $z \in \mathbb{C}^g$ and
$a, b \in \{0, 1\}^g$, viewed as column vectors,

$$
\theta_{a,b}(z,\tau) = \sum_{n\in \Z^{g} + \tfrac a2} \exp\bigl(\pi i n^T \tau n + 2\pi i n^T (z + \tfrac b2)\bigr).
$$

The $2^{2g}$ functions are the theta functions of level 2. Shifting $z$
by half-periods permutes them up to explicit exponential factors:

$$
\theta_{a,b}(z,\tau) = \exp\bigl(\pi i a^T\tau a/4 + \pi i a^T(z + \tfrac b2)\bigr) \, \theta_{0,0}\bigl(z + \tau\tfrac{a}{2} + \tfrac{b}{2}, \tau\bigr).
$$

The characteristic is *even* if $a^T b$ is even and *odd* otherwise.
$\theta_{a,b}(-z, \tau) = (-1)^{a^T b} \theta_{a,b}(z, \tau)$, so the odd
theta functions vanish at $z = 0$; there are $2^{g-1}(2^g + 1)$ even and
$2^{g-1}(2^g - 1)$ odd characteristics.

For $g = 1$, in the notation of [`acb_modular`](acb_modular),
$\theta_{0,0} = \theta_3$, $\theta_{0,1} = \theta_4$, $\theta_{1,0} = \theta_2$
and $\theta_{1,1} = -\theta_1$.

## Numbering of characteristics {#numbering}

A characteristic $(a, b)$ is stored as the `ulong` whose binary
expansion, most significant bit first, is
$a_0 a_1 \ldots a_{g-1} b_0 b_1 \ldots b_{g-1}$; vectors of $2^{2g}$
theta values are ordered by this number. For $g = 2$, $(a, b) = ((1,0), (0,0))$
is number $8 = 1000_2$, and the even characteristics are those with
numbers 0, 1, 2, 3, 4, 6, 8, 9, 12, 15. When $a$ or $b$ is handled alone
(vectors of $2^g$ values $\theta_{a,0}$, say), it is numbered on $g$
bits in the same way. `acb_theta_char_set_slong_vec` converts a vector
of bits to this number and `acb_theta_char_bit` extracts a bit.

## Derivatives {#derivatives}

Only derivatives with respect to $z$ are computed: those with respect to
$\tau$ follow from the heat equation
$2\pi i (1 + \delta_{j,k}) \, \partial\theta_{a,b}/\partial\tau_{j,k} = \partial^2\theta_{a,b}/\partial z_j \partial z_k$.
A *derivation tuple* $k = (k_0, \ldots, k_{g-1})$ of nonnegative
integers corresponds to the normalised derivative

$$
\frac{1}{k_0! \cdots k_{g-1}!} \frac{\partial^{|k|}\theta_{a,b}}{\partial z_0^{k_0}\cdots \partial z_{g-1}^{k_{g-1}}}, \qquad |k| = k_0 + \cdots + k_{g-1},
$$

which is the coefficient of $h_0^{k_0} \cdots h_{g-1}^{k_{g-1}}$ in the
Taylor expansion of $\theta_{a,b}(z + h, \tau)$, as in
`acb_modular_theta_jet`. The derivatives up to total order `ord`, of
which there are $\binom{\mathit{ord}+g}{g}$ (`acb_theta_jet_nb`), are
ordered by total order and then reverse-lexicographically; for $g = 2$:
$(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), (3,0), \ldots$.
`acb_theta_jet_tuples` lists them and `acb_theta_jet_index` finds the
position of a tuple. The functions `acb_theta_jet_mul`,
`acb_theta_jet_compose`, `acb_theta_jet_exp_pi_i` and
`acb_theta_jet_exp_qf` compute with such vectors for arbitrary
functions.

## The symplectic group {#symplectic}

$\mathrm{Sp}_{2g}(\Z)$ is the group of integer matrices
$m = \begin{pmatrix} \alpha & \beta \\ \gamma & \delta \end{pmatrix}$
with $m^T J m = J$, $J = \begin{pmatrix} 0 & I_g \\ -I_g & 0 \end{pmatrix}$.
It acts on $\mathbb{C}^g \times \mathcal{H}_g$ by

$$
m \cdot (z, \tau) = \bigl((\gamma\tau + \delta)^{-T} z, \; (\alpha\tau + \beta)(\gamma\tau + \delta)^{-1}\bigr),
$$

and the theta functions transform as

$$
\theta_{a,b}(m\cdot(z,\tau)) = \kappa(m) \zeta_8^{e(m, a, b)} \det(\gamma\tau + \delta)^{1/2} e^{\pi i z^T (\gamma\tau + \delta)^{-1} \gamma z} \theta_{a',b'}(z,\tau),
$$

where $\zeta_8 = e^{\pi i/4}$, $(a', b')$ and $e(m, a, b)$ are given by
`acb_theta_char_table`, and $\kappa(m)$, an 8th root of unity that
depends on the branch of the square root, by `acb_siegel_kappa`
[@Igu1972], [@Mum1983]. Every $\tau$ is equivalent to a *reduced* matrix,
for which the theta series converge fast and uniformly; see
[](algorithms/special-functions/theta-functions#reduction).

## Precision {#precision}

The magnitude of $\theta_{a,b}(z, \tau)$ is, barring cancellation, about

$$
|\theta_{a,b}(z,\tau)| \approx e^{\pi y^T Y^{-1} y} e^{-d^2},
$$

where $y = \operatorname{Im}(z)$ and $d$ is the distance from
$-Y^{-1}y$ to the shifted lattice $\Z^g + \tfrac a2$ for the norm
$\lVert x \rVert_\tau^2 = \pi x^T Y x$: the largest term of the series is
the one with $n$ closest to $-Y^{-1}y$. The module therefore works with
the normalised function

$$
\widetilde{\theta}_{a,b}(z,\tau) = e^{-\pi y^T Y^{-1} y} \, \theta_{a,b}(z,\tau),
$$

which is not holomorphic but is bounded on $\mathbb{C}^g$ for fixed $\tau$, and
satisfies the same duplication formulas. A precision argument `prec`
means that $\widetilde{\theta}_{a,b}(z,\tau)$ is computed with an
absolute error of about $2^{-\mathit{prec}}$, so the outputs of
`acb_theta_jet` and `acb_theta_all` have absolute errors of about
$e^{\pi y^T Y^{-1} y} 2^{-\mathit{prec}}$ (for the reduced arguments).
This avoids unbounded work when $y$ is far from zero. Some internal
functions (with a `distances` argument or a `shifted_prec` flag) also
take the factor $e^{-d^2}$ into account and compute
$\widetilde\theta_{a,b}$ with absolute error about
$e^{-d^2} 2^{-\mathit{prec}}$, a *shifted absolute precision*; the
squared distances are computed by `acb_theta_eld_distances`.

All outputs are certified enclosures. If the imaginary part of $\tau$ is
not certainly positive definite at the working precision, the outputs
have infinite radius.

## Testing and profiling programs {#profiling}

The test programs in `src/acb_theta/test` check, among others, that the
transformation formula and the duplication formulas hold on random
input, that `acb_siegel_reduce` produces matrices accepted by
`acb_siegel_is_reduced`, and that the quasi-linear algorithm agrees with
summation. The programs in `src/acb_theta/profile` are used for tuning:

- `p-sum g prec all_a all_b tilde skew` times `acb_theta_sum` (distances,
  contexts, summation) on random input whose last row and column of
  $\tau$ are multiplied by `skew`;
- `p-ql_setup g nb_steps all skew` checks that `acb_theta_ql_setup`
  stays cheap and succeeds as the precision and number of steps grow;
- `p-ql_exact g pmin pmax cst s exp` compares, for precisions from
  `pmin` to `pmax`, the pattern suggested by `acb_theta_ql_nb_steps`
  with neighbouring patterns in `acb_theta_ql_exact`, on matrices whose
  last `s` rows and columns are scaled by $2^{\mathit{exp}}$;
- `p-ql_jet_fd g prec ord all` compares `acb_theta_ql_jet_fd` with
  `acb_theta_sum_jet`, which decides the cutoffs of `acb_theta_ql_jet`.
