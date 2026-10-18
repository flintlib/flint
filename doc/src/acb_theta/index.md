---
title: acb_theta.h
header: acb_theta.h
description: Riemann theta functions and their derivatives in any genus, with a uniform quasi-linear algorithm, the Siegel modular group, and Siegel modular forms of degree 2.
category: Special functions
order: 5
---

This module evaluates the Riemann theta functions with characteristics

$$
\theta_{a,b}(z,\tau) = \sum_{n\in \Z^{g} + \tfrac a2} \exp\bigl(\pi i n^T \tau n + 2\pi i n^T (z + \tfrac b2)\bigr),
\qquad a, b \in \{0, 1\}^g,
$$

and their partial derivatives with respect to $z$, for any genus
$g \ge 1$, with rigorous error bounds. Here $\tau$ is in the Siegel upper
half-space $\mathcal{H}_g$ of symmetric complex $g \times g$ matrices
with positive definite imaginary part, and $z \in \mathbb{C}^g$. The main
function `acb_theta_jet` reduces $(z, \tau)$ with the action of the
symplectic group $\mathrm{Sp}_{2g}(\Z)$, evaluates the theta functions at
the reduced point with the algorithm of Elkies and Kieffer [@EK2025],
which combines duplication formulas with short summations and runs in
time quasi-linear in the precision for every fixed $g$, and applies the
theta transformation formula. The module also provides arithmetic in
$\mathrm{Sp}_{2g}(\Z)$, reduction in $\mathcal{H}_g$, the summation and
duplication building blocks, and, for $g = 2$, Siegel modular forms and
covariants of binary sextics computed from theta values.

For $g = 1$, the theta functions are those of
[`acb_modular`](acb_modular), whose summation code this module uses; in
that notation,
$(\theta_{0,0}, \theta_{0,1}, \theta_{1,0}, \theta_{1,1}) = (\theta_3, \theta_4, \theta_2, -\theta_1)$.

Start with the guide [](acb_theta/conventions), which explains the
numbering of characteristics and derivatives and the meaning of the
precision argument, and the function page [](acb_theta/jet). The
algorithms are described in [](algorithms/special-functions/theta-functions)
and [](algorithms/special-functions/theta-quasi-linear).

::: summary algorithms/special-functions/theta-quasi-linear :::

## Conventions {#conventions}

- $\tau$ is an `acb_mat_t` of size $g \times g$; $z$ is an `acb_ptr` of
  length $g$, and a list of $\mathit{nb}$ vectors `zs` has length
  $\mathit{nb} \cdot g$. $Y = \operatorname{Im}(\tau)$ and
  $y = \operatorname{Im}(z)$.
- A characteristic $(a, b)$ is encoded as the integer whose $2g$ binary
  digits, most significant first, are $a_0, \ldots, a_{g-1}, b_0, \ldots, b_{g-1}$,
  and vectors of $2^{2g}$ theta values are ordered by this integer.
  Characteristics $a$ or $b$ alone are encoded on $g$ bits.
- Derivatives are taken with respect to $z_0, \ldots, z_{g-1}$ and
  normalised as Taylor coefficients; derivation tuples are ordered by
  total order, then reverse-lexicographically ([](#jets)).
- The working precision is an *absolute* precision for the normalised
  function $\widetilde\theta_{a,b}(z, \tau) = e^{-\pi y^T Y^{-1} y} \theta_{a,b}(z,\tau)$;
  see [](acb_theta/conventions#precision).
- Outputs are certified: if $Y$ is not certainly positive definite at the
  working precision, or the reduction fails, the outputs are
  indeterminate.
- Symplectic matrices are `fmpz_mat_t` of size $2g \times 2g$, written
  $m = \begin{pmatrix} \alpha & \beta \\ \gamma & \delta \end{pmatrix}$
  with $g \times g$ blocks; they act by
  $m \cdot \tau = (\alpha\tau + \beta)(\gamma\tau + \delta)^{-1}$ and
  $m \cdot (z, \tau) = ((\gamma\tau + \delta)^{-T} z, m \cdot \tau)$.

::: macro
ACB_THETA_LOW_PREC

The low precision (32 bits) used internally for estimates such as
Cholesky factors and choices of parameters.
:::

## Theta functions {#main}

::: function
void acb_theta_jet(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong ord, ulong ab, int all, int sqr, slong prec)

Sets `th` to the partial derivatives with respect to $z$, up to total
order `ord`, of $\theta_{a,b}(z, \tau)$ at each of the `nb` vectors $z$
in `zs`: for the single characteristic `ab` if `all` is zero, or for all
$2^{2g}$ characteristics if `all` is nonzero. With `ord` = 0 and `sqr`
nonzero, it computes the squares $\theta_{a,b}(z,\tau)^2$ instead, which
is faster when `all` is set. See [](acb_theta/jet).
:::

::: function
void acb_theta_one(acb_ptr th, acb_srcptr z, const acb_mat_t tau, ulong ab, slong prec)
void acb_theta_all(acb_ptr th, acb_srcptr z, const acb_mat_t tau, int sqr, slong prec)

Set `th` to the single value $\theta_{a,b}(z,\tau)$ for the characteristic
`ab`, respectively to the vector of $2^{2g}$ values $\theta_{a,b}(z,\tau)$
(or their squares if `sqr` is nonzero) for all characteristics. These are
`acb_theta_jet` with `nb` = 1 and `ord` = 0; computing all values costs
little more than computing one.
:::

::: function
void acb_theta_jet_notransform(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong ord, ulong ab, int all, int sqr, slong prec)

Same as `acb_theta_jet`, but without reducing $(z, \tau)$; the input
should already be reduced.

::: details
Before calling `acb_theta_ql_jet`, the function uses these formulas:

- if `ord` is zero and `all` and `sqr` are nonzero, the duplication
  formula
  $\theta_{a,b}(z,\tau)^2 = \sum_{a'} (-1)^{a'^T b} \theta_{a',0}(0,2\tau) \theta_{a+a',0}(2z,2\tau)$,
  which needs only the $2^g$ functions $\theta_{a,0}$ at $2\tau$;
- if `all` is zero and `ab` is zero,
  $\theta_{0,0}(z,\tau) = \sum_{a} \theta_{a,0}(2z, 4\tau)$;
- if `all` is zero and `ab` is nonzero,
  $\theta_{a,b}(z,\tau) = \exp\bigl(\pi i a^T\tau a/4 + \pi i a^T(z + \tfrac b2)\bigr) \theta_{0,0}(z + \tau\tfrac{a}{2} + \tfrac{b}{2}, \tau)$,
  which reduces to the previous case (for derivatives, the exponential
  factor is expanded with `acb_theta_jet_exp_pi_i`); with `sqr`, the
  value is squared at the end.
:::
:::

::: function
int acb_theta_reduce_tau(acb_ptr new_zs, acb_mat_t new_tau, fmpz_mat_t mat, acb_mat_t N, acb_mat_t ct, acb_ptr exps, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)

Reduces $\tau$ under $\mathrm{Sp}_{2g}(\Z)$ and transforms the `nb`
vectors in `zs` accordingly. Returns 1 on success and 0 on failure
(malformed input, or a precision too low to see that $Y$ is positive
definite), in which case the outputs are undefined.

::: details
On success, `new_tau` $= \mathit{mat} \cdot \tau$ is the reduced matrix
and `mat` the symplectic matrix found by `acb_siegel_reduce`, `ct` is the
transpose of $(\gamma\tau + \delta)^{-1}$, $N$ is the matrix
$\pi i \gamma \cdot \mathit{ct}$ of the exponential factor of the
transformation formula, and `new_zs` holds the vectors
$\mathit{ct} \cdot z$. The vector `exps` (of length `nb`) is used as
workspace.
:::
:::

::: function
int acb_theta_reduce_z(acb_ptr new_zs, arb_ptr rs, acb_ptr cs, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)

Reduces the `nb` vectors in `zs` with respect to the lattice
$\Z^g + \tau \Z^g$. For each $z$, $Y^{-1} y$ is rounded to the nearest
even integral vector $r$ (stored in `rs`), and $z' = z - \tau r$ minus
the nearest even integral vector to its real part is stored in `new_zs`;
then $\theta_{a,b}(z',\tau) = e^{-\pi i r^T (z + z')} \theta_{a,b}(z,\tau)$
for all characteristics, and the exponential factor is stored in `cs`.
Returns 0 (with undefined outputs) if the rounding fails because of
extreme values, and 1 otherwise.
:::

## The Siegel modular group {#sp2gz}

Symplectic matrices are ordinary `fmpz_mat_t`, so the functions of
[`fmpz_mat`](fmpz_mat) (such as `fmpz_mat_equal` and `fmpz_mat_mul`)
apply to them. Except for `sp2gz_is_correct`, the functions assume that
the input is square of even size $2g$.

::: function
slong sp2gz_dim(const fmpz_mat_t mat)

Returns $g$, half the number of rows of `mat`.
:::

::: function
void sp2gz_set_blocks(fmpz_mat_t mat, const fmpz_mat_t alpha, const fmpz_mat_t beta, const fmpz_mat_t gamma, const fmpz_mat_t delta)

Sets `mat` to $\begin{pmatrix} \alpha & \beta \\ \gamma & \delta \end{pmatrix}$.
The dimensions must match.
:::

::: function
void sp2gz_j(fmpz_mat_t mat)
void sp2gz_block_diag(fmpz_mat_t mat, const fmpz_mat_t U)
void sp2gz_trig(fmpz_mat_t mat, const fmpz_mat_t S)

Set `mat` to $J = \begin{pmatrix} 0 & I_g \\ -I_g & 0 \end{pmatrix}$, to
the block-diagonal matrix $\begin{pmatrix} U & 0 \\ 0 & U^{-T} \end{pmatrix}$
for $U \in \mathrm{GL}_g(\Z)$, and to the trigonal matrix
$\begin{pmatrix} I_g & S \\ 0 & I_g \end{pmatrix}$ for a symmetric
integer matrix $S$, respectively.
:::

::: function
void sp2gz_embed(fmpz_mat_t res, const fmpz_mat_t mat)
void sp2gz_restrict(fmpz_mat_t res, const fmpz_mat_t mat)

`sp2gz_embed` sets the $2g \times 2g$ matrix `res` to the embedding of
the $2r \times 2r$ symplectic matrix `mat` ($r \le g$) acting on the first
$r$ coordinates:
$\begin{pmatrix} \alpha & 0 & \beta & 0 \\ 0 & I_{g-r} & 0 & 0 \\ \gamma & 0 & \delta & 0 \\ 0 & 0 & 0 & I_{g-r} \end{pmatrix}$.
`sp2gz_restrict` sets the $2r \times 2r$ matrix `res` to the matrix whose
$r \times r$ blocks are the upper left corners of the blocks of `mat`;
the result need not be symplectic.
:::

::: function
slong sp2gz_nb_fundamental(slong g)
void sp2gz_fundamental(fmpz_mat_t mat, slong j)

`sp2gz_nb_fundamental` returns the number of fundamental symplectic
matrices used in the reduction algorithm, and `sp2gz_fundamental` sets
`mat` to the $j$-th of them. For $g = 1$ there is one ($J$) and for
$g = 2$ there are 19, which define the boundary of Gottschling's
fundamental domain [@Got1959]. For $g > 2$ no such complete set is
known; the module uses the $19 g(g-1)/2$ matrices obtained by applying
the $g = 2$ matrices to any pair of coordinates and the $2^g$ embeddings
of a lower-dimensional $J$ on any subset of coordinates.
:::

::: function
int sp2gz_is_correct(const fmpz_mat_t mat)
int sp2gz_is_j(const fmpz_mat_t mat)
int sp2gz_is_block_diag(const fmpz_mat_t mat)
int sp2gz_is_trig(const fmpz_mat_t mat)

Return whether `mat` is a symplectic matrix, respectively whether the
symplectic matrix `mat` is $J$, block-diagonal or trigonal as above.
:::

::: function
int sp2gz_is_embedded(fmpz_mat_t res, const fmpz_mat_t mat)

Given a $2g \times 2g$ symplectic matrix `mat` and a square matrix `res`
of size $2r$ with $r \le g$, returns whether `mat` is the result of
`sp2gz_embed` on a $2r \times 2r$ symplectic matrix, and if so stores
that matrix in `res` (otherwise `res` is undefined).
:::

::: function
void sp2gz_inv(fmpz_mat_t inv, const fmpz_mat_t mat)

Sets `inv` to the inverse
$\begin{pmatrix} \delta^T & -\beta^T \\ -\gamma^T & \alpha^T \end{pmatrix}$
of the symplectic matrix `mat`, without any arithmetic.
:::

::: function
fmpz_mat_struct * sp2gz_decompose(slong * nb, const fmpz_mat_t mat)

Returns a newly allocated vector of symplectic matrices whose product
from left to right is `mat`, and stores its length in `*nb`. Each matrix
is block-diagonal, trigonal, $J$, an embedded $J$ from a lower dimension,
or an embedded matrix from dimension 1. The caller frees the result by
calling `fmpz_mat_clear` on each entry and then `flint_free`.
:::

::: function
void sp2gz_randtest(fmpz_mat_t mat, flint_rand_t state, slong bits)

Sets `mat` to a random symplectic matrix with entries of about `bits`
bits, a product of block-diagonal and trigonal matrices and $J$.
:::

## The Siegel half-space {#siegel}

::: function
void acb_siegel_cocycle(acb_mat_t c, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)
void acb_siegel_transform_cocycle_inv(acb_mat_t w, acb_mat_t c, acb_mat_t cinv, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)
void acb_siegel_transform(acb_mat_t w, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)

`acb_siegel_cocycle` sets `c` to the cocycle $\gamma\tau + \delta$;
`acb_siegel_transform` sets `w` to
$\mathit{mat} \cdot \tau = (\alpha\tau + \beta)(\gamma\tau + \delta)^{-1}$;
`acb_siegel_transform_cocycle_inv` computes both and also
$(\gamma\tau + \delta)^{-1}$.
:::

::: function
void acb_siegel_cho_yinv(arb_mat_t cho, arb_mat_t yinv, const acb_mat_t tau, slong prec)

Sets `yinv` to $Y^{-1}$ and `cho` to the upper triangular Cholesky
matrix $C$ with positive diagonal and $\pi Y = C^T C$. If $Y$ is not
certainly positive definite at the working precision, both are set to
indeterminate matrices.
:::

::: function
void acb_siegel_reduce(fmpz_mat_t mat, const acb_mat_t tau, slong prec)

Sets `mat` to a symplectic matrix such that $\mathit{mat} \cdot \tau$ is
as reduced as possible, by repeatedly reducing $Y$ with LLL
(block-diagonal matrices), the real part of $\tau$ (trigonal matrices),
and applying fundamental matrices that increase $|\det Y|$. If the
entries of $\tau$ are not of reasonable size or $\det Y$ is vanishingly
small, `mat` is set to the identity.
:::

::: function
int acb_siegel_is_reduced(const acb_mat_t tau, slong tol_exp, slong prec)

Returns nonzero if $\tau$ is certainly in the reduced domain with
tolerance $\varepsilon = 2^{\mathit{tol\_exp}}$:
$|\operatorname{Re}(\tau_{j,k})| \le \frac12 + \varepsilon$ for all $j, k$,
$Y$ passes `arb_mat_spd_is_lll_reduced` with parameters
$(0.99 - \varepsilon, 0.51 + \varepsilon)$, and
$|\det(\gamma\tau + \delta)| \ge 1 - \varepsilon$ for every fundamental
matrix of `sp2gz_fundamental`.
:::

::: function
slong acb_siegel_kappa(acb_t sqrtdet, const fmpz_mat_t mat, const acb_mat_t tau, int sqr, slong prec)
slong acb_siegel_kappa2(const fmpz_mat_t mat)

If `sqr` is zero, `acb_siegel_kappa` returns $0 \le r < 8$ with
$\kappa(\mathit{mat}) = \zeta_8^r$ and sets `sqrtdet` to the matching
square root of $\det(\gamma\tau + \delta)$ in the theta transformation
formula; if `sqr` is nonzero, it returns $0 \le r < 4$ with
$\kappa(\mathit{mat})^2 = i^r$ and sets `sqrtdet` to
$\det(\gamma\tau + \delta)$. `acb_siegel_kappa2` returns $0 \le r < 4$
with $\kappa(\mathit{mat})^2 = i^r$, which does not depend on a choice of
square root.

::: details
By [@Igu1972, p. 176] and [@Mum1983, p. 189], for any symplectic $m$,
$(z,\tau)\in \mathbb{C}^g\times \mathcal{H}_g$ and characteristic $(a,b)$,

$$
\theta_{a,b}(m\cdot(z,\tau)) = \kappa(m) \zeta_8^{e(m, a, b)} \det(\gamma\tau + \delta)^{1/2} e^{\pi i z^T (\gamma\tau + \delta)^{-1} \gamma z} \theta_{a',b'}(z,\tau),
$$

where $\zeta_8 = e^{i\pi/4}$, $(a', b')$ and the integer $e(m, a, b)$
(Igusa's $\phi_m$) are given by `acb_theta_char_table`, and the 8th root
of unity $\kappa(m)$ is determined up to sign by the choice of branch of
$\det(\gamma\tau+\delta)^{1/2}$.

The function decomposes `mat` with `sp2gz_decompose`. For trigonal and
block-diagonal matrices the action on $\theta_{0,0}$ is explicit; for
embedded matrices from $\mathrm{SL}_2(\Z)$ it uses
`acb_modular_theta_transform`; and for an embedded $J$ from dimension
$2 \le r \le g$, $\kappa(m) \zeta_8^{e(m,0,0)} i^{r/2} \det(\tau_0)^{1/2} = 1$,
where $\tau_0$ is the upper left $r \times r$ block of $\tau$ and the
square root is the branch that equals $i^{r/2}\det(Y_0)^{1/2}$ for purely
imaginary $\tau_0 = iY_0$. To find this branch (when `sqr` is zero), it
picks a purely imaginary matrix $A$ such that
$P(t) = \det(A + \tfrac{t+1}{2} (\tau_0 - A))$ has degree $r$ and no roots
near the segment $[-1,1]$, and follows the branch of $P(t)^{1/2}$ from
$t = -1$ to $t = 1$ as in [@MN2019].
:::
:::

::: function
void acb_siegel_randtest(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits)
void acb_siegel_randtest_reduced(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits)
void acb_siegel_randtest_compact(acb_mat_t tau, flint_rand_t state, int exact, slong prec)

Set `tau` to a random matrix in $\mathcal{H}_g$: possibly far from
reduced; reduced, with an imaginary part that may have large entries; or
reduced with bounded imaginary part. With `exact` nonzero, the entries
are exact dyadic numbers.
:::

::: function
void acb_siegel_randtest_vec(acb_ptr z, flint_rand_t state, slong g, slong prec)
void acb_siegel_randtest_vec_reduced(acb_ptr zs, flint_rand_t state, slong nb, const acb_mat_t tau, int exact, slong prec)

Set `z` to a random vector of length $g$, possibly with large entries,
respectively `zs` to `nb` random vectors in $[-1,1]^g + \tau[-1,1]^g$,
which are close to reduced with respect to $\tau$; with `exact` nonzero,
the entries are exact dyadic numbers.
:::

## Theta characteristics {#characteristics}

::: function
ulong acb_theta_char_set_slong_vec(const slong * vec, slong len)
int acb_theta_char_bit(ulong ch, slong j, slong n)

`acb_theta_char_set_slong_vec` returns the integer whose binary digits,
most significant first, are the `len` entries (0 or 1) of `vec`; with
`len` $= 2g$ and `vec` $= (a, b)$ this is the index of $\theta_{a,b}$ in
vectors of theta values. `acb_theta_char_bit` returns bit $j$ (counted
from 0, most significant first) of `ch` seen as an element of
$\{0,1\}^n$.
:::

::: function
void acb_theta_char_get_arb(arb_ptr v, ulong a, slong g)
void acb_theta_char_get_acb(acb_ptr v, ulong a, slong g)

Set `v` to $a/2 \in \mathbb{R}^g$ or $\mathbb{C}^g$.
:::

::: function
slong acb_theta_char_dot(ulong a, ulong b, slong g)
slong acb_theta_char_dot_slong(ulong a, const slong * n, slong g)

Return $\sum_{i} a_i b_i$, respectively $\sum_i a_i n_i$, modulo 4, as an
integer in $[0, 3]$.
:::

::: function
int acb_theta_char_is_even(ulong ab, slong g)

Returns whether the characteristic $(a,b)$ is even, that is, $a^T b$ is
even. For odd characteristics, $\theta_{a,b}(0, \tau)$ vanishes
identically.
:::

::: function
void acb_theta_char_table(ulong * ch, slong * e, const fmpz_mat_t mat, ulong ab, int all)

Sets `ch` to the characteristic $(a', b')$ and `e` to $e(\mathit{mat}, a, b)$
of the transformation formula (see `acb_siegel_kappa`) for the
characteristic `ab`, or if `all` is nonzero, sets the vectors `ch` and
`e` of length $2^{2g}$ to these data for all characteristics (`ab` is
then ignored).
:::

::: function
void acb_theta_char_shuffle(acb_ptr res, const fmpz_mat_t mat, acb_srcptr th, int sqr, slong prec)

Applies the transformation formula for `mat` to the vector `th` of
$2^{2g}$ theta values, omitting the factors $\kappa$, the determinant and
the exponential, and stores the result in `res`. With `sqr` nonzero,
$\zeta_8$ is replaced by $i$, as for squared theta values. Used for
testing.
:::

## Derivatives {#jets}

Derivatives with respect to $\tau$ follow from those with respect to $z$
by the heat equation

$$
\frac{\partial\theta_{a,b}}{\partial \tau_{j,k}} = \frac{1}{2\pi i(1 +\delta_{j,k})}
\frac{\partial^2\theta_{a,b}}{\partial z_j \partial z_k},
$$

so only partial derivatives with respect to $z$ are computed. A
*derivation tuple* $(k_0, \ldots, k_{g-1})$ (a vector of `slong` of
length $g$) stands for the operator
$\frac{1}{k_0! \cdots k_{g-1}!} \frac{\partial^{|k|}}{\partial z_0^{k_0}\cdots \partial z_{g-1}^{k_{g-1}}}$,
as in a Taylor expansion, where $|k| = \sum k_i$. All tuples of total
order at most `ord` are considered, ordered by total order and then
reverse-lexicographically: for $g = 2$, $(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), \ldots$.
These functions work with such vectors of derivatives of any holomorphic
functions.

::: function
slong acb_theta_jet_nb(slong ord, slong g)

Returns the number $\binom{\mathit{ord} + g}{g}$ of derivation tuples of
total order at most `ord` (0 if `ord` is negative).
:::

::: function
slong acb_theta_jet_total_order(const slong * tup, slong g)
void acb_theta_jet_tuples(slong * tups, slong ord, slong g)
slong acb_theta_jet_index(const slong * tup, slong g)

Return the total order $|k|$ of `tup`; set `tups` to the concatenation
of all derivation tuples of total order at most `ord`, in order; and
return the index of `tup` in that order.
:::

::: function
void acb_theta_jet_mul(acb_ptr res, acb_srcptr v1, acb_srcptr v2, slong ord, slong g, slong prec)

Sets `res` to the derivatives of $fh$ given the derivatives `v1` and `v2`
of $f$ and $h$, up to total order `ord`.
:::

::: function
void acb_theta_jet_compose(acb_ptr res, acb_srcptr v, const acb_mat_t N, slong ord, slong prec)

Sets `res` to the derivatives of $z \mapsto f(Nz)$, given the derivatives
`v` of $f$ at $Nz$; $g$ is the size of $N$.
:::

::: function
void acb_theta_jet_exp_pi_i(acb_ptr res, arb_srcptr a, slong ord, slong g, slong prec)
void acb_theta_jet_exp_qf(acb_ptr res, acb_srcptr z, const acb_mat_t N, slong ord, slong prec)

Set `res` to the derivatives at $z = 0$ of $\exp(\pi i (a_0 z_0 + \cdots + a_{g-1} z_{g-1}))$,
respectively to the derivatives at the given $z$ of $\exp(\pi i z^T N z)$.
:::

## Ellipsoids {#ellipsoids}

The summation algorithms sum the theta series over the lattice points in
an ellipsoid, as in [@DHBHS2004]; see
[](algorithms/special-functions/theta-functions#summation). Given an
upper triangular matrix $C$ with positive diagonal (a *Cholesky matrix*),
a radius $R \ge 0$, a center $v \in \mathbb{R}^g$ and $1 \le d \le g$, the
ellipsoid $E$ is the set of $n \in \Z^g$ with
$(v + Cn)^T (v + Cn) \le R^2$ and fixed last coordinates
$n_d, \ldots, n_{g-1}$. It is stored recursively: the range and midpoint
of the allowed values of $n_{d-1}$, and for each value a
$(d-1)$-dimensional child, split into left and right children by the
midpoint. For $d = g$ and fixed $C$, this takes $O(R^{g-1})$ space for
about $O(R^g)$ points.

::: type
acb_theta_eld_struct
acb_theta_eld_t

An ellipsoid as described above. An `acb_theta_eld_t` is an array of
length 1, so that it is passed by reference.
:::

::: function
void acb_theta_eld_init(acb_theta_eld_t E, slong d, slong g)
void acb_theta_eld_clear(acb_theta_eld_t E)

Initialise `E` as a $d$-dimensional ellipsoid in ambient dimension $g$,
with $1 \le d \le g$, and free it with all its children.
:::

::: function
int acb_theta_eld_set(acb_theta_eld_t E, const arb_mat_t C, const arf_t R2, arb_srcptr v)

Sets `E` to the ellipsoid with Cholesky matrix $C$, squared radius
$R^2$ = `R2` and center $v$, and returns 1. Returns 0 (leaving `E`
undefined) if the coordinates of its points do not fit in `slong` or the
ellipsoid is unreasonably large.
:::

::: function
slong acb_theta_eld_nb_pts(const acb_theta_eld_t E)
void acb_theta_eld_points(slong * pts, const acb_theta_eld_t E)

Return the number of points of `E`, and set `pts` (allocated by the
caller) to the concatenation of the points, as vectors of length $g$.
:::

::: function
slong acb_theta_eld_box(const acb_theta_eld_t E, slong j)

Returns an upper bound for $|n_j|$ over the points $n$ of `E`, for
$0 \le j < g$.
:::

::: function
slong acb_theta_eld_nb_border(const acb_theta_eld_t E)
void acb_theta_eld_border(slong * pts, const acb_theta_eld_t E)

Return the number of points in the *border* of `E`, a set of points just
outside it, and set `pts` (allocated by the caller) to these points.
Used for testing.
:::

::: function
int acb_theta_eld_contains(const acb_theta_eld_t E, const slong * pt)

Returns whether the point `pt` of length $g$ is in `E`.
:::

::: function
void acb_theta_eld_print(const acb_theta_eld_t E)

Prints a complete description of `E`, which can be long in high
dimension.
:::

::: function
void acb_theta_eld_distances(arb_ptr ds, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)

Sets `ds` to the concatenation of `nb` vectors of length $2^g$: for each
$z$, the squared distances $\mathrm{Dist}_\tau(-Y^{-1}y, \Z^g + \tfrac a2)^2$
for all $a \in \{0,1\}^g$, for the norm $\lVert x \rVert_\tau^2 = \pi x^T Y x$.

::: details
The coordinates of $-Y^{-1}y$ are rounded to a point of
$\Z^g + \tfrac a2$, which gives an upper bound for the distance, and the
points of the ellipsoid of that radius are enumerated to find any closer
point.
:::
:::

## Summation algorithms {#summation}

These functions evaluate theta functions by summing the series over an
ellipsoid and adding a bound for the tail. They are used at low and
moderate precision, and inside the quasi-linear algorithm. The inputs
are assumed to be reduced, which allows one ellipsoid to serve several
vectors $z$. For $g = 1$ they call `acb_modular_theta_sum`.

::: function
void acb_theta_sum_radius(arf_t R2, arf_t eps, const arb_mat_t cho, slong ord, slong prec)

Sets `eps` to $2^{-\mathit{prec}}$ and `R2` to a squared radius
$R^2 \ge \max(4, 2\,\mathit{ord})$ such that
$(1 + \sqrt{8/\pi}) R^{g - 1 + 2\,\mathit{ord}} e^{-R^2} \prod_{j} (1 + \sqrt{2\pi}/\gamma_j) \le 2^{-\mathit{prec}}$,
where $\gamma_j$ are the diagonal entries of $C$ = `cho`. For `ord` = 0,
summing over the ellipsoid of radius $R$ computes
$\widetilde{\theta}_{a,b}$ with a tail of absolute value at most `eps`.
:::

::: function
void acb_theta_sum_jet_radius(arf_t R2, arf_t eps, const arb_mat_t cho, arb_srcptr v, slong ord, slong prec)

Sets `R2` and `eps` as `acb_theta_sum_radius` does, but for the
derivatives up to order `ord`, where $v = -C Y^{-1} y$.

::: details
The differentiated series is

$$
\frac{\partial^{|k|}\theta_{a,b}}{\partial z_0^{k_0}\cdots \partial z_{g-1}^{k_{g-1}}}(z,\tau)
= (2\pi i)^{|k|} e^{\pi y^T Y^{-1} y} \sum_{n\in \Z^g + \tfrac a2}
n_0^{k_0} \cdots n_{g-1}^{k_{g-1}} \xi_n e^{-\pi (n + Y^{-1}y)^T Y (n + Y^{-1}y)}
$$

with $|\xi_n| = 1$, and with $m = Cn + v$,
$|n_0^{k_0}\cdots n_{g-1}^{k_{g-1}}| \le (\lVert C^{-1}\rVert_\infty \lVert m\rVert_2 + \lVert Y^{-1}y\rVert_\infty)^{|k|}$.
The function first calls `acb_theta_sum_radius` with `ord` = 0. If
$R \le \lVert Y^{-1}y\rVert_\infty/\lVert C^{-1}\rVert_\infty$, it
multiplies `eps` by $\max(1, 2 \lVert Y^{-1}y \rVert_\infty)^{\mathit{ord}}$.
Otherwise it calls `acb_theta_sum_radius` with the given `ord`, takes the
maximum of `R2` and $(\lVert Y^{-1}y \rVert_\infty /\lVert C^{-1} \rVert_\infty)^2$,
and multiplies `eps` by $\max(1, 2\lVert C^{-1}\rVert_\infty)^{\mathit{ord}}$.
:::
:::

::: function
void acb_theta_sum_term(acb_t res, acb_srcptr z, const acb_mat_t tau, const slong * tup, const slong * n, slong prec)

Sets `res` to $n_0^{k_0} \cdots n_{g-1}^{k_{g-1}}\exp(\pi i(n^T \tau n + 2 n^T z))$,
where $k$ = `tup` (or the zero tuple if `tup` is `NULL`). Used for
testing.
:::

::: function
slong acb_theta_sum_addprec(const arb_t d)

Returns an integer close to $d/\log 2$ if $d$ is finite and of
reasonable size, and 0 otherwise: the number of extra bits needed for a
term of size $e^{-d}$.
:::

::: type
acb_theta_sum_worker_t

The type of the function that processes one line (a one-dimensional
slice) of an ellipsoid:
`void worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2, const slong * precs, slong len, const acb_t c, const slong * coords, slong ord, slong g, slong prec, slong fullprec)`.

::: details
For a line with fixed $n_1, \ldots, n_{g-1}$ and
$n_{\min} \le n_0 \le n_{\max}$, the arguments are: the output `th` to
which the terms are added; $v_1 = (e^{\pi i j^2 \tau_{0,0}})_j$ and
$v_2 = (x^j)_j$ for $n_{\min} \le j \le n_{\max}$ with
$x = e^{2\pi i z_0} \prod_{k \ge 1} e^{2\pi i n_k \tau_{0,k}}$; the
cofactor
$c = \prod_{k \ge 1} e^{2\pi i n_k z_k} \prod_{1\leq j\leq k < g} e^{\pi i (2 - \delta_{j,k}) n_j n_k \tau_{j,k}}$;
the working precision of each term in `precs`; the common length
`len` $= n_{\max} - n_{\min} + 1$ of `v1`, `v2` and `precs`;
`coords` $= (n_{\min}, n_1, \ldots, n_{g-1})$; the derivation order
`ord`; the precision `prec` of the line and `fullprec` for the sum.
This lets the worker use `acb_dot`, with fewer than two complex
multiplications per lattice point on average and with lower precision
for points far from the center.
:::
:::

::: function
void acb_theta_sum_sqr_pow(acb_ptr * sqr_pow, const acb_mat_t exp_tau, const acb_theta_eld_t E, slong prec)

For $0 \le j < g$, sets `sqr_pow[j]` (preallocated with length $B_j + 1$)
to the powers $e_j^{n^2}$, $0 \le n \le B_j$, where $e_j$ is the $j$-th
diagonal entry of `exp_tau` and $B_j$ = `acb_theta_eld_box(E, j)`.
:::

::: function
void acb_theta_sum_work(acb_ptr th, slong len, acb_srcptr exp_z, acb_srcptr exp_z_inv, const acb_mat_t exp_tau, const acb_mat_t exp_tau_inv, const acb_ptr * sqr_pow, const acb_theta_eld_t E, slong ord, slong prec, acb_theta_sum_worker_t worker)

Runs the summation over the ellipsoid `E` ($g \ge 2$), calling `worker`
on each line, without adding a tail bound.

::: details
The entries $(j, k)$, $j \le k$, of `exp_tau` and `exp_tau_inv` must
contain $\exp(\pi i (2 - \delta_{j,k}) \tau_{j,k})$ and its inverse,
`sqr_pow` the output of `acb_theta_sum_sqr_pow`, and `exp_z`,
`exp_z_inv` the values $\exp(2\pi i z_j)$ and their inverses for one
vector $z$. The parameters `len`, `ord` and the output `th` are passed
to the worker. These data are normally taken from contexts of type
`acb_theta_ctx_tau_t` and `acb_theta_ctx_z_t`.
:::
:::

::: function
void acb_theta_sum(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb, const acb_theta_ctx_tau_t ctx_tau, arb_srcptr distances, int all_a, int all_b, int tilde, slong prec)

Evaluates theta functions by summation at the `nb` pairs $(z, \tau)$
given by the contexts. With `all_a` (respectively `all_b`) zero, only
$a = 0$ (respectively $b = 0$) is computed, otherwise all values; with
`tilde` nonzero, the normalised $\widetilde{\theta}_{a,b}$ are computed.
The vector `distances` of length $2^g$ (common to all $z$) sets the
precision added for each $a$; it may be zero, or the output of
`acb_theta_eld_distances`, which is meaningful when the vectors $z$ differ
by real vectors.
:::

::: function
void acb_theta_sum_jet(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb, const acb_theta_ctx_tau_t ctx_tau, slong ord, int all_a, int all_b, slong prec)

Sets `th` to the derivatives up to total order `ord` of the theta
functions at the `nb` pairs given by the contexts, by summation, with
the characteristics selected by `all_a` and `all_b` as in
`acb_theta_sum`.
:::

## Context structures {#contexts}

The summation functions only need exponentials of the entries of $\tau$
and $z$. Context structures store them so that they can be reused for
several summations and updated cheaply when $(z, \tau)$ is replaced by
$(2z, 2\tau)$ in the quasi-linear algorithm.

::: type
acb_theta_ctx_tau_struct
acb_theta_ctx_tau_t

The data for summation at a matrix $\tau$: in particular the matrix
`exp_tau_div_4` with entries $\exp(\pi i (2 - \delta_{j,k}) \tau_{j,k}/4)$
for $j \le k$ (and its squares, fourth powers and inverses), and for
$g > 1$ the Cholesky matrix of $\pi Y$ and $Y^{-1}$. With `allow_shift`,
additional fields allow the evaluation of $\theta_{a,0}$ for $a \ne 0$.
An `acb_theta_ctx_tau_t` is an array of length 1.
:::

::: type
acb_theta_ctx_z_struct
acb_theta_ctx_z_t

The data for summation at a vector $z$ (given a $\tau$ context): the
values $\exp(2\pi i z_j)$, their squares and inverses, and for $g > 1$
the center $v$ of the ellipsoids. An `acb_theta_ctx_z_t` is an array of
length 1.
:::

::: function
void acb_theta_ctx_tau_init(acb_theta_ctx_tau_t ctx, int allow_shift, slong g)
void acb_theta_ctx_tau_clear(acb_theta_ctx_tau_t ctx)
void acb_theta_ctx_z_init(acb_theta_ctx_z_t ctx, slong g)
void acb_theta_ctx_z_clear(acb_theta_ctx_z_t ctx)

Initialise contexts for dimension $g$ (with the extra fields for
$\theta_{a,0}$, $a \ne 0$, if `allow_shift` is nonzero), and free them.
:::

::: function
acb_theta_ctx_z_struct * acb_theta_ctx_z_vec_init(slong nb, slong g)
void acb_theta_ctx_z_vec_clear(acb_theta_ctx_z_struct * vec, slong nb)

Allocate and initialise a vector of `nb` contexts for dimension $g$, and
free it.
:::

::: function
void acb_theta_ctx_tau_set(acb_theta_ctx_tau_t ctx, const acb_mat_t tau, slong prec)
void acb_theta_ctx_z_set(acb_theta_ctx_z_t ctx, acb_srcptr z, const acb_theta_ctx_tau_t ctx_tau, slong prec)

Compute the context data for $\tau$, respectively for $z$ given the
context of $\tau$. The dimensions must match.
:::

::: function
void acb_theta_ctx_tau_dupl(acb_theta_ctx_tau_t ctx, slong prec)
void acb_theta_ctx_z_dupl(acb_theta_ctx_z_t ctx, slong prec)

Replace the context data for $\tau$ by those for $2\tau$, respectively
those for $(z, \tau)$ by those for $(2z, 2\tau)$, by squaring, which is
much cheaper than recomputing them.
:::

::: function
void acb_theta_ctx_z_add_real(acb_theta_ctx_z_t res, const acb_theta_ctx_z_t ctx, const acb_theta_ctx_z_t ctx_real, slong prec)

Given contexts for $(z, \tau)$ and $(t, \tau)$ with $t$ real, sets `res`
to a context for $(z + t, \tau)$.
:::

::: function
void acb_theta_ctx_z_common_v(arb_ptr v, const acb_theta_ctx_z_struct * vec, slong nb, slong prec)

Sets `v` to an ellipsoid center valid for all `nb` contexts (the union
of their centers), for one summation over all of them.
:::

::: function
void acb_theta_ctx_exp_inv(acb_t exp_inv, const acb_t exp, const acb_t x, int is_real, slong prec)
void acb_theta_ctx_sqr_inv(acb_t sqr_inv, const acb_t inv, const acb_t sqr, int is_real, slong prec)

Given $e = \exp(\pi i x)$, `acb_theta_ctx_exp_inv` sets `exp_inv` to
$\exp(-\pi i x)$: by conjugation if `is_real` is nonzero, else by
inversion, recomputing from $x$ if the inverse is not finite. Given
$\exp(-\pi i x)$ and $\exp(2\pi i x)$, `acb_theta_ctx_sqr_inv` sets
`sqr_inv` to $\exp(-2\pi i x)$, by conjugation if `is_real` is nonzero
and by squaring otherwise.
:::

::: function
int acb_theta_ctx_tau_overlaps(const acb_theta_ctx_tau_t ctx1, const acb_theta_ctx_tau_t ctx2)
int acb_theta_ctx_z_overlaps(const acb_theta_ctx_z_t ctx1, const acb_theta_ctx_z_t ctx2)

Return whether the data of the two contexts overlap in the sense of
`acb_overlaps`. Used for testing.
:::

## Duplication steps {#agm}

The quasi-linear algorithm uses the duplication formula
$\theta_{a,0}(z,\tau) \theta_{a,0}(z',\tau) = \sum_{a'} \theta_{a',0}(z+z',2\tau) \theta_{a+a',0}(z-z',2\tau)$,
which also holds for the normalised functions; see
[](algorithms/special-functions/theta-quasi-linear#duplication).

::: function
void acb_theta_agm_sqrt(acb_ptr res, acb_srcptr a, acb_srcptr roots, slong nb, slong prec)

Sets each of the `nb` entries of `res` to the square root of the entry of
`a` that overlaps the corresponding entry of `roots` (a low-precision
approximation) and whose opposite does not. If both square roots overlap
`roots`, the entry is their union; if neither does, it is indeterminate.
:::

::: function
void acb_theta_agm_mul(acb_ptr res, acb_srcptr a1, acb_srcptr a2, slong g, int all, slong prec)

For $0 \le k < 2^g$, sets entry $k$ of `res` to
$\sum_{b \in \{0,1\}^g} a_{1,b}\, a_{2, b + k}$, where $b + k$ is the
bitwise exclusive or. If `all` is nonzero, `res` has length $2^{2g}$ and
also contains, for $1 \le a < 2^g$, the vectors with entries
$\sum_{b} (-1)^{a^T b} a_{1,b}\, a_{2, b+k}$. This is computed with two
Hadamard transforms and $2^g$ products [@LT2016], which loses precision
when the entries have different orders of magnitude; it is faster when
`a1` and `a2` are the same pointer (squarings).
:::

::: function
void acb_theta_agm_mul_tight(acb_ptr res, acb_srcptr a0, acb_srcptr a, arb_srcptr d0, arb_srcptr d, slong g, int all, slong prec)

The same as `acb_theta_agm_mul` on `a0` and `a`, but with error bounds
that are much tighter when $|a_k| \approx e^{-d_k}$ and
$|a_{0,k}| \approx e^{-d_{0,k}}$, where `d0` and `d` are the output of
`acb_theta_eld_distances` at $(0, \tau)$ and $(z, \tau)$, as for
normalised theta values.

::: details
The function computes $m$ and $\varepsilon$ such that
$|a_k| \le m e^{-d_k}$ and $\operatorname{rad}(a_k) \le \varepsilon e^{-d_k}$
for all $k$, and $m_0, \varepsilon_0$ for `a0` similarly, calls
`acb_theta_agm_mul` on the midpoints at a higher precision, and adds
$2^g e^{-d_k} (m_0 \varepsilon + m \varepsilon_0 + \varepsilon\varepsilon_0)$
to the radius of entry $k$. This is valid because for $v = -Y^{-1} y$,

$$
\mathrm{Dist}_\tau(0, \Z^g + \tfrac b2)^2 + \mathrm{Dist}_\tau(v, \Z^g + \tfrac{b + k}{2})^2 \ge \mathrm{Dist}_\tau(v, \Z^g + \tfrac{k}{2})^2:
$$

if $x \in \Z^g + \frac b2$ and $w \in v + \Z^g + \frac{b+k}{2}$ realise
the two distances, then so does $-x$, and choosing the sign with
$\langle x, w \rangle_\tau \le 0$ gives
$\lVert x + w \rVert_\tau^2 \le \lVert x \rVert_\tau^2 + \lVert w \rVert_\tau^2$
with $x + w \in v + \Z^g + \frac k2$.
:::
:::

## Quasi-linear algorithms on reduced input {#ql}

These functions implement the algorithm of [@EK2025] described in
[](algorithms/special-functions/theta-quasi-linear). The inputs are
assumed to be reduced.

::: function
int acb_theta_ql_nb_steps(slong * pattern, const acb_mat_t tau, int cst, slong prec)

Decides how many duplication steps to perform at precision `prec` and
when to switch to lower dimensions, and stores the decision in `pattern`
(of length $g$): roughly, entry $j$ is the number $m$ of steps such that
$2^m \gamma_j^2$ is of the order of `prec`, where $\gamma_j$ is the
$j$-th diagonal entry of the Cholesky matrix of $\pi Y$. The flag `cst`
should be nonzero when only theta constants ($z = 0$) are needed. Returns
0 if the Cholesky matrix is indeterminate or too extreme, and 1
otherwise.

::: details
After $m$ duplications the ellipsoid of the summation at $2^m\tau$ has
width $O(1)$ in direction $j$. For reduced $\tau$, `pattern` is expected
to be roughly decreasing. The rough values are adjusted by rules obtained
from timings with the profiling program `p-ql_exact` (in
`src/acb_theta/profile`): dimension-lowering is forced or avoided
depending on how many steps have already been done, a few steps are
removed in genus 1, and the result is made nonincreasing and
nonnegative. This function is the main place to tune the algorithm.
:::
:::

::: function
int acb_theta_ql_lower_dim(acb_ptr * new_zs, acb_ptr * cofactors, slong ** pts, slong * nb, arf_t err, slong * fullprec, acb_srcptr z, const acb_mat_t tau, arb_srcptr distances, slong s, ulong a, slong prec)

Implements the dimension-lowering step: writes the normalised theta
values $\widetilde{\theta}_{a',0}(z, \tau)$ whose last $g - s$ bits of
$a'$ equal `a` as short linear combinations of theta values in dimension
$s$ ($1 \le s < g$). The pair $(z,\tau)$ should be exact (and ideally
reduced), and `distances` its output of `acb_theta_eld_distances`.

::: details
The outputs are:

- `fullprec`: the precision at which the theta values in dimension $s$
  must be computed, `prec` plus guard bits from the largest relevant
  distance;
- `err`: the tail bound of `acb_theta_sum_radius` at precision
  `fullprec` for a squared radius $R^2$ (the ellipsoid itself is not
  computed);
- `pts`, `nb`: writing $a' = (a_0, a)$ and $n = (n_0, n_1)$ with
  $n_0 \in \Z^s + \tfrac{a_0}{2}$, $n_1 \in \Z^{g - s} + \tfrac{a}{2}$,
  the possible $n_1$ lie in an ellipsoid of squared radius $R^2$ in
  dimension $g - s$ whose Cholesky matrix is the lower right block of
  that of $\pi Y$; its points $n_1 - \tfrac a2$ are listed in `pts`
  (allocated by the function) and their number in `nb`;
- `new_zs`, `cofactors`: writing
  $\tau = \begin{pmatrix} \tau_0 & x\\x^T & \tau_1\end{pmatrix}$ and
  $z = (z_0, z_1)$, the terms of the series with a given $n_1$ sum to
  $e^{\pi i (n_1^T \tau_1 n_1 + 2 n_1^T z_1)} e^{\pi (y_0^T Y_0^{-1} y_0 - y^T Y^{-1} y)} \widetilde{\theta}_{a_0,0}(z_0 + x n_1, \tau_0)$;
  `new_zs` (allocated) receives the `nb` vectors $z_0 + x n_1$ (still
  exact) and `cofactors` (allocated) the exponential factors.

The caller frees `new_zs`, `cofactors` and `pts`. The return value is 0
if the computation fails (for instance if the ellipsoid cannot be
computed), and 1 otherwise.
:::
:::

::: function
void acb_theta_ql_recombine(acb_ptr th, acb_srcptr th0, acb_srcptr cofactors, const slong * pts, slong nb, const arf_t err, slong fullprec, slong s, ulong a, int all, slong g, slong prec)

The converse of `acb_theta_ql_lower_dim`: recovers the normalised theta
values at $(z, \tau)$ from the values in dimension $s$ at the vectors
`new_zs` and $\tau_0$. If `all` is nonzero, `th0` holds `nb` vectors of
the $2^{2s}$ values $\widetilde{\theta}_{a_0,b_0}$ and `th` receives all
$2^{2g}$ values $\widetilde{\theta}_{a,b}(z,\tau)$; otherwise `th0` holds
`nb` vectors of the $2^s$ values $\widetilde{\theta}_{a_0,0}$ and `th`
receives the $2^g$ values $\widetilde{\theta}_{a,0}(z,\tau)$. The other
arguments are the outputs of `acb_theta_ql_lower_dim`.
:::

::: function
int acb_theta_ql_setup(acb_ptr rts, acb_ptr rts_all, acb_ptr t, slong * guard, slong * easy_steps, acb_srcptr zs, slong nb, const acb_mat_t tau, arb_srcptr distances, slong nb_steps, int all, slong prec)

Prepares `nb_steps` $\ge 1$ duplication steps for the `nb` $\ge 1$ exact,
reduced pairs $(z, \tau)$ given by `zs`, whose first vector must be zero:
it computes low-precision approximations of all the square roots to be
taken, chooses an auxiliary real vector $t$ if needed, and returns 1 on
success. It returns 0 if no suitable $t$ is found, which should happen
only with negligible probability. Inexact input gives correct but
possibly useless results.

::: details
`distances` holds `nb` vectors of length $2^g$ from
`acb_theta_eld_distances`. With `all` nonzero, the values
$\widetilde{\theta}_{a,b}(z,\tau)$ for all characteristics are prepared,
otherwise only $\widetilde{\theta}_{a,0}$. The vectors `rts`, `rts_all`,
`t` and `easy_steps` must be allocated with lengths
$2^g \cdot 3 \cdot \mathit{nb} \cdot \mathit{nb\_steps}$,
$2^{2g} \cdot \mathit{nb}$ (used only if `all` is set), $g$ and `nb`.

Starting at a very low precision (such as 8 bits):

1. For each $z$, approximations of $\widetilde{\theta}_{a,b}(z,\tau)$
   (stored in `rts_all`) or $\widetilde{\theta}_{a,0}(z,\tau)$ (in `rts`)
   are computed by summation, with error bounds adjusted by the
   distances so that they likely exclude zero. If none contains zero, the
   last step can use the simplest duplication formula, and approximations
   of $\widetilde{\theta}_{a,0}(2^k z, 2^k\tau)$ are computed for
   $k = 1, 2, \ldots$ until `nb_steps` $- 1$ or until one contains zero;
   the number of such *easy steps* is stored in `easy_steps`.
2. If all entries of `easy_steps` equal `nb_steps`, this is done.
   Otherwise a random real vector $t$ is chosen, and the first entry of
   `easy_steps` (for $z = 0$) is set to the minimum of all entries.
3. For each $z$ with $m$ easy steps, approximations of
   $\widetilde{\theta}_{a,0}(2^k (z + t), 2^k\tau)$ and
   $\widetilde{\theta}_{a,0}(2^k(z + 2t), 2^k\tau)$ are computed for
   $m \le k <$ `nb_steps` (only the second for $k = 0$, or all
   $\widetilde{\theta}_{a,b}(z + 2t, \tau)$ if `all` is set). If one
   contains zero, another $t$ is tried, a few times.
4. If no suitable $t$ was found, the precision is doubled and the process
   restarts from step 1, up to precision `prec`.

Finally `guard` is set to the number of bits expected to be lost in the
duplication steps.
:::
:::

::: function
void acb_theta_ql_exact(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, const slong * pattern, int all, int shifted_prec, slong prec)

Runs the quasi-linear algorithm for the `nb` exact, reduced pairs
$(z, \tau)$, with the duplication and dimension-lowering pattern
`pattern`, and sets `th` to the normalised values
$\widetilde{\theta}_{a,b}(z,\tau)$ for all characteristics (`all`
nonzero) or $\widetilde{\theta}_{a,0}(z,\tau)$. With `shifted_prec`
nonzero, the precision of each value takes the distances into account,
as in `acb_theta_sum`. Inexact input gives correct but possibly useless
results.

::: details
With no duplication steps, `acb_theta_sum` is called. Otherwise
`acb_theta_ql_setup` prepares the steps; the theta values at $2^m \tau$
are computed by summation, or with `acb_theta_ql_lower_dim`, a recursive
call in dimension $s$ on a longer list of vectors, and
`acb_theta_ql_recombine`; and the duplication formulas with square roots
are applied $m$ times.
:::
:::

::: function
void acb_theta_ql_local_bound(arb_t c, arb_t rho, acb_srcptr z, const acb_mat_t tau, slong ord)

Sets $c \ge 1$ and $0 < \rho \le 1$ such that $|\theta_{a,b}(\cdot,\tau)| \le c$
for all characteristics on every ball of radius $\rho$ (for the
$\infty$-norm) centered at a point of `z`, with $\rho$ chosen to make the
Cauchy bounds for derivatives up to order `ord` small.

::: details
One can take $c = c_0\exp((c_1 + c_2\rho)^2)$ for any $\rho$, with

$$
c_0 = 2^g \prod_{j=0}^{g-1} (1 + 2\gamma_j^{-1}), \qquad
c_1 = \sqrt{\pi y^T Y^{-1} y}, \qquad
c_2 = \sup_{\lVert x \rVert_\infty\leq 1} \sqrt{\pi x^T Y^{-1} x},
$$

where an upper bound for $c_2$ comes from the Cholesky decomposition of
$\pi Y^{-1}$. Then $\rho$ minimises $\exp((c_1 + c_2\rho)^2)/\rho^{2\,\mathit{ord}+1}$:
it is the minimum of 1 and the positive root of
$2c_2\rho (c_1 + c_2\rho) = 2\,\mathit{ord}+1$.
:::
:::

::: function
void acb_theta_ql_jet_error(arb_ptr err, acb_srcptr z, const acb_mat_t tau, acb_srcptr dth, slong ord, slong prec)

Given approximations `dth` of the derivatives of some $\theta_{a,b}$ up
to total order $\mathit{ord} + 2$ at $(z, \tau)$, sets `err` to a vector
such that the derivatives up to order `ord` at the midpoint of
$(z,\tau)$ and at any point of the balls differ by at most `err`
entrywise. This uses the heat equation and a Lipschitz-type bound.
:::

::: function
void acb_theta_ql_jet_fd(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong ord, int all, slong prec)
void acb_theta_ql_jet(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong ord, int all, slong prec)

Set `th` to the derivatives up to order `ord` of $\theta_{a,b}$ at the
`nb` pairs $(z, \tau)$, for all $a$ and for all $b$ (`all` nonzero) or
$b = 0$; the output is the concatenation of `nb` times $2^{2g}$ (or $2^g$)
vectors of derivatives. `acb_theta_ql_jet_fd` uses finite differences on
the output of `acb_theta_ql_exact`, in time quasi-linear in
$\mathit{prec} \cdot \mathit{ord}^g$ for fixed $g$; `acb_theta_ql_jet`
uses direct summation instead when `ord` is positive and the precision is
low (few duplication steps would be needed).

::: details
The finite difference method evaluates $\theta_{a,b}(z + h_n, \tau)$ at
$h_n = (\varepsilon \zeta^{n_0},\ldots, \varepsilon \zeta^{n_{g-1}})$ for
a primitive $m$-th root of unity $\zeta$, $m = \mathit{ord} + 1$, and all
$n \in \{0,\ldots, m - 1\}^g$; a discrete Fourier transform gives
$m^g \sum_{k \equiv p \bmod m} a_k \varepsilon^{|k|}$ for each $p$, where
$a_k$ are the Taylor coefficients. By the Cauchy formula with the bound of
`acb_theta_ql_local_bound`, the error is at most
$2cg\varepsilon^{|p|+m}/\rho^{m}$ before division by $\varepsilon^{|p|}$.
The function chooses $\varepsilon$, strips $(z, \tau)$ of their radii,
raises the precision to compensate for the division by
$\varepsilon^{\mathit{ord}} (\mathit{ord}+1)^g$, and restores certified
radii with `acb_theta_ql_jet_error` and derivatives to order
$\mathit{ord} + 2$ computed at low precision; see
[](algorithms/special-functions/theta-quasi-linear#derivatives).
:::
:::

## Siegel modular forms in genus 2 {#g2}

For $g = 2$, many Siegel modular forms can be computed from theta values,
in analogy with `acb_modular_delta` and `acb_modular_eisenstein` for
$g = 1$. A Siegel modular form of weight $\det^k\otimes \mathrm{Sym}^j$
is a holomorphic function $f : \mathcal{H}_2 \to \mathbb{C}_j[X]$ (polynomials of
degree at most $j$) with

$$
f((\alpha\tau + \beta)(\gamma\tau + \delta)^{-1}) = \det(\gamma\tau + \delta)^k\cdot \mathrm{Sym}^j(\gamma\tau + \delta)(f(\tau))
$$

for all $m \in \mathrm{Sp}_4(\Z)$, where for
$r = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$,
$\mathrm{Sym}^j(r)$ maps $P(X)$ to $(bX + d)^j P\bigl(\tfrac{aX + c}{bX + d}\bigr)$.
A nonzero $f$ exists only for even $j$. If poles along the diagonal
$\{\operatorname{diag}(\tau_1, \tau_2)\}$ are allowed, the ring of such
forms is finitely generated, and corresponds to the ring of covariants
of binary sextics: a covariant of degree $k$ and index $j$ corresponds to
a Siegel modular function of weight $\det^{k - j/2}\otimes \mathrm{Sym}^j$
[@CFG2017]. Writing
$\tau = \begin{pmatrix} \tau_1 & \tau_2 \\ \tau_2 & \tau_3\end{pmatrix}$,
$q_j = e^{2\pi i \tau_j}$.

::: function
void acb_theta_g2_detk_symj(acb_poly_t res, const acb_mat_t m, const acb_poly_t f, slong k, slong j, slong prec)

Sets `res` to $\det(m)^k \mathrm{Sym}^j(m)(f)$ for a $2 \times 2$ matrix
$m$ and a polynomial $f$ of degree at most $j$ (higher coefficients are
ignored).
:::

::: function
void acb_theta_g2_transvectant(acb_poly_t res, const acb_poly_t g, const acb_poly_t h, slong m, slong n, slong k, int lead, slong prec)

Sets `res` to the $k$-th transvectant of the polynomials $g$ and $h$,
seen as binary forms of degrees $m$ and $n$ (higher coefficients are
ignored), or to its leading coefficient (as a constant polynomial) if
`lead` is nonzero.

::: details
With $g$, $h$ homogeneous in $x_1, x_2$,

$$
(g,h)_k = \frac{(m-k)!(n-k)!}{m!n!}  \sum_{j=0}^{k} (-1)^{k-j} \binom{k}{j}
\frac{\partial^k g}{\partial x_1^{k-j}\partial x_2^j}
\frac{\partial^k h}{\partial x_1^{j}\partial x_2^{k-j}}.
$$
:::
:::

::: function
slong acb_theta_g2_character(const fmpz_mat_t mat)

Returns the value (0 or 1, in $\Z/2\Z$) at `mat` of the unique nontrivial
character of $\mathrm{Sp}_4(\Z)$ [@CFG2019, §12].
:::

::: function
void acb_theta_g2_even_weight(acb_t psi4, acb_t psi6, acb_t chi10, acb_t chi12, acb_srcptr th2, slong prec)

Sets the outputs to the Eisenstein series $\psi_4$, $\psi_6$ and the cusp
forms $\chi_{10}$, $\chi_{12}$, given the vector `th2` of the 16 squared
theta constants $\theta_{a,b}(0, \tau)^2$.

::: details
The formulas are those of [@Str2014, §7.1], normalised as
$\psi_4 = h_4/4$, $\psi_6 = h_6/4$, $\chi_{10} = -2^{-12} h_{10}$ and
$\chi_{12} = 2^{-15}h_{12}$, so that the Fourier expansions begin

$$
\begin{aligned}
\psi_4(\tau) &= 1 + 240(q_1 + q_3) + \cdots, &
\psi_6(\tau) &= 1 - 504(q_1 + q_3) + \cdots, \\
\chi_{10}(\tau) &= (q_2 - 2 + q_2^{-1}) q_1 q_3 + \cdots, &
\chi_{12}(\tau) &= (q_2 + 10 + q_2^{-1}) q_1 q_3 + \cdots.
\end{aligned}
$$

$\chi_{10}$ and $\chi_{12}$ differ from Igusa's [@Igu1979] by constant
factors.
:::
:::

::: function
void acb_theta_g2_chi5(acb_t res, acb_srcptr th, slong prec)

Sets `res` to $\chi_5 = - 2^{-6} \prod_{(a,b)\text{ even}} \theta_{a,b}(0, \tau)$,
given the 16 theta constants `th`. $\chi_5$ is a Siegel cusp form with
character [@CFG2019], and $\chi_5^2 = \chi_{10}$.
:::

::: function
void acb_theta_g2_chi35(acb_t res, acb_srcptr th, slong prec)

Sets `res` to the cusp form $\chi_{35}$, the unique scalar Siegel modular
form of weight 35 up to scalars, normalised as
$\chi_{35}(\tau) = q_1^2 q_3^2 (q_1 - q_3 )(q_2 - q_2^{-1}) + \cdots$,
given the 16 theta constants `th`, by the formula of Bolza [@Bol1887]
(see [@Mum1984, Prop. 6.2, p. 98] for Bolza's notation in terms of
characteristics).
:::

::: function
void acb_theta_g2_chi3_6(acb_poly_t res, acb_srcptr dth, slong prec)

Sets `res` to the vector-valued cusp form $\chi_{3,6}$ of weight
$\det^3\otimes \mathrm{Sym}^6$ with character, given the first
derivatives `dth` of all 16 theta functions at $z = 0$ (as computed by
`acb_theta_jet` with `ord` = 1 and `all` set) [@CFG2017]:
$\chi_{3,6}(\tau) = \frac{1}{64\pi^6} \prod_{(a,b) \text{ odd}} \left(\frac{\partial \theta_{a,b}}{\partial z_1}(0,\tau) x_1 + \frac{\partial\theta_{a,b}}{\partial z_2}(0,\tau) x_2\right)$.
:::

::: function
void acb_theta_g2_sextic_chi5(acb_poly_t res, acb_t chi5, const acb_mat_t tau, slong prec)

Sets `chi5` to $\chi_5(\tau)$ and `res` to $\chi_{-2,6}(\tau) = \chi_{3,6}/\chi_5$,
the meromorphic form that corresponds to the binary sextic itself under
the correspondence with covariants. The matrix $\tau$ is reduced, the
theta values are computed with `acb_theta_jet_notransform`, and the
transformation law of Siegel modular forms is applied.
:::

::: function
void acb_theta_g2_covariants(acb_poly_struct * res, const acb_poly_t f, int lead, slong prec)

Sets the 26 entries of `res` to the generators of the ring of covariants
of binary sextics evaluated at $f$ (terms of degree above 6 are
ignored), or to their leading coefficients if `lead` is nonzero.

::: details
The order and normalisation are:
$C_{1,6}=f$,
$C_{2,0}= 60(f,f)_6$,
$C_{2,4}= 75(f,f)_4$,
$C_{2,8}= 90(f,f)_2$,
$C_{3,2}= 30(f,C_{2,4})_4$,
$C_{3,6}= 30(f,C_{2,4})_2$,
$C_{3,8}= 6(f,C_{2,4})_1$,
$C_{3,12}= 6 (f,C_{2,8})_1$,
$C_{4,0}= 2 (C_{2,4},C_{2,4})_4$,
$C_{4,4}= 30 (f,C_{3,2})_2$,
$C_{4,6}= 6(f,C_{3,2})_1$,
$C_{4,10}= 2(C_{2,8},C_{2,4})_1$,
$C_{5,2}=(C_{2,4},C_{3,2})_2$,
$C_{5,4}=\frac 25 (C_{2,4},C_{3,2})_1$,
$C_{5,8}= 2(C_{2,8},C_{3,2})_1$,
$C_{6,0}= 2(C_{3,2},C_{3,2})_2$,
$C_{6,6}^{(1)}= \frac 25(C_{3,6},C_{3,2})_1$,
$C_{6,6}^{(2)}= \frac 83(C_{3,8},C_{3,2})_2$,
$C_{7,2}= 30(f,C_{3,2}^2)_4$,
$C_{7,4}= 12(f,C_{3,2}^2)_3$,
$C_{8,2}= \frac 25(C_{2,4},C_{3,2}^2)_3$,
$C_{9,4}= 4(C_{3,8},C_{3,2}^2)_4$,
$C_{10,0}= 20(f,C_{3,2}^3)_6$,
$C_{10,2}= \frac 65(f,C_{3,2}^3)_5$,
$C_{12,2}= \frac 85(C_{3,8},C_{3,2}^3)_6$,
$C_{15,0}= \frac{1}{30000} (C_{3,8},C_{3,2}^4)_8$.
The factors make the covariants of the generic sextic
$f = \sum a_i x_1^{6-i}x_2^i$ integral and primitive as polynomials in
$a_0,\ldots,a_6,x_1,x_2$.
:::
:::
