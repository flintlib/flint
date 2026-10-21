---
title: Arithmetic in number fields
summary: How nf_elem represents elements of Q[x]/(f) with special cases for degrees 1 and 2 and for Q(i), multiplies with reduction by precomputed powers of the generator, inverts by the extended Euclidean algorithm, and computes norms as resultants and traces from power sums of the roots.
complexity: Multiplication in degree $d$: one product of degree-$(d-1)$ integer polynomials plus $O(d^2)$ coefficient operations for the reduction; inverse: one extended GCD over $\mathbf{Q}$; norm: one resultant; trace: $O(d)$ operations.
order: 70
---

A number field $K = \mathbf{Q}(\theta)$ of degree $d$ is represented as
$\mathbf{Q}[x]/(f)$ for the minimal polynomial $f \in \mathbf{Q}[x]$ of
$\theta$, and its elements as polynomials of degree less than $d$. The
modules `nf` and `nf_elem` originate in the ANTIC library [@Har2015]; the
standard reference for the mathematics and the general algorithms is
[@Coh1996, Ch. 4].

## Representation {#representation}

An element $\alpha = (a_0 + a_1\theta + \cdots + a_{d-1}\theta^{d-1})/q$
is stored as an integer numerator polynomial and one positive
denominator $q$, in lowest terms (an `fmpq_poly_t`). A single
denominator keeps additions and multiplications in $\Z[x]$, where FLINT's
fast polynomial arithmetic applies, at the cost of a content computation
for canonicalisation.

Small degrees get special representations, selected by `nf_init`:

- for $d = 1$ the element is a rational number, a numerator and a
  denominator (`lnf_elem_t`);
- for $d = 2$ the numerator is a fixed array of three integers
  (`qnf_elem_t`), so that no memory management is needed and a product
  $(a_0 + a_1\theta)(b_0 + b_1\theta)$ can be stored before reduction;
- for $f = x^2 + 1$ ($K = \mathbf{Q}(i)$) a dedicated multiplication is
  used.

## Multiplication and reduction {#mul}

To multiply $\alpha = A(\theta)/q$ and $\beta = B(\theta)/r$, `nf_elem_mul`
computes $C = AB \in \Z[x]$ of degree at most $2d - 2$ with
`_fmpz_poly_mul`, the denominator $qr$, and then reduces $C$ modulo $f$
and removes the common content of numerator and denominator.

For $3 \le d < 30$ the reduction uses the precomputed remainders
$P_k = x^k \bmod f$ for $d \le k \le 2d - 2$, stored in the `nf_t`.

::: algorithm Reduction with precomputed powers {#alg-rem-powers}
**Input:** $C = \sum_{k \le 2d - 2} c_k x^k$, the remainders $P_k$.

**Output:** $C \bmod f$.

1. For $k = d, \ldots, \deg C$: add $c_k P_k$ to the low part
   $\sum_{k < d} c_k x^k$.
2. Return the low part.
:::

The algorithm is correct since $C \equiv \sum_{k < d} c_k x^k +
\sum_{k \ge d} c_k P_k \pmod f$. It costs $(d - 1)^2$ multiplications of
coefficients, like division with remainder, but has no dependency
between the steps: the remainders do not depend on $C$. If $f$ is monic
with integer coefficients, the $P_k$ are integral and the reduction
stays in $\Z[x]$ (`_fmpz_poly_rem_powers_precomp`); otherwise the $P_k$
have denominators that are powers of the leading coefficient of $f$ and
the rational version `_fmpq_poly_rem_powers_precomp` is used. For
$d \ge 30$ the product is reduced by division with remainder,
`_fmpz_poly_divrem` for monic integral $f$ and `_fmpq_poly_rem` with a
precomputed inverse of the leading coefficient otherwise.

For $d = 2$, a single step suffices: with $f = f_2x^2 + f_1x + f_0$,
$c_2\theta^2 = -c_2(f_1\theta + f_0)/f_2$, and for monic $f$ this is two
multiply-subtracts. For $\mathbf{Q}(i)$ the product
$(a + bi)(c + di) = (ac - bd) + (ad + bc)i$ is formed with two fused
multiply-add operations (`fmpz_fmms`, `fmpz_fmma`), and the denominator
is reduced by a three-term GCD.

`nf_elem_mul_red` with `red` $= 0$ skips the reduction, so that a sum
$\sum_j \alpha_j\beta_j$ can be accumulated in unreduced form and
reduced once at the end, which saves $d - 1$ reductions per term.

## Inversion and division {#inverse}

Since $f$ is irreducible, any nonzero $A$ of degree less than $d$ is
coprime to $f$, and the extended Euclidean algorithm gives $S, T \in
\mathbf{Q}[x]$ with $SA + Tf = 1$, so $A(\theta)^{-1} = S(\theta)$.
`nf_elem_inv` uses `fmpq_poly_xgcd` (and `_fmpq_poly_xgcd` for
$d = 2$), which takes out contents and denominators and calls the
multimodular `_fmpz_poly_xgcd_modular`; for coprime integer polynomials
that function returns the resultant $r$ and integer cofactors with
$sA + tf = r$ ([](algorithms/polynomials/gcd-integers)). For $d = 1$ the
inverse is
the reciprocal of a rational number. Division is multiplication by the
inverse. If $f$ is reducible, which `nf_init` does not exclude, an
element sharing a factor with $f$ has no inverse.

## Norm {#norm}

::: proposition {#prop-norm-resultant}
Let $f = c\prod_{j=1}^d (x - \theta_j)$ and let $A \in \mathbf{Q}[x]$ be
nonzero of degree $m < d$. Then
$N_{K/\mathbf{Q}}(A(\theta)) = \operatorname{Res}(f, A) / c^m$.
:::

::: proof
The norm is the product of the conjugates $A(\theta_j)$, and by the
product formula for resultants $\operatorname{Res}(f, A) =
c^m \prod_j A(\theta_j)$.
:::

`_nf_elem_norm` computes the resultant of the integer numerator of $f$
and the numerator of $\alpha$ with `_fmpq_poly_resultant`, which takes
out the contents and denominators and calls `_fmpz_poly_resultant`
(multimodular for large inputs), and divides by $c^m$ for the leading
coefficient $c$ of the numerator of $f$; the denominator of $f$ does not
change the roots. `_nf_elem_norm_div` divides the resultant by a known
factor during the multimodular computation, so that only enough primes
for the quotient (with a given bound on its bit size) are needed; this
is useful when the norm is known to be divisible by a large number, as
in ideal arithmetic.

## Trace {#trace}

The trace of $\alpha = \sum_k a_k\theta^k / q$ is
$\sum_k a_k S_k / q$, where $S_k = \sum_j \theta_j^k$ are the power sums
of the roots of $f$. `nf_init` computes $S_0 = d, S_1, \ldots, S_{d-1}$
once, by Newton's identities

$$
c_d S_k + c_{d-1} S_{k-1} + \cdots + c_{d-k+1} S_1 + k\,c_{d-k} = 0,
\qquad 1 \le k < d,
$$

for $f = \sum_k c_k x^k$ [@Coh1996, p. 163]. Since $S_k$ has
denominator dividing $c_d^k$, the values are stored with the common
denominator $c_d^{d-1}$. A trace then costs $d$ multiplications.

## Representation matrices and modular reduction {#matrices}

The matrix of multiplication by $\alpha$ in the basis
$1, \theta, \ldots, \theta^{d-1}$ has as row $j$ the coordinates of
$\theta^j\alpha$; `nf_elem_rep_mat` computes the rows with
`nf_elem_mul_gen`, one shift and one reduction step each. Its
characteristic polynomial is the characteristic polynomial of $\alpha$,
whose constant term and second coefficient give the norm and trace
again.

The functions `nf_elem_mod_fmpz`, `nf_elem_smod_fmpz` and
`nf_elem_coprime_den` replace an element by one congruent to it modulo
an integer $m$ with reduced coefficients, as needed in multimodular and
$p$-adic algorithms over $K$. For `nf_elem_coprime_den`, write the
denominator as $q = q_1q_2$ where $q_1$ contains exactly the primes
dividing $m$; then $q_2$ is invertible modulo $mq_1$, and replacing the
numerator $A$ by $A q_2^{-1} \bmod mq_1$ and the denominator by $q_1$
changes $\alpha$ by an element whose denominator divides $q_2$ and whose
numerator is divisible by $m$.

## Use in FLINT {#use}

- Representation and precomputation: `nf_init`, `nf_elem_init`,
  `nf_elem_canonicalise`.
- Multiplication and reduction: `nf_elem_mul`, `nf_elem_mul_red`,
  `_nf_elem_mul_red`, `nf_elem_reduce`, `_nf_elem_reduce`,
  `nf_elem_mul_gen`, `nf_elem_pow`.
- Inverse and division: `nf_elem_inv`, `nf_elem_div`.
- Norm and trace: `nf_elem_norm`, `nf_elem_norm_div`, `nf_elem_trace`.
- Matrices and modular reduction: `nf_elem_rep_mat`,
  `nf_elem_rep_mat_fmpz_mat_den`, `nf_elem_mod_fmpz`,
  `nf_elem_smod_fmpz`, `nf_elem_coprime_den`.
- The generic ring `gr_ctx_init_nf` uses these functions; `qqbar` and
  `ca` use `nf_elem` for elements of number fields.
