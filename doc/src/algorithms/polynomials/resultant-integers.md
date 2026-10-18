---
title: Resultants and discriminants over Z
summary: Resultants of integer polynomials by the subresultant pseudo-remainder sequence and by the multimodular algorithm with Hadamard's bound, discriminants as resultants with the derivative, and resultants with a known divisor, as implemented in fmpz_poly.
complexity: The multimodular algorithm computes $\operatorname{res}(f, g)$ for degrees $m \ge n$ and $b$-bit coefficients with $O((m + n)(b + \log(m + n)) / \bits)$ primes, each costing a reduction and a modular resultant in $O(M(m) \log m)$ word operations; the subresultant algorithm costs $O(mn)$ operations on integers of the same size as the result.
order: 37
---

The resultant of $f = a \prod_{i=1}^m (x - \alpha_i)$ and $g = b
\prod_{j=1}^n (x - \beta_j)$ in $\Z[x]$,

$$
\operatorname{res}(f, g) = a^n b^m \prod_{i, j} (\alpha_i - \beta_j) =
a^n \prod_{i=1}^m g(\alpha_i),
$$

is the determinant of the $(m + n) \times (m + n)$ Sylvester matrix,
whose first $n$ rows hold shifted copies of the coefficients of $f$
and whose last $m$ rows those of $g$. It is an integer, it vanishes
exactly when $f$ and $g$ have a common root, and it is a polynomial in
the coefficients, which makes it the tool of elimination: for example,
$\operatorname{res}_y(f(y), g(x - y))$ is a polynomial whose roots are
the sums $\alpha_i + \beta_j$. The discriminant

$$
\operatorname{disc}(f) = \frac{(-1)^{m(m-1)/2}}{a} \operatorname{res}(f, f') =
a^{2m - 2} \prod_{i < j} (\alpha_i - \alpha_j)^2
$$

vanishes exactly when $f$ has a repeated root, and for a monic
irreducible $f$ it is, up to the square of an index, the discriminant
of the number field $\mathbb{Q}[x]/(f)$. Over a field the resultant can
be computed with the Euclidean or half-GCD algorithm
([](algorithms/polynomials/half-gcd)); over $\Z$ the questions are
again how to avoid coefficient growth, and how many primes a modular
method needs.

## The size of the resultant {#bound}

::: theorem Hadamard's bound {#thm-hadamard}
For $f, g \in \Z[x]$ of degrees $m$ and $n$,

$$
|\operatorname{res}(f, g)| \le \|f\|_2^{\,n} \, \|g\|_2^{\,m}.
$$
:::

::: proof
Hadamard's inequality [@Had1893] bounds the absolute value of a
determinant by the product of the Euclidean norms of its rows. The
Sylvester matrix has $n$ rows of norm $\|f\|_2$ and $m$ rows of norm
$\|g\|_2$.
:::

Removing contents first gives a better bound, since
$\operatorname{res}(cf, dg) = c^n d^m \operatorname{res}(f, g)$ for
integers $c, d$.

## The subresultant algorithm {#subresultant}

The subresultant pseudo-remainder sequence of
[](algorithms/polynomials/gcd-integers#subresultant) also gives the
resultant: if $f$ and $g$ have no common factor, the sequence ends in
a constant, and the resultant is that constant times known powers of
$g$ and $h$ and a sign, by the multiplicativity of resultants and
$\operatorname{res}(f, g) = (-1)^{mn} \operatorname{res}(g, f)$
(Algorithm 3.3.7 of [@Coh1996]). `_fmpz_poly_resultant_euclidean`
implements it after removing the contents, which it multiplies back at
the end. The intermediate polynomials have coefficients of the size
of subresultants, at most that of the resultant itself, so the
algorithm performs $O(mn)$ operations on integers of $O((m + n)(b +
\log(m + n)))$ bits. It is the fastest method for short polynomials
with small coefficients. The generic `_gr_poly_resultant_subresultant`
is adapted from it and works over any integral domain.

## The multimodular algorithm {#modular}

::: lemma {#lem-mod}
If the prime $p$ divides neither $\operatorname{lc}(f)$ nor
$\operatorname{lc}(g)$, then $\operatorname{res}(f \bmod p, g \bmod p)
= \operatorname{res}(f, g) \bmod p$.
:::

::: proof
The degrees are preserved, so the Sylvester matrix of the reductions
is the reduction of the Sylvester matrix, and reduction modulo $p$
commutes with determinants.
:::

`_fmpz_poly_resultant_modular` [@Col1971] removes the contents,
computes the bound of [](#thm-hadamard) for the primitive parts, and
takes enough primes $p > 2^{63}$ not dividing the product of the
leading coefficients for their product to exceed twice the bound. For
each prime it reduces the polynomials and computes the resultant with
`_nmod_poly_resultant`, which uses the half-GCD algorithm for long
polynomials, and it combines the residues with one multi-CRT
(`fmpz_multi_CRT_ui`) into the symmetric range. Unlike the GCD, the
resultant has no unlucky primes, and the number of primes is fixed in
advance; Hadamard's bound is usually close to the actual size for
random input.

When a factor $d$ of the resultant and a bound for the quotient are
known, `_fmpz_poly_resultant_modular_div` computes
$\operatorname{res}(f, g) / d$ modulo primes that do not divide $d$
(dividing by $d$ modulo each prime) with $\lceil (\mathit{nbits} + 1)/63
\rceil$ primes on 64-bit machines, where `nbits` is the bit bound passed
for the quotient; the extra bit is needed by the reconstruction in the
symmetric range. It also splits the powers of the
contents between the known divisor and the remaining factor. This is
useful when the resultant is known to be divisible by a large number,
for example in discriminant computations for number fields.

## Use in FLINT {#use}

- `fmpz_poly_resultant` uses `fmpz_poly_resultant_euclidean` when the
  shorter polynomial has length $\ell_2 \le 144$ and $\ell_2^3 (b_1 +
  b_2) \le 6 \cdot 10^6$ for coefficient sizes $b_1$, $b_2$ bits, and
  `fmpz_poly_resultant_modular` otherwise.
- `fmpz_poly_discriminant` computes $\operatorname{res}(f, f')$ with
  `_fmpz_poly_resultant`, adjusts the sign and divides by the leading
  coefficient.
- `fmpz_poly_xgcd_modular` starts with the resultant
  ([](algorithms/polynomials/gcd-integers#xgcd)).
- `fmpq_poly_resultant`, `fmpq_poly_discriminant` and
  `fmpq_poly_resultant_modular_div` clear denominators and call the
  `fmpz_poly` functions.

::: functions fmpz_poly_resultant fmpz_poly_resultant_euclidean fmpz_poly_resultant_modular fmpz_poly_resultant_modular_div fmpz_poly_discriminant :::
