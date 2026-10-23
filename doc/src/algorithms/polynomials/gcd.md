---
title: Polynomial GCDs and resultants
kind: topic
summary: The Euclidean algorithm, the half-gcd algorithm and the subresultant pseudo-remainder sequence as implemented generically in gr_poly, how the generic code chooses between them from the properties of the coefficient ring, how resultants are derived from the same sequences, and the specialised algorithms of nmod_poly, fmpz_mod_poly, fq_poly and fmpz_poly.
order: 35
---

A greatest common divisor of two univariate polynomials $A$ and $B$ can
be computed over any field by the Euclidean algorithm, and over any
unique factorisation domain by a pseudo-remainder sequence, which
avoids divisions. These algorithms are implemented once in
[`gr_poly`](gr_poly#gcd) for any coefficient ring. Which one is fast
depends on the ring: the Euclidean algorithm and its asymptotically
fast form, the half-gcd algorithm, are the methods of choice over finite
fields; over $\mathbb{Q}$ and $\Z$ their intermediate coefficients grow, and the
best algorithms reduce the problem modulo primes instead.

Throughout, $A$ and $B$ are nonzero with $m = \deg A \ge n = \deg B$,
and $\operatorname{lc}$ denotes the leading coefficient.

## The Euclidean algorithm {#euclidean}

Over a field $K$, the Euclidean remainder sequence is $R_0 = A$, $R_1 =
B$ and $R_{i+1} = R_{i-1} \bmod R_i$, until $R_{k+1} = 0$; then $R_k$ is a
GCD of $A$ and $B$, and $\gcd(A, B)$ is its monic associate. Each step
costs $O(\deg R_i \cdot (\deg R_{i-1} - \deg R_i + 1))$ operations, and
the whole sequence $O(mn)$ operations [@GG2003, Ch. 3]. Keeping track of
cofactors $S_i$, $T_i$ with $S_i A + T_i B = R_i$ (starting from $S_0 =
1$, $T_0 = 0$, $S_1 = 0$, $T_1 = 1$ and applying the quotients) gives the
extended GCD at the same asymptotic cost.

The resultant follows from the remainder sequence. If $R = A \bmod B$
has degree $k$, then

$$
\operatorname{res}(A, B) = (-1)^{mn} \operatorname{lc}(B)^{m - k}
\operatorname{res}(B, R),
$$

while $\operatorname{res}(A, B) = 0$ if $R = 0$ and $n > 0$, and
$\operatorname{res}(A, B) = \operatorname{lc}(B)^m$ if $n = 0$. The
Euclidean resultant algorithm accumulates these factors along the
sequence ([](gr_poly#resultant)).

In `gr_poly` these are `_gr_poly_gcd_euclidean`,
`_gr_poly_xgcd_euclidean` and `_gr_poly_resultant_euclidean`. They need
to invert the leading coefficient of each remainder and to decide when
a remainder vanishes. Over a field with decidable equality this always
works. Over $\Z/n\Z$ with composite $n$ pretending to be a field
(`gr_ctx_set_is_pretend_field`), a non-invertible leading coefficient
makes the algorithm fail with an error status, and the factor of $n$
it has found can be read with `gr_ctx_recover_zero_divisor`, which is
the basis of dynamic evaluation [@DDDD1985]. Over inexact rings such as
real balls the algorithm succeeds when the zero tests are decided,
typically for exactly represented inputs, and returns `GR_UNABLE` when
a remainder cannot be shown to be zero or nonzero: the GCD of
polynomials known only approximately is not a well-posed problem.

## The half-gcd algorithm {#half-gcd}

The quotients of the first half of the Euclidean remainder sequence
depend only on the top halves of $A$ and $B$. The half-gcd algorithm
computes the $2 \times 2$ polynomial matrix that maps $(A, B)$ to a
pair of consecutive remainders of about half the degree, recursively,
and repeats until the remaining polynomials are short, giving the GCD,
the extended GCD and the resultant with $O(M(n) \log n)$ operations
[@Moe1973; @BGY1980; @ThullYap1990]. The generic implementation
`_gr_poly_hgcd` uses the Euclidean algorithm below an inner cutoff, and
`_gr_poly_gcd_hgcd`, `_gr_poly_xgcd_hgcd` and `_gr_poly_resultant_hgcd`
call it repeatedly above an outer cutoff. Its matrix products are done
by the ring method `_gr_poly_hgcd_mat_mul`, with Strassen's seven
products and, when the ring has an FFT representation, products of
transformed entries ([](gr_poly#gcd)). Over infinite fields such as $\mathbb{Q}$
the half-gcd algorithm is not used by default, because its
intermediate coefficients are no smaller than those of the Euclidean
algorithm while its overhead is larger, and over approximate fields
because of numerical instability.

::: summary algorithms/polynomials/half-gcd :::

## Pseudo-remainder sequences {#prs}

Over an integral domain $D$ that is not a field, $A \bmod B$ is not
defined, but the *pseudo-remainder*

$$
\operatorname{prem}(A, B) = \operatorname{lc}(B)^{m - n + 1} A \bmod B
$$

has coefficients in $D$. Repeated pseudo-division leads to a sequence
whose remainders are associate (over the fraction field) to the
Euclidean remainders, but whose coefficients grow exponentially in the
number of steps. The subresultant pseudo-remainder sequence of Collins
and of Brown and Traub [@Col1967; @Bro1971; @BT1971] divides each
pseudo-remainder by a factor known in advance to divide it exactly, so
that its elements are, up to sign, the subresultants of $A$ and $B$:
determinants of submatrices of the Sylvester matrix, whose size is
bounded by Hadamard's inequality and grows only linearly in the step
number.

The generic `_gr_poly_gcd_subresultant` follows Cohen's formulation
[@Coh1996, Algorithm 3.3.1] over a unique factorisation domain:

1. Compute the contents $a = \operatorname{cont}(A)$ and $b =
   \operatorname{cont}(B)$ (GCDs of the coefficients, using the
   ring's GCD) and $d = \gcd(a, b)$; replace $A$ and $B$ by their
   primitive parts, and set $g = h = 1$.
2. Let $\delta = \deg A - \deg B$ and $R = \operatorname{prem}(A, B)$.
   If $R = 0$, stop with $B$; if $\deg R = 0$, the GCD is $d$.
3. Set $A \leftarrow B$, $B \leftarrow R / (g h^\delta)$ (an exact
   division), $g \leftarrow \operatorname{lc}(A)$ and $h \leftarrow
   h^{1 - \delta} g^\delta$ (also exact), and go to step 2.
4. The GCD is $d$ times the primitive part of the final $B$.

The extended version `_gr_poly_xgcd_subresultant` and the resultant
version `_gr_poly_resultant_subresultant` follow the same sequence. The
resultant version needs only exact divisions and works over any
integral domain with exact division.

As a last resort, `_gr_poly_resultant_sylvester` evaluates the
determinant of the Sylvester matrix with `gr_mat_det`, which uses a
division-free algorithm by default and therefore works over any
commutative ring, including rings with zero divisors.

## Generic algorithm selection {#generic}

`gr_poly_gcd`, `gr_poly_xgcd` and `gr_poly_resultant` call ring
methods, whose defaults choose from the properties of the ring:

| Ring | GCD and extended GCD | Resultant |
|---|---|---|
| finite field (or finite ring pretending to be a field) | Euclidean below length 200, half-gcd above (inner cutoff 100) | small formulas, then Euclidean or half-gcd with the same cutoffs |
| other field ($\mathbb{Q}$, number fields, `qqbar`, `ca`, balls, floating point) | Euclidean | small formulas, then Euclidean |
| unique factorisation domain ($\Z$, $\Z[y]$, $\Z[i]$) | subresultant | small formulas, then subresultant |
| other rings ($\Z/6\Z$) | `GR_UNABLE` | Sylvester determinant |

For the resultant, each failing algorithm falls back to the next one in
the list subresultant, Sylvester determinant. The "small formulas"
(`_gr_poly_resultant_small`) are division-free closed forms for a
divisor of length at most 3 or a dividend of length at most 2. The
non-underscore functions canonicalise the GCD: monic over a field,
with a canonical unit as leading coefficient over other rings (positive
over $\Z$); see [](gr_poly/gcd) for the conventions and an example.

## Specialised algorithms {#specialisations}

Several rings override the methods, mostly to tune the cutoffs between
the Euclidean and half-gcd algorithms, which depend on the cost of the
ring's multiplication:

| Ring or module | GCD algorithm |
|---|---|
| `nmod_poly` (word-size prime $p$) | Euclidean below a cutoff tabulated by the bit size of $p$ (470 to 2199), with a Montgomery-form inner loop for some moduli; half-gcd above. The resultant uses the same cutoff, the extended GCD its own table. |
| `fmpz_mod_poly` (multi-word prime) | generic Euclidean below length 256, half-gcd above with inner cutoff 128 |
| `mpn_mod` | Euclidean below length 240, half-gcd above with inner cutoff 80 |
| `fq_poly`, `fq_nmod_poly`, `fq_zech_poly` | generic Euclidean below 80 to 120 (depending on the variant and on whether the characteristic has at most 8 bits), half-gcd above with inner cutoff 25 to 35 |
| `fmpz_poly` | subresultant below length 6; otherwise the heuristic GCD, which evaluates at a large integer and reconstructs from one integer GCD, when the coefficients are small, and otherwise (or if the heuristic fails) the multimodular algorithm ([](algorithms/polynomials/gcd-integers)) |
| `fmpq_poly` | GCD of the primitive parts with `fmpz_poly_gcd` |
| multivariate polynomials | Brown's, Zippel's and Hensel-lifting algorithms ([](algorithms/multivariate/gcd)) |

`gr_poly_gcd` over `fmpz` calls `_fmpz_poly_gcd` through the ring
method, but over `fmpq` it runs the generic Euclidean algorithm, which
is much slower than `fmpq_poly_gcd` for large inputs. Resultants are
better served: polynomial rings over $\Z/p\Z$, $\Z$ and $\mathbb{Q}$ override
the resultant method with bivariate algorithms by evaluation and
interpolation and by Chinese remaindering
([](algorithms/multivariate/modp-bivariate-resultant),
[](algorithms/polynomials/resultant-integers)).

::: summary algorithms/polynomials/gcd-integers :::

::: summary algorithms/polynomials/resultant-integers :::
