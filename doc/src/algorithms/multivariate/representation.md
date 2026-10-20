---
title: Monomial orderings and exponent packing
summary: How FLINT stores sparse multivariate polynomials: terms sorted by a monomial ordering, exponent vectors packed into machine words so that comparing, multiplying and dividing monomials are word operations, with guard bits for overflow detection; and the deflations, repackings and changes of variables that the algorithms apply to this representation.
complexity: Comparison, multiplication and divisibility test of packed monomials cost $O(N)$ word operations for $N$ words per monomial; with $k$ fields of width $b \le 64$ on a 64-bit machine, $N = \lceil k / \lfloor 64/b \rfloor \rceil$.
order: 10
---

Every sparse multivariate module of FLINT (`fmpz_mpoly`, `fmpq_mpoly`,
`nmod_mpoly`, `fmpz_mod_mpoly`, `fq_nmod_mpoly`, `fq_zech_mpoly` and the
generic `gr_mpoly`) uses the same *distributed* representation: a
polynomial is a list of terms $c_i x^{e_i}$ with nonzero coefficients,
sorted so that $x^{e_0} > x^{e_1} > \cdots$ for a fixed monomial ordering,
with the coefficients in one array and the exponent vectors in another.
The exponent vectors are *packed*: several exponents share a machine
word, in a layout chosen so that the monomial operations needed by the
algorithms become a few word operations. This design follows Monagan and
Pearce [@MonPea2007], who showed that packing is a main factor in the
speed of sparse polynomial arithmetic, and Bachmann and Schönemann
[@BacSch1998] for the treatment of orderings. The functions are in the
module [](mpoly), and the notes of Schultz [@Sch-mpoly] give further
details.

## Monomial orderings {#orderings}

A *monomial ordering* is a total order $<$ on the monomials $x^e = x_0^{e_0}
\cdots x_{n-1}^{e_{n-1}}$, $e \in \N^n$, such that $1 \le x^e$ and $x^a <
x^b \Rightarrow x^{a+c} < x^{b+c}$. Such an order is a well-order
(Dickson's lemma), which is what makes the division algorithm and the
heap algorithms terminate [@CLO2015, §2.2]. FLINT supports three, given by
`ordering_t`, with $|e| = e_0 + \cdots + e_{n-1}$:

- `ORD_LEX`: $x^a > x^b$ if the first nonzero entry of $a - b$ is positive;
- `ORD_DEGLEX`: $x^a > x^b$ if $|a| > |b|$, or $|a| = |b|$ and $x^a > x^b$
  in lex order;
- `ORD_DEGREVLEX`: $x^a > x^b$ if $|a| > |b|$, or $|a| = |b|$ and the last
  nonzero entry of $a - b$ is negative.

The variable $x_0$ is the most significant one. The *leading term*
$\operatorname{lt}(A)$, *monomial* $\operatorname{lm}(A)$ and
*coefficient* $\operatorname{lc}(A)$ refer to the first term in this order.
The GCD and factoring code converts its inputs to `ORD_LEX`, in which the
leading coefficient with respect to $x_0$ is a polynomial in the other
variables, as the recursive view $R[x_1, \ldots][x_0]$ requires.

## Packed monomials {#packing}

Fix a *field width* $b$. The exponent vector $e$ is stored as a sequence of
*fields* $f_0, f_1, \ldots, f_{k-1}$, each a $b$-bit unsigned integer, where
$k = n$ for `ORD_LEX` and $k = n + 1$ for the degree orderings:

- for `ORD_LEX` and `ORD_DEGLEX`, field $n - 1 - i$ holds $e_i$, so that
  $e_0$ is in the most significant position among the variables;
- for `ORD_DEGREVLEX`, field $i$ holds $e_i$;
- for the degree orderings, field $n$ holds the total degree $|e|$.

The fields are concatenated, field 0 in the least significant bits, into
the integer
$$
\Phi(e) = \sum_{j < k} f_j \, 2^{b j},
$$
which is stored in $N$ words, least significant word first. Two layouts
are used:

- *single precision*, $8 \le b \le$ `FLINT_BITS`: a word holds $\lfloor
  \mathrm{FLINT\_BITS}/b \rfloor$ fields, fields do not straddle words, and
  unused high bits of a word are zero; $N = \lceil k / \lfloor
  \mathrm{FLINT\_BITS}/b \rfloor \rceil$;
- *multiprecision*, $b$ a multiple of `FLINT_BITS`: each field takes
  $b/\mathrm{FLINT\_BITS}$ words and $N = k b / \mathrm{FLINT\_BITS}$.

In both cases the fields must satisfy $f_j < 2^{b-1}$: the top bit of
each field is a *guard bit*, which is zero in every stored monomial. When
$b$ does not divide `FLINT_BITS` the words are not literally the digits of
$\Phi(e)$, because of the unused bits, but all statements below hold word
by word, and we keep writing $\Phi$ for the packed form.

The width is chosen by `mpoly_fix_bits`, which rounds a required number of
bits up to the largest $b$ with the same number of words: three
variables in lex order with 64-bit words give $b = 21$, two give $b = 32$.
A polynomial starts with $b = 8$ (`MPOLY_MIN_BITS`), and an operation whose
result needs larger fields repacks with a larger width. All monomials of a
polynomial use the same width, recorded in its `bits` field.

### Comparison {#comparison}

Let `cmpmask` be zero for `ORD_LEX` and `ORD_DEGLEX`, and for
`ORD_DEGREVLEX` the pattern with all bits of fields $0, \ldots, n-1$ set and
those of field $n$ clear (`mpoly_get_cmpmask`).

::: lemma Comparison of packed monomials {#lem-cmp}
For each of the three orderings, $x^a > x^b$ if and only if
$\Phi(a) \oplus \mathit{cmpmask} > \Phi(b) \oplus \mathit{cmpmask}$ as
unsigned integers, where $\oplus$ is the bitwise exclusive or.
:::

::: proof
Comparing the integers is comparing the fields lexicographically from the
most significant one. For `ORD_LEX` the fields from the top are $e_0,
e_1, \ldots$, which is the definition. For `ORD_DEGLEX` the total degree
comes first, then $e_0, e_1, \ldots$. For `ORD_DEGREVLEX` the top field is
the total degree, and the remaining fields from the top are $e_{n-1},
\ldots, e_0$, each complemented by the mask. Complementing a $b$-bit field
reverses the order of its values, so for equal degrees $\Phi(a) \oplus
\mathit{cmpmask} > \Phi(b) \oplus \mathit{cmpmask}$ if and only if, at the
last index $i$ with $a_i \ne b_i$, $a_i < b_i$. That is the definition of
reverse lexicographic order.
:::

So comparing monomials costs at most $N$ word comparisons
(`mpoly_monomial_cmp`), and a single one in the frequent case $N = 1$.

The layout relies only on four properties: multiplying monomials adds
fields, comparing them compares the packed words after an exclusive or
with a mask, every exponent is one of the fields, and all fields are
nonnegative.

### Multiplication and overflow {#arithmetic}

Since every field, including the degree field, is additive, $\Phi(a + c) =
\Phi(a) + \Phi(c)$ as long as no field overflows, so monomials are
multiplied by adding words (`mpoly_monomial_add`). The guard bits make
overflow detectable.

::: lemma Overflow {#lem-overflow}
Let $a, c$ have fields $< 2^{b-1}$ and let $M$ be the pattern with the
guard bit of each field set (`mpoly_overflow_mask_sp`). In single
precision the word sum $\Phi(a) + \Phi(c)$ involves no carry between
fields, and it is a valid packed monomial (all fields of $a + c$ are
$< 2^{b-1}$) if and only if $(\Phi(a) + \Phi(c)) \mathbin{\&} M = 0$.
:::

::: proof
Each field of the sum is $< 2^{b-1} + 2^{b-1} = 2^b$, so it fits in its
$b$ bits and nothing is carried into the next field. The field of $a + c$
is then the field of the sum, and it is $< 2^{b-1}$ exactly when its top
bit is zero.
:::

In multiprecision the guard bit is the top bit of the most significant
word of each field (`mpoly_monomial_overflows_mp`). Algorithms usually
avoid testing: before multiplying $B$ by $C$ they compute the maximum of
each field over $B$ and over $C$ (`mpoly_max_fields_fmpz`), add them and
choose $b$ large enough for the sum, so that no product of monomials can
overflow. Overflow tests are needed where results are not known in
advance, as in exact division and in powering.

### Divisibility, halving, minima and maxima {#divisibility}

::: lemma Divisibility test {#lem-divides}
With $a, c$ and $M$ as above, $x^c$ divides $x^a$ if and only if the
word difference $\Phi(a) - \Phi(c)$ satisfies
$(\Phi(a) - \Phi(c)) \mathbin{\&} M = 0$, and then $\Phi(a) - \Phi(c) =
\Phi(a - c)$.
:::

::: proof
If $a_j \ge c_j$ for all fields, the subtraction proceeds field by field
without borrows and gives $\Phi(a - c)$, whose guard bits are zero.
Otherwise let $j$ be the lowest field with $a_j < c_j$. All lower fields
are subtracted without borrow, and field $j$ becomes $a_j - c_j + 2^b >
2^b - 2^{b-1} = 2^{b-1}$, borrowing from the field above, so its guard bit
is set.
:::

`mpoly_monomial_divides` performs this test (word by word, so a borrow
never crosses a word boundary unnoticed), and the exact division code uses
it to detect a leading monomial of the dividend that is not divisible by
that of the divisor. The same trick gives, for square roots
(`mpoly_monomial_halves`): $x^a$ is a square if and only if the lowest bit
of each word is zero and $(\Phi(a) \gg 1) \mathbin{\&} M = 0$, since
shifting a word right by one moves the lowest bit of each field into the
guard bit of the field below. Field-wise minima and maxima are computed
without unpacking: in $s = M + \Phi(a) - \Phi(c)$ each field equals $2^{b-1}
+ a_j - c_j \in [1, 2^b)$, so there are no borrows and the guard bit of
field $j$ is set exactly when $a_j \ge c_j$; spreading it over the field
gives a mask selecting $a_j$ or $c_j$ (`mpoly_monomial_min1`,
`mpoly_monomial_max`).

## Other representations {#other}

The distributed representation suits arithmetic, but some algorithms want
a recursive or dense view, to which FLINT converts temporarily:

- *Recursive sparse*: `nmod_mpolyu_t` and its relatives are sparse
  polynomials in one or two main variables whose coefficients are sparse
  polynomials in the others, and `nmod_mpolyn_t` is a sparse polynomial
  whose coefficients are dense univariate polynomials in the last variable.
  Brown's GCD algorithm works with the latter, Zippel's with the former
  ([](algorithms/multivariate/gcd)).
- *Univariate with multivariate coefficients*: `fmpz_mpoly_univar_t` and
  its relatives, public types used for resultants, discriminants and
  contents.
- *Dense*: `fmpz_mpolyd_t` stores all coefficients in a box of given
  degree bounds and is used by `fmpz_mpoly_mul_dense`; the array methods of
  multiplication and division use temporary dense arrays of coefficients
  indexed by the exponents ([](algorithms/multivariate/multiplication)).
- *Dense bivariate*: `fmpz_bpoly_t`, `n_bpoly_t` and their relatives,
  arrays of univariate polynomials, used in bivariate factorisation and
  resultants.

## Deflation, compression and repacking {#deflation}

Many algorithms first simplify the exponents. Write $s_i$ for the minimal
exponent of $x_i$ in $A$ and $t_i$ for the GCD of the differences $e_i -
s_i$ over the terms. Then $A = x^s \, \iota_t(A')$ where $\iota_t$ is the
substitution $x_i \mapsto x_i^{t_i}$ and $A'$ is the *deflated*
polynomial, with exponents $(e_i - s_i)/t_i$ (`mpoly_gcd_info_stride`,
`fmpz_mpoly_deflation`, `fmpz_mpoly_deflate`, `fmpz_mpoly_inflate`).

::: lemma GCDs commute with substitution of powers {#lem-inflate}
Let $D$ be $\Z$ or a field, let $t_i \ge 1$, and let $\iota$ be the ring
homomorphism of $S = D[x_0, \ldots, x_{n-1}]$ to itself given by $x_i
\mapsto x_i^{t_i}$. Then $\gcd(\iota(A), \iota(B)) = \iota(\gcd(A, B))$ up
to units.
:::

::: proof
Let $R$ be the image of $\iota$, the ring $D[y_0, \ldots, y_{n-1}]$ with
$y_i = x_i^{t_i}$; $\iota$ is an isomorphism onto it. The ring $S$ is a
free $R$-module with basis the monomials $x^j$, $0 \le j_i < t_i$, hence
flat over $R$, and for a flat extension $(aS) \cap (bS) = (aR \cap bR) S$.
In the unique factorisation domains $R$ and $S$ the intersection of two
principal ideals is generated by the least common multiple, so
$\operatorname{lcm}_S(a, b)$ and $\operatorname{lcm}_R(a, b)$ agree up to
units of $S$, and so do the greatest common divisors $ab /
\operatorname{lcm}$.
:::

Thus the GCD of $A = x^{s} \iota(A')$ and $B = x^{s'} \iota(B')$ with
common strides $t$ is $x^{\min(s, s')} \iota(\gcd(A', B'))$, since
$\iota(A')$ and $\iota(B')$ are not divisible by any variable. The GCD
code uses this to reduce degrees. The factoring code goes further and
applies an invertible integer affine map to the exponents
(*compression*, `mpoly_compression_set`), described in
[](algorithms/multivariate/factorisation#compression).

When exponents need more than one word per field, most algorithms first
try to *repack* them into one word (`fmpz_mpoly_repack_bits`), which
succeeds whenever all fields are below $2^{\mathrm{FLINT\_BITS}-1}$; the GCD
code also deflates. Algorithms that only exist for single-word fields
(array multiplication, dense methods, GCD and factoring internals) decline
or fail otherwise.

## Use in FLINT {#use}

- Contexts and layout: `mpoly_ctx_init` (field order, `nfields`),
  `mpoly_words_per_exp`, `mpoly_fix_bits`, `mpoly_get_cmpmask`,
  `mpoly_overflow_mask_sp`; packing and unpacking with `mpoly_pack_vec_ui`,
  `mpoly_set_monomial_ui`, `mpoly_get_monomial_ui` and their `fmpz`
  versions; repacking with `mpoly_repack_monomials`.
- Monomial operations ([](#lem-cmp), [](#lem-overflow),
  [](#lem-divides)): `mpoly_monomial_cmp`, `mpoly_monomial_add`,
  `mpoly_monomial_overflows`, `mpoly_monomial_divides`,
  `mpoly_monomial_halves`, `mpoly_monomial_max`, `mpoly_monomial_min1` and
  their `_mp` versions.
- Deflation ([](#lem-inflate)): `mpoly_gcd_info_stride`,
  `mpoly_to_mpolyl_perm_deflate`, `mpoly_from_mpolyl_perm_inflate`,
  `fmpz_mpoly_deflation`, `fmpz_mpoly_deflate`, `fmpz_mpoly_inflate`.

::: functions mpoly_monomial_cmp mpoly_monomial_divides mpoly_fix_bits mpoly_get_cmpmask fmpz_mpoly_repack_bits fmpz_mpoly_deflate :::
