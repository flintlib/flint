---
title: Choosing types
summary: How to choose among FLINT's types for integers, residues, finite fields, polynomials, matrices, p-adic and algebraic numbers, real and complex numbers and power series, with their limits, their performance and the generic alternative.
order: 1
---

FLINT usually has more than one type for a given mathematical object.
The integers modulo $n$ alone can be `ulong` residues with an `nmod_t`,
`fmpz` residues with an `fmpz_mod_ctx_t`, fixed-length limb arrays in an
`mpn_mod` context, or elements of any of these through the generic
interface [`gr`](gr). The types differ in the range of parameters they
accept, in speed (often by large factors for small inputs), in the
functionality that exists for them, and in how they are called. The
sections below give the general principles and tables of the types by
domain, with their limits.

## General principles {#principles}

*Use the narrowest representation that fits.* A type whose elements
have a bounded size, such as a word-size residue, needs no memory
allocation, no branches on the size of the operands and no pointer
chasing, and can be packed into vectors and matrices on which
vectorised kernels operate. If the modulus is below $2^{64}$, `nmod`
types are several times faster than `fmpz_mod` types for the same
computation; for small non-prime finite fields, Zech logarithms
are much faster than polynomial representations. The price is a
limit on the parameters, which the tables below give for each type.

*Use the type of the module that does what you need.* The classic
modules (`fmpz_poly`, `nmod_mat`, `fmpz_mpoly`, ...) contain algorithms
that exist for one ring only: factorisation of integer polynomials, LLL
reduction, Hermite normal forms, multivariate GCDs. Before choosing a
representation for speed, check that the module of that type has the
operations your program needs.

*Decide between exact and approximate arithmetic.* Exact types (`fmpz`,
`fmpq`, residues, finite fields, number fields, `qqbar`, `ca`) give
exact results whose size may grow. Ball types (`arb`, `acb` and the
`dfloat` and `decimal` balls) give rigorous error bounds at a fixed
working precision. Floating-point types (`arf`, `nfloat`, the `dfloat`
numbers, `double`) give approximations without automatic error tracking. See
[](guide/types/ball-arithmetic) and [](#exact).

*Use the generic interface for rings that have no module of their own
and for code that must work over several rings.* Polynomials over
$\Z[i]$, matrices over a number field, power series over $\Z/n\Z$ for a
large $n$ or polynomials whose coefficients are power series have no
dedicated type, but are available as `gr_poly`, `gr_mat` and so on over
the right context. Where a classic module exists, its functions are at
least as fast as the generic ones (the generic contexts call them for
the critical operations) and usually more complete. The trade-offs are
described in [](guide/generic).

*Mind the context.* Types whose parameters are not stored in the
object (`fmpz_mod`, `fq`, `padic`, the multivariate types, `gr`) take a
context argument, which must outlive the objects. For historical reasons,
the `nmod_poly` and `nmod_mat` types store their modulus in each object,
which is simpler to call but makes objects larger.

## Integers {#integers}

| Type | Represents | Limits | Use it for |
|---|---|---|---|
| `ulong`, `slong` with [`ulong_extras`](ulong_extras) (`n_` functions) and [`long_extras`](long_extras) (`z_`) | machine words | $0 \le n < 2^{\bits}$, resp. $\lvert n \rvert < 2^{\bits - 1}$; overflow wraps silently | inner loops; number theory on words (`n_is_prime`, `n_factor`, `n_gcd`, `n_mulmod_preinv`); exponents and moduli |
| `fmpz_t` in [`fmpz`](fmpz) | an integer of any size | memory; values with $\lvert n \rvert < 2^{\bits - 2}$ are stored in the word itself, larger ones as GMP integers ([](fmpz/representation)) | the default integer type; fast for small values, asymptotically fast multiplication, division and GCD for large ones |
| `nn_ptr` limb arrays with [`mpn_extras`](mpn_extras) and GMP's `mpn` layer | nonnegative integers as limb arrays | the caller manages lengths and memory | fixed-size multiprecision kernels and implementing other types; the fastest code for operands of a few limbs |
| `radix_integer_t` in [`radix`](radix) | integers in a radix such as $10^{19}$ | | decimal or base-$p$ digit manipulation and conversion in linear time |
| `fmpzi_t` in [`fmpzi`](fmpzi) | Gaussian integers $\Z[i]$ | memory | Euclidean division, GCD and primality in $\Z[i]$ |
| `fmpq_t` in [`fmpq`](fmpq) | a rational number in lowest terms | memory | exact rational arithmetic; vectors and matrices in [`fmpq_vec`](fmpq_vec) and [`fmpq_mat`](fmpq_mat) |

## Integers modulo n {#mod-n}

| Type | Modulus | Elements | Use it for |
|---|---|---|---|
| `ulong` with `nmod_t`, in [`nmod`](nmod) | $1 \le n < 2^{\bits}$ | reduced residues in one word | the default whenever $n$ fits in a word; vectors, polynomials and matrices in [`nmod_vec`](nmod_vec), [`nmod_poly`](nmod_poly), [`nmod_mat`](nmod_mat), [`nmod_mpoly`](nmod_mpoly) |
| `nmod8`, `nmod32`, `nmod_redc`, `nmod_redc_fast` contexts of `gr` ([](nmod)) | $n \le 255$, $n < 2^{32}$, odd $n \ge 3$, odd $3 \le n < 2^{\bits - 2}$ | `uint8_t`, `uint32_t`, Montgomery forms | saving memory in large generic vectors and matrices; mostly less optimised than `nmod` |
| `mpn_mod` contexts of `gr`, in [`mpn_mod`](mpn_mod) | $2^{\bits} \le n < 2^{16 \bits}$ (2 to 16 limbs) | fixed-length limb arrays, no heap memory | generic polynomials and matrices modulo moduli of 65 to 1024 bits; faster than `fmpz_mod` there ([](mpn_mod/usage)) |
| `fmpz_t` with `fmpz_mod_ctx_t`, in [`fmpz_mod`](fmpz_mod) | any $n \ge 1$ | `fmpz` residues in $[0, n)$ | moduli of any size, with the full classic set of modules [`fmpz_mod_poly`](fmpz_mod_poly), [`fmpz_mod_mat`](fmpz_mod_mat), [`fmpz_mod_mpoly`](fmpz_mod_mpoly) |

Operations that divide (polynomial division and GCD, matrix inversion
and solving) generally require a prime modulus, or at least invertible
pivots and leading coefficients. With the classic types the requirement
is a documented precondition, and violating it raises an exception or
gives undefined results; the generic `gr` interface instead returns
`GR_DOMAIN` or `GR_UNABLE`.

## Finite fields {#finite-fields}

| Type | Field $\mathbb{F}_q$, $q = p^d$ | Elements | Use it for |
|---|---|---|---|
| [`fq_zech`](fq_zech) | $q < 2^{\bits}$; context tables of $2q + p$ words | discrete logarithms in one word | small fields: multiplication and division are additions of logarithms ([](algorithms/finite-fields/zech-logarithms)) |
| [`fq_nmod`](fq_nmod) | $p < 2^{\bits}$, any $d$ | `nmod_poly` of length $< d$ | extension fields of word-size characteristic |
| [`fq`](fq) | any $p$, any $d$ | `fmpz_poly` of length $< d$ | large characteristic |
| [`fq_default`](fq_default) | any | a union of the above, `ulong` or `fmpz` | not knowing in advance: the representation is chosen when the context is created |
| `nmod`, `fmpz_mod` and `mpn_mod` with a prime modulus | prime fields | residues | $d = 1$: faster than an `fq` type of degree 1 |

`fq_default_ctx_init` chooses `fq_zech` if $d > 1$ and $d$ times the bit
length of $p$ is at most 16 (so $q < 2^{16}$), `fq_nmod` if $d > 1$ and
$p < 2^{\bits}$, `nmod` or `fmpz_mod` if $d = 1$, and `fq` otherwise
(`fq_default_ctx_init_type` forces a choice); see
[](fq_default/getting-started). Each representation has modules for
vectors, matrices, polynomials, polynomial factorisation and embeddings
(`fq_nmod_poly`, `fq_zech_mat`, `fq_default_poly_factor`, ...), and
`fq_nmod` and `fq_zech` also multivariate polynomials
([`fq_nmod_mpoly`](fq_nmod_mpoly), [`fq_zech_mpoly`](fq_zech_mpoly)).

## Univariate polynomials and power series {#polynomials}

| Type | Coefficients | Use it for |
|---|---|---|
| [`nmod_poly`](nmod_poly) | $\Z/n\Z$, $n < 2^{\bits}$ | the fastest dense polynomials: multiplication by Kronecker substitution or the small-prime FFT ([](nmod_poly/mul)); factorisation in [`nmod_poly_factor`](nmod_poly_factor) |
| [`fmpz_poly`](fmpz_poly) | $\Z$ | integer polynomials; factorisation ([`fmpz_poly_factor`](fmpz_poly_factor)), real roots, resultants, cyclotomic and other special polynomials |
| [`fmpq_poly`](fmpq_poly) | $\mathbb{Q}$ | rational polynomials, stored as an `fmpz_poly` numerator with one common denominator, so that the algorithms are those of `fmpz_poly` |
| [`fmpz_mod_poly`](fmpz_mod_poly) | $\Z/n\Z$, any $n$ | polynomials over residue rings with a multiprecision modulus; factorisation in [`fmpz_mod_poly_factor`](fmpz_mod_poly_factor) |
| [`fq_poly`](fq_poly), [`fq_nmod_poly`](fq_nmod_poly), [`fq_zech_poly`](fq_zech_poly), [`fq_default_poly`](fq_default_poly) | $\mathbb{F}_q$ | as for the field elements above |
| [`arb_poly`](arb_poly), [`acb_poly`](acb_poly) | real or complex balls | rigorous polynomial arithmetic and roots, and truncated power series: the `_series` functions are FLINT's automatic differentiation of special functions |
| [`padic_poly`](padic_poly), [`ca_poly`](ca_poly) | $\mathbb{Q}_p$, exact complex numbers | |
| [`gr_poly`](gr_poly) | any `gr` ring | rings without a dedicated module ($\Z[i]$, number fields, `mpn_mod`, `nfloat`, polynomials over polynomials), and generic code |
| [`gr_series`](gr_series) | any `gr` ring | power series as ring elements, either with an error term $O(x^n)$ or exactly in $R[x]/(x^n)$ |
| [`fmpz_poly_q`](fmpz_poly_q) | $\mathbb{Q}(x)$ | univariate rational functions |

Lengths and degrees are `slong`; the practical limit is memory. Most
classic polynomial modules have `_series` functions (inverse, square
root, exponential, logarithm, composition and reversion of power series
truncated to a length).

## Multivariate polynomials {#multivariate}

| Type | Coefficients | Use it for |
|---|---|---|
| [`fmpz_mpoly`](fmpz_mpoly), [`fmpq_mpoly`](fmpq_mpoly) | $\Z$, $\mathbb{Q}$ | sparse polynomial arithmetic, exact division, GCD and factorisation ([`fmpz_mpoly_factor`](fmpz_mpoly_factor)) |
| [`nmod_mpoly`](nmod_mpoly), [`fmpz_mod_mpoly`](fmpz_mod_mpoly) | $\Z/n\Z$ | the same over residue rings (GCD and factorisation need a prime modulus) |
| [`fq_nmod_mpoly`](fq_nmod_mpoly), [`fq_zech_mpoly`](fq_zech_mpoly) | $\mathbb{F}_q$, word-size $p$ | the same over finite fields |
| [`fmpz_mpoly_q`](fmpz_mpoly_q), [`fmpz_mod_mpoly_q`](fmpz_mod_mpoly_q) | | multivariate rational functions |
| [`gr_mpoly`](gr_mpoly) | any `gr` ring | arithmetic and division over other rings; no GCD or factorisation |

All these types store terms sparsely, sorted in a monomial ordering of
the context, with exponent vectors packed into fields of `bits` bits
([](mpoly)). Exponents are not bounded: up to `FLINT_BITS` bits, several
fields share a word (at least `MPOLY_MIN_BITS` = 8 bits per field), and
larger exponents take several words per field; the packing is enlarged
automatically when a result needs it. `fmpq_mpoly` stores a rational
content times a primitive `fmpz_mpoly`. For univariate problems,
it is generally better to use the univariate types which store polynomials
densely.

## Matrices {#matrices}

| Type | Entries | Use it for |
|---|---|---|
| [`nmod_mat`](nmod_mat) | $\Z/n\Z$, $n < 2^{\bits}$ | fast linear algebra modulo word-size $n$, with blocked and vectorised multiplication ([](nmod_mat/mul), [](nmod_mat/linear-algebra)); the base of the multimodular algorithms over $\Z$ |
| [`fmpz_mat`](fmpz_mat) | $\Z$ | exact integer linear algebra: multimodular and $p$-adic solving, determinants, Hermite and Smith normal forms; LLL in [`fmpz_lll`](fmpz_lll) |
| [`fmpq_mat`](fmpq_mat) | $\mathbb{Q}$ | rational matrices; nontrivial operations clear denominators and use `fmpz_mat` |
| [`fmpz_mod_mat`](fmpz_mod_mat), [`fq_mat`](fq_mat), [`fq_nmod_mat`](fq_nmod_mat), [`fq_zech_mat`](fq_zech_mat), [`fq_default_mat`](fq_default_mat) | residues, finite fields | linear algebra over these rings |
| [`arb_mat`](arb_mat), [`acb_mat`](acb_mat) | real or complex balls | rigorous numerical linear algebra, eigenvalues, matrix exponentials |
| [`ca_mat`](ca_mat) | exact real and complex numbers | exact linear algebra over $\overline{\mathbb{Q}}$ and transcendental extensions |
| [`nmod_poly_mat`](nmod_poly_mat), [`fmpz_poly_mat`](fmpz_poly_mat), [`padic_mat`](padic_mat), [`bool_mat`](bool_mat) | polynomials, $\mathbb{Q}_p$, Booleans | |
| [`gr_mat`](gr_mat) | any `gr` ring | other rings (number fields, `mpn_mod`, `nfloat`, polynomial rings), approximate floating-point linear algebra, and generic code |
| [`d_mat`](d_mat) | `double` | a small support module with basic operations only |

For approximate (non-rigorous) linear algebra at or above double
precision, `gr_mat` over an `nfloat` or `dfloat` context is the fastest
choice; for certified results, use `arb_mat` or `acb_mat`.

## p-adic numbers and number fields {#padic}

| Type | Represents | Limits | Use it for |
|---|---|---|---|
| [`padic`](padic) | $\mathbb{Q}_p$: $p^v u$ modulo $p^N$ | any prime $p$; precision $N$ per variable | the classic $p$-adic type, with polynomials and matrices ([`padic_poly`](padic_poly), [`padic_mat`](padic_mat)) |
| [`qadic`](qadic) | unramified extensions $\mathbb{Q}_q$ of degree $d$ | any $p$, any $d$ | Frobenius, norms and traces, Teichmüller lifts |
| [`padic_radix`](padic_radix) | $\mathbb{Q}_p$ with propagated error terms | word-size $p$ | a `gr` type with digits in radix $p^e$, faster than `padic` at high precision ([](padic_radix/getting-started)) |
| [`nf_elem`](nf_elem) with [`nf`](nf) | number fields $\mathbb{Q}[x]/(f)$ | $f$ irreducible of degree $\ge 1$ (not checked) | arithmetic, norms and traces; special fixed-size representations for degree 1 and 2 |

## Real and complex numbers {#real-complex}

| Type | Represents | Precision and range | Use it for |
|---|---|---|---|
| [`arb`](arb), [`acb`](acb) | balls $[m \pm r]$ with an `arf` midpoint and a `mag` radius | any precision in bits, unbounded exponents | rigorous numerics and special functions (the default; see [](guide/types/ball-arithmetic)) |
| [`arf`](arf), [`mag`](mag), [`acf`](acf) | binary floating-point numbers; radii with a 30-bit mantissa | `arf`: any precision, unbounded exponents, mantissas of up to two limbs stored inline | midpoints and radii of balls; approximate computation with correctly rounded basic operations |
| [`nfloat`](nfloat) | floating-point numbers of a fixed number of words, real or complex | 1 to 66 words (64 to 4224 bits on 64-bit machines); exponents below about $2^{61}$ in magnitude | fast approximate linear algebra and root finding through `gr_mat` and `gr_poly` ([](nfloat/getting-started)) |
| [`dfloat`](dfloat) | sums of $N = 1, \ldots, 4$ doubles, and balls of them, real or complex | about $53N$ bits; the exponent range of `double` | fast approximate or rigorous arithmetic up to about 200 bits; requires correctly rounded `fma` (`dfloat_is_supported`) |
| [`decimal`](decimal) | decimal floating-point numbers and balls | any number of decimal digits, seven rounding modes | results specified in decimal digits; cheap conversion to and from decimal strings |
| [`mp_real`](mp_real) | balls and fixed-point numbers in the radix $2^{\bits}$ | limb granularity | a low-level backend for ball arithmetic at high precision |
| `double` with [`double_extras`](double_extras), [`d_vec`](d_vec), [`d_mat`](d_mat), and `di_t` in [`double_interval`](double_interval) | IEEE doubles; intervals with double endpoints | 53 bits | cheap estimates and cheap rigorous bounds |
| [`arb_fpwrap`](arb_fpwrap) | special functions of `double` arguments | results accurate to double precision | calling FLINT's special functions from floating-point code |

The ball types are the only real and complex types whose results are
guaranteed in all cases (the `dfloat` and `decimal` balls are rigorous
as well, but cover less functionality). The floating-point types are
for approximate computations, typically followed by a rigorous
verification, and are not rings in the sense of `gr`
([](algorithms/generic/inexact)).

## Exact real and complex numbers {#exact}

Ball arithmetic can prove that two numbers differ, by computing them
accurately enough that their balls are disjoint, but never that they are
equal: a ball that contains zero is compatible with a tiny nonzero
number. Computations that need equality (is a pivot zero, does an
expression reproduce a known constant) can use FLINT's exact real and
complex numbers, which represent numbers symbolically, so that
arithmetic is exact and numerical values are available to any
precision. This part of FLINT was developed as the library Calcium
[@Joh2021calcium].

| Type | Represents | Use it for |
|---|---|---|
| [`qqbar`](qqbar) | algebraic numbers, stored canonically by a minimal polynomial and an isolating ball | individual algebraic numbers of moderate degree; every predicate is decided ([](qqbar/using)) |
| [`ca`](ca) | elements of fields $\mathbb{Q}(a_1, \ldots, a_n)$ generated by algebraic and transcendental numbers ($\mathbb{Q}(\sqrt 2, \pi)$, ...), constructed automatically; also infinities and *Undefined* | exact computation with algebraic and elementary transcendental numbers, with vectors, polynomials and matrices in [`ca_vec`](ca_vec), [`ca_poly`](ca_poly), [`ca_mat`](ca_mat) ([](ca/using)) |
| [`fexpr`](fexpr) | symbolic expressions such as `Add(Sqrt(2), Pi)` | constructing numbers from formulas, printing them as text or LaTeX; not a number type ([](fexpr/using)) |

The predicates of `ca` return `T_TRUE`, `T_FALSE` or `T_UNKNOWN`
([](calcium#truth)): equality of algebraic numbers is decided within
configurable limits, and equality involving transcendental numbers in
many but not all cases, since it is not known to be decidable in
general. A program must handle `T_UNKNOWN`
([](ca/check-is-zero#unknown), [](algorithms/exact/zero-testing#in-practice)).
`qqbar` decides every predicate, but its minimal polynomials can grow
quickly; computations whose degrees explode, such as sums of many square
roots, can be cheaper with `ca`, which keeps them in a field with
several generators. For numbers in $\mathbb{Q}$ or in one fixed number
field, [`fmpq`](fmpq) and [`nf_elem`](nf_elem) are much faster.

Exact numbers are far too slow to replace floating-point or ball
arithmetic in numerical computations. A mostly numerical computation
should use balls, and fall back on exact numbers only where an exact
comparison is needed; exact numbers are also useful for testing
numerical code. The contexts `gr_ctx_init_real_qqbar`,
`gr_ctx_init_complex_qqbar`, `gr_ctx_init_real_ca`,
`gr_ctx_init_complex_ca` and related ones make these numbers available
to the generic polynomial and matrix algorithms ([](gr_domains#algebraic),
[](gr_domains#real-complex)). [](algorithms/exact/overview) compares the
representations in detail, and [](tutorials/exact) has worked examples.

## Generic rings {#generic}

The generic interface [`gr`](gr) makes all the ring types above, and
rings built from them, available through one set of functions over a
context `gr_ctx_t` ([](gr_domains)). Use it for rings without a module
of their own, for the types that exist only as generic rings (`mpn_mod`,
`nfloat`, `dfloat`, `decimal`, `padic_radix`), and for code that must
work over several rings. [](gr/getting-started) is the tutorial, and
[](guide/generic) discusses what generic code costs.

