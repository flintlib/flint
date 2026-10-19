---
title: gr_domains
header: gr.h
source: gr
description: The built-in rings and other structures of the generic interface, and their constructors.
category: Generic rings
order: 2
---

The constructors of the structures that come with FLINT for use with the
generic interface [`gr`](gr) are listed below: integers and rationals,
residue rings and finite fields, number fields and algebraic numbers,
real and complex numbers in ball, exact and floating-point
representations, and the rings built from other rings (polynomials,
series, matrices, fractions, quotients). Each constructor initialises a
`gr_ctx_t`, which must be freed with `gr_ctx_clear` after all its
elements; a ring built on a base ring stores a pointer to the base
context, which must stay alive (and must not be moved) until the outer
ring has been cleared. The constructors are declared in `gr.h`; some
rings are documented with the module that implements them, as listed in
each section.

Most constructors allocate nothing: they store the parameters in the
context and point it to a method table that is built once per type, so
creating a context costs a few nanoseconds and multimodular algorithms
can create one for each prime ([](algorithms/generic/method-tables#cheap)).
Exceptions are noted (finite fields, number fields, Calcium fields and
quotient rings precompute data).

The element types are the usual FLINT types, so a context can be used
with variables of the concrete type: an element of
`gr_ctx_init_fmpz` is an `fmpz`, of `gr_ctx_init_real_arb` an
`arb_struct`, of `gr_ctx_init_gr_poly` a `gr_poly_struct`.

| Structure | Constructors | Exact | Notes |
|---|---|---|---|
| $\Z$, $\mathbb{Q}$, $\Z[i]$ | `gr_ctx_init_fmpz`, `gr_ctx_init_fmpq`, `gr_ctx_init_fmpzi` | yes | |
| $\Z/n\Z$ | `gr_ctx_init_nmod`, `gr_ctx_init_mpn_mod`, `gr_ctx_init_fmpz_mod` | yes | word-size, $2^{64} \le n < 2^{1024}$, any $n$ |
| $\mathbb{F}_{p^d}$ | `gr_ctx_init_fq`, `gr_ctx_init_fq_nmod`, `gr_ctx_init_fq_zech` | yes | |
| number fields | `gr_ctx_init_nf` | yes | |
| $\mathbb{Q}_p$ | `gr_ctx_init_padic_radix` | at infinite precision | see [`padic_radix`](padic_radix) |
| $\overline{\mathbb{Q}} \cap \mathbb{R}$, $\overline{\mathbb{Q}}$ | `gr_ctx_init_real_qqbar`, `gr_ctx_init_complex_qqbar` | yes | |
| $\mathbb{R}$, $\mathbb{C}$ | `gr_ctx_init_real_ca`, `gr_ctx_init_complex_ca` | yes | some predicates undecidable |
| $\mathbb{R}$, $\mathbb{C}$ (balls) | `gr_ctx_init_real_arb`, `gr_ctx_init_complex_acb` | no | rigorous enclosures |
| floating point | `gr_ctx_init_real_float_arf`, `gr_ctx_init_complex_float_acf`, `nfloat_ctx_init`, `gr_ctx_init_dfloat`, `gr_ctx_init_mpf` | no | not rings |
| $R[x]$, $R[x_1, \ldots, x_n]$ | `gr_ctx_init_gr_poly`, `gr_ctx_init_gr_mpoly`, ... | as $R$ | |
| $R[[x]]$, $R[x]/(x^n)$ | `gr_series_ctx_init`, `gr_series_mod_ctx_init` | no, yes | |
| $R[x]/(m)$ | `gr_ctx_init_gr_poly_quotient` | as $R$ | |
| $\operatorname{Frac}(R)$, $R[i]$ | `gr_ctx_init_gr_fraction`, `gr_ctx_init_gr_complex` | as $R$ | |
| $R^n$, $M_{m \times n}(R)$ | `gr_ctx_init_vector_space_gr_vec`, `gr_ctx_init_matrix_space` | as $R$ | |
| groups | `gr_ctx_init_perm`, `gr_ctx_init_psl2z`, `gr_ctx_init_dirichlet_group` | yes | multiplicative |

## Integers and rational numbers {#integers}

::: function
void gr_ctx_init_fmpz(gr_ctx_t ctx)

Initialises `ctx` to the ring of integers $\Z$ with elements of type
`fmpz`.
:::

::: function
void gr_ctx_init_fmpq(gr_ctx_t ctx)

Initialises `ctx` to the field of rational numbers $\mathbb{Q}$ with elements of
type `fmpq`.
:::

::: function
void gr_ctx_init_fmpzi(gr_ctx_t ctx)

Initialises `ctx` to the ring of Gaussian integers $\Z[i]$ with elements
of type `fmpzi_struct` (see [`fmpzi`](fmpzi)).
:::

The ring of integers in radix representation (`gr_ctx_init_radix_integer`)
is documented in [`radix`](radix), and the field of $p$-adic numbers
$\mathbb{Q}_p$ in radix representation (`gr_ctx_init_padic_radix`,
declared in `padic_radix.h`) in [`padic_radix`](padic_radix).

## Residue rings {#residue-rings}

The contexts for $\Z/n\Z$ with a word-size modulus
(`gr_ctx_init_nmod`, `gr_ctx_init_nmod8`, `gr_ctx_init_nmod32`,
`gr_ctx_init_nmod_redc`, `gr_ctx_init_nmod_redc_fast`) are documented in
[`nmod`](nmod), and the context for moduli of 2 to 16 limbs with flat
limb-array elements, `gr_ctx_init_mpn_mod`, in [`mpn_mod`](mpn_mod)
(with the guide [](mpn_mod/usage)).
Whether the modulus is prime is decided by `gr_ctx_is_field` for word-size
moduli and must be declared with `gr_ctx_set_is_field` for the others;
all these rings can also pretend to be fields
(`gr_ctx_set_is_pretend_field`).

::: function
void gr_ctx_init_fmpz_mod(gr_ctx_t ctx, const fmpz_t n)

Initialises `ctx` to $\Z/n\Z$ for any modulus $n \ge 1$, with elements of
type `fmpz` reduced to $[0, n)$. The context allocates an
`fmpz_mod_ctx_t`. Primality of $n$ is not tested: `gr_ctx_is_field`
returns `T_UNKNOWN` until it is declared with `gr_ctx_set_is_field`.
For moduli between $2^{64}$ and $2^{1024}$, `gr_ctx_init_mpn_mod` is
faster.
:::

::: function
void _gr_ctx_init_fmpz_mod_from_ref(gr_ctx_t ctx, const void * fmod_ctx)

Initialises `ctx` to $\Z/n\Z$ using an existing `fmpz_mod_ctx_t`, which is
not copied and must outlive `ctx`. Used by `fmpz_mod_poly` and related
modules to call generic algorithms.
:::

## Finite fields {#finite-fields}

::: function
void gr_ctx_init_fq(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var)
void gr_ctx_init_fq_nmod(gr_ctx_t ctx, ulong p, slong d, const char * var)
void gr_ctx_init_fq_zech(gr_ctx_t ctx, ulong p, slong d, const char * var)

Initialise `ctx` to the finite field $\mathbb{F}_q$ with $q = p^d$,
using a Conway polynomial if one is in the database and otherwise an
irreducible polynomial of minimal weight (see
[](algorithms/finite-fields/defining-polynomials)).
The prime $p$ is assumed and not checked. The generator is printed as
`var`, or `a` if `var` is `NULL`. The element types are `fq_struct`,
`fq_nmod_struct` and `fq_zech_struct`; the `fq_zech` field stores
discrete logarithm tables of size $q$, so it is only practical for small
$q$ (and requires $q < 2^{\bits}$). The contexts allocate the underlying
`fq_ctx_t`. The context of [`fq_default`](fq_default), which chooses one
of these representations from the size of $q$, is itself a `gr`
context.
:::

::: function
void _gr_ctx_init_fq_from_ref(gr_ctx_t ctx, const void * fq_ctx)
void _gr_ctx_init_fq_nmod_from_ref(gr_ctx_t ctx, const void * fq_nmod_ctx)
void _gr_ctx_init_fq_zech_from_ref(gr_ctx_t ctx, const void * fq_zech_ctx)

Initialise `ctx` to the finite field described by an existing
`fq_ctx_t`, `fq_nmod_ctx_t` or `fq_zech_ctx_t`, which is not copied and
must outlive `ctx`. The constructors with a given modulus
(`gr_ctx_init_fq_nmod_modulus_nmod_poly` and similar) are documented in
the finite field modules.
:::

## Number fields and algebraic numbers {#algebraic}

::: function
void gr_ctx_init_nf(gr_ctx_t ctx, const fmpq_poly_t poly)
void gr_ctx_init_nf_fmpz_poly(gr_ctx_t ctx, const fmpz_poly_t poly)

Initialise `ctx` to the number field $\mathbb{Q}[a]/(f)$ defined by the
polynomial $f$ = `poly`, with elements of type `nf_elem_struct` (see
[`nf_elem`](nf_elem)). The
polynomial must be irreducible; this is not checked. The context
allocates an `nf_t` with precomputed data, so it is more expensive to
create than most contexts. Declared only if `fmpq_poly.h` is included
before `gr.h`.
:::

::: function
void _gr_ctx_init_nf_from_ref(gr_ctx_t ctx, const void * nfctx)

Initialises `ctx` to the number field described by an existing `nf_t`,
which is not copied and must outlive `ctx`.
:::

::: function
void gr_ctx_init_real_qqbar(gr_ctx_t ctx)
void gr_ctx_init_complex_qqbar(gr_ctx_t ctx)

Initialise `ctx` to the field of real algebraic numbers
$\overline{\mathbb{Q}} \cap \mathbb{R}$, respectively all algebraic numbers
$\overline{\mathbb{Q}}$, with elements of type `qqbar_struct` (a minimal
polynomial and an isolating ball; see [`qqbar`](qqbar) and
[the guide to qqbar](qqbar/using#gr)). Arithmetic is exact and equality is
decidable, but the degree of results can grow quickly.
:::

::: function
void _gr_ctx_qqbar_set_limits(gr_ctx_t ctx, slong deg_limit, slong bits_limit)

Limits the degrees of intermediate algebraic numbers of a *qqbar*
context to `deg_limit` and their bit sizes to about `bits_limit`; a
negative value removes the limit (the default). An operation that would
exceed a limit returns `GR_UNABLE`. The limits refer to the sizes of
resultants before factorisation (see `qqbar_binop_within_limits`), so
adding two numbers of degree 100 needs a degree limit of at least 10000.
Not all methods respect the limits.
:::

## Real and complex numbers {#real-complex}

FLINT has three kinds of real and complex numbers: balls (`arb`, `acb`),
which represent each number by a rigorous enclosure; exact numbers
(`ca`), for which arithmetic is exact but some predicates are
undecidable; and floating-point numbers, which are approximations
without error bounds and are not rings in the sense of the interface.
Generic algorithms over balls and exact numbers return rigorous results
([](algorithms/generic/inexact)).

::: function
void gr_ctx_init_real_arb(gr_ctx_t ctx, slong prec)
void gr_ctx_init_complex_acb(gr_ctx_t ctx, slong prec)

Initialise `ctx` to the field of real or complex numbers represented by
balls of type `arb_struct` and `acb_struct`, computing with a working
precision of `prec` bits. The precision can be changed with
`gr_ctx_set_real_prec`. Every result is an enclosure of the exact
result; predicates return `T_UNKNOWN` when the balls do not decide them.
:::

::: function
void gr_ctx_init_real_ca(gr_ctx_t ctx)
void gr_ctx_init_complex_ca(gr_ctx_t ctx)
void gr_ctx_init_real_algebraic_ca(gr_ctx_t ctx)
void gr_ctx_init_complex_algebraic_ca(gr_ctx_t ctx)

Initialise `ctx` to the field of real numbers, complex numbers, real
algebraic numbers or complex algebraic numbers, represented exactly by
elements of type `ca_struct` (see [`ca`](ca) and
[the guide to ca](ca/using#gr)). Operations whose result
is not in the field (for example $\sqrt{-1}$ in the reals, $\pi$ in the
algebraic numbers) return `GR_DOMAIN`, and predicates that cannot be
decided return `T_UNKNOWN`. The context allocates a `ca_ctx_t`; it is
not thread-safe (`gr_ctx_is_threadsafe` returns `T_FALSE`).
:::

::: function
void gr_ctx_init_complex_extended_ca(gr_ctx_t ctx)

As `gr_ctx_init_complex_ca`, also allowing the special values of `ca`:
signed and unsigned infinities, *Undefined* and *Unknown*. This is not a
field.
:::

::: function
void _gr_ctx_init_ca_from_ref(gr_ctx_t ctx, int which_ring, void * ca_ctx)

Initialises `ctx` to one of the `ca` structures, selected by
`which_ring` (`GR_CTX_RR_CA`, `GR_CTX_CC_CA`,
`GR_CTX_REAL_ALGEBRAIC_CA`, `GR_CTX_COMPLEX_ALGEBRAIC_CA` or
`GR_CTX_COMPLEX_EXTENDED_CA`), using an existing `ca_ctx_t` that is not
copied and must outlive `ctx`.
:::

::: function
void gr_ctx_ca_set_option(gr_ctx_t ctx, slong option, slong value)
slong gr_ctx_ca_get_option(gr_ctx_t ctx, slong option)

Set or get an option (such as `CA_OPT_PREC_LIMIT`) of the `ca_ctx_t` of
a `ca` context, as `ca_ctx_set_option` and `ca_ctx_get_option`.
:::

::: function
void gr_ctx_init_gr_complex(gr_ctx_t ctx, gr_ctx_t real_ctx)

Initialises `ctx` to the algebra $R[i]$ with $i^2 = -1$ over the ring $R$
given by `real_ctx`. An element $a + bi$ is stored as the pair $(a, b)$
of elements of $R$, contiguously.

::: details
Typically $R$ is an implementation of the real numbers or of a subring
of them (for example `nfloat`), and the result is the complex numbers or
a subring. The construction also makes sense over real vector spaces and
over general rings, in which $i$ is a formal element with $i^2 = -1$,
but then operations beyond basic arithmetic (absolute value, argument)
may not make sense.
:::
:::

## Floating-point numbers {#floating-point}

Floating-point numbers approximate the real and complex numbers but do
not form rings: arithmetic rounds, so the ring axioms hold only
approximately, and the context predicates (`gr_ctx_is_ring`,
`gr_ctx_is_field`) return `T_FALSE`, while
`gr_ctx_is_approx_commutative_ring` returns `T_TRUE`. Floating-point
contexts can still be used with most generic algorithms, which then
compute approximations without error bounds. Equality compares
floating-point values, with the convention that NaN is not equal to
itself. Correct rounding is guaranteed only for the basic operations
(add, sub, mul, div, sqrt) of real numbers. The fixed-precision types
[`nfloat`](nfloat) (`nfloat_ctx_init`, `nfloat_complex_ctx_init`),
[`dfloat`](dfloat) (`gr_ctx_init_dfloat`) and `decimal` (`gr_ctx_init_decfloat` and
related) are documented in their own modules and are much faster than
the contexts below at low and medium precision.

::: function
void gr_ctx_init_real_float_arf(gr_ctx_t ctx, slong prec)
void gr_ctx_init_complex_float_acf(gr_ctx_t ctx, slong prec)

Initialise `ctx` to binary floating-point numbers with elements of type
`arf_struct`, respectively complex floating-point numbers of type
`acf_struct`, with precision `prec` bits.
:::

::: function
void gr_ctx_init_mpf(gr_ctx_t ctx, slong prec)

Initialises `ctx` to the GMP floating-point numbers `mpf_t` with
precision `prec` bits. Only basic arithmetic, comparisons and a few
conversions are implemented; this context exists mainly for comparison
with the `mpf`-based algorithms of `fmpz_lll`.
:::

## Vectors and matrices {#vectors-matrices}

::: type
vector_ctx_t

The data of a vector context, stored in the `data` field of the
context: `base_ring` (a pointer to the context of the entries),
`all_sizes` (nonzero if vectors of every length are elements) and `n`
(the length otherwise).
:::

::: macro
VECTOR_CTX(ring_ctx)

The `vector_ctx_t` stored in a vector context.
:::

::: function
void gr_ctx_init_vector_gr_vec(gr_ctx_t ctx, gr_ctx_t base_ring)

Initialises `ctx` to the set of vectors of any length over `base_ring`,
with elements of type `gr_vec_struct` and componentwise arithmetic.
Operations on vectors of different lengths return `GR_DOMAIN`. This is
not a ring.
:::

::: function
void gr_ctx_init_vector_space_gr_vec(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)

Initialises `ctx` to the space $R^n$ of vectors of length $n$ over
`base_ring`, with elements of type `gr_vec_struct`. With componentwise
multiplication this is the ring $R^n$, a product of rings.
:::

::: type
matrix_ctx_t

The data of a matrix context: `base_ring`, `all_sizes`, `nrows` and
`ncols`.
:::

::: macro
MATRIX_CTX(ring_ctx)

The `matrix_ctx_t` stored in a matrix context.
:::

::: function
void gr_ctx_init_matrix_domain(gr_ctx_t ctx, gr_ctx_t base_ring)

Initialises `ctx` to the set of matrices of any shape over `base_ring`,
with elements of type `gr_mat_struct`. Operations on matrices of
incompatible shapes return `GR_DOMAIN`, as does `gr_one`, since there is
no single identity.
:::

::: function
void gr_ctx_init_matrix_space(gr_ctx_t ctx, gr_ctx_t base_ring, slong nrows, slong ncols)
void gr_ctx_init_matrix_ring(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)

Initialise `ctx` to the space of `nrows` by `ncols` matrices over
`base_ring`, respectively the ring $M_n(R)$ of $n \times n$ matrices,
with elements of type `gr_mat_struct`.
:::

::: function
truth_t gr_ctx_matrix_is_fixed_size(gr_ctx_t ctx)

Returns `T_TRUE` for a matrix space or ring of fixed shape and `T_FALSE`
for the domain of matrices of all shapes. `ctx` must be a matrix
context.
:::

## Polynomial rings {#polynomials}

::: function
void gr_ctx_init_fmpz_poly(gr_ctx_t ctx)
void gr_ctx_init_fmpq_poly(gr_ctx_t ctx)

Initialise `ctx` to the ring of univariate polynomials $\Z[x]$,
respectively $\mathbb{Q}[x]$, with elements of type `fmpz_poly_struct` and
`fmpq_poly_struct`, using the arithmetic of `fmpz_poly` and `fmpq_poly`.
:::

::: type
polynomial_ctx_t

The data of a `gr_poly` context: `base_ring`, `degree_limit` and `var`
(the name of the generator).
:::

::: macro
POLYNOMIAL_CTX(ring_ctx)
POLYNOMIAL_ELEM_CTX(ring_ctx)

The `polynomial_ctx_t` of a polynomial context and its `base_ring`.
:::

::: function
void gr_ctx_init_gr_poly(gr_ctx_t ctx, gr_ctx_t base_ring)

Initialises `ctx` to the ring $R[x]$ of densely represented univariate
polynomials over `base_ring`, with elements of type `gr_poly_struct`;
see [`gr_poly`](gr_poly). The generator is named `x` by default
(`gr_ctx_set_gen_name` changes it). Over inexact rings the length of a
polynomial is only an upper bound for its degree
([](algorithms/generic/inexact#weak-normalisation)).
:::

::: function
void gr_ctx_init_fmpz_mpoly(gr_ctx_t ctx, slong nvars, const ordering_t ord)
void gr_ctx_init_fmpq_mpoly(gr_ctx_t ctx, slong nvars, const ordering_t ord)
void gr_ctx_init_gr_mpoly(gr_ctx_t ctx, gr_ctx_t base_ring, slong nvars, const ordering_t ord)

Initialise `ctx` to the ring of sparse multivariate polynomials in
`nvars` variables with the monomial ordering `ord`, over $\Z$
(`fmpz_mpoly_struct` elements), over $\mathbb{Q}$ (`fmpq_mpoly_struct`), or over
`base_ring` (`gr_mpoly_struct`, see [`gr_mpoly`](gr_mpoly)). The
variables are named `x1`, `x2`, ... by default. Declared only if
`mpoly.h` is included before `gr.h`.
:::

Ore polynomial rings (`gr_ctx_init_gr_ore_poly`) are documented in
[`gr_ore_poly`](gr_ore_poly).

## Power series {#series}

The power series rings are documented in [`gr_series`](gr_series), with
the constructors `gr_series_ctx_init` and `gr_series_mod_ctx_init`.

::: function
void gr_ctx_init_gr_series(gr_ctx_t ctx, gr_ctx_t base_ring, slong prec)
void gr_ctx_init_series_mod_gr_poly(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)

Aliases of `gr_series_ctx_init` (power series $R[[x]]$ with an error
term $O(x^{\mathit{prec}})$ tracked per element) and
`gr_series_mod_ctx_init` (the exact ring $R[x]/(x^n)$).
:::

## Quotient rings of polynomial rings {#poly-quotient}

::: function
void gr_ctx_init_gr_poly_quotient(gr_ctx_t ctx, gr_ctx_t base, const gr_poly_t modulus)

Initialises `ctx` to the quotient ring $F[x] / (m)$ of the polynomial
ring over the commutative ring $F$ = `base` by the monic polynomial
$m$ = `modulus`. Elements have type `gr_poly_struct` and are kept
reduced modulo $m$. Throws an exception if $m$ is zero or not monic. The
context stores precomputed data for reduction modulo $m$, a
`gr_poly_preinv_t` whose representation is chosen by the method
`GR_METHOD_POLY_PREINV_SET` of the base ring, which is used for
multiplication and for powers.

::: details
If $F$ is a field and $m$ is irreducible, the quotient ring is a field,
but this is not assumed: `gr_ctx_is_field` returns `T_UNKNOWN` unless
the user declares $m$ irreducible with `gr_ctx_set_is_field` (it then
returns `T_TRUE` if $F$ is known to be a field).

The ring can also pretend to be a field (`gr_ctx_set_is_pretend_field`)
when $m$ is only conjecturally irreducible, provided that $F$ is a field
or pretends to be one. All operations remain correct as ring
operations, and an inversion of a zero divisor fails with `GR_UNABLE`
while recording a nontrivial monic factor of $m$ in a table in the
context. The first recorded factor can be retrieved as an element with
`gr_ctx_recover_zero_divisor`, and the whole table with the functions
below, from which a refined modulus can be obtained. A failure in a base
ring that pretends to be a field propagates as `GR_UNABLE` as well, and
the zero divisor is then recorded in the base ring. A modulus known to
be reducible can also be used without the pretence, in which case the
context represents a genuine (non-integral) quotient ring in which zero
divisors are simply not invertible (`GR_DOMAIN`).
:::
:::

::: function
const gr_poly_struct * gr_poly_quotient_ctx_modulus(gr_ctx_t ctx)
gr_ctx_struct * gr_poly_quotient_ctx_base(gr_ctx_t ctx)
slong gr_poly_quotient_ctx_degree(gr_ctx_t ctx)

Return the current modulus, the base ring, and the degree of the
modulus.
:::

::: function
int gr_poly_quotient_ctx_refine(gr_ctx_t ctx, const gr_poly_t new_modulus)

Replaces the modulus by `new_modulus`, which must be monic of degree at
least 1, and returns `GR_DOMAIN` otherwise. If the new modulus divides
the old one (for instance a factor recovered from a zero divisor),
existing elements remain valid representatives of their images in the
new ring: they are reduced by subsequent operations, so no conversion is
needed. Previous moduli are kept until the context is cleared, so that a
pointer obtained from `gr_poly_quotient_ctx_modulus` stays valid.
:::

::: function
ulong gr_poly_quotient_ctx_version(gr_ctx_t ctx)

Returns a counter that is incremented by each refinement.
:::

::: function
slong gr_poly_quotient_ctx_num_zero_divisors(gr_ctx_t ctx)
const gr_poly_struct * gr_poly_quotient_ctx_zero_divisor(gr_ctx_t ctx, slong i)
void gr_poly_quotient_ctx_clear_zero_divisors(gr_ctx_t ctx)

Access the table of nontrivial monic factors of the modulus discovered
by failed inversions while the ring pretends to be a field: the number
of entries, entry $i$, and clearing the table. Recording is thread-safe,
but these accessors must not be used concurrently with operations that
may record a new factor.
:::

::: function
int gr_poly_quotient_get_poly(gr_poly_t res, gr_srcptr x, gr_ctx_t ctx)
int gr_poly_quotient_set_poly(gr_ptr res, const gr_poly_t x, gr_ctx_t ctx)

Convert between elements of the quotient ring and polynomials over the
base ring; the result is reduced modulo the current modulus.
:::

## Fraction fields {#fractions}

::: function
void gr_ctx_init_gr_fraction(gr_ctx_t ctx, gr_ctx_t domain, int flags)

Initialises `ctx` to the fraction field of the integral domain `domain`.
A fraction is stored as two consecutive elements of the domain
(numerator, denominator). By default fractions are reduced by removing
the common content with `gr_gcd` and the denominator is normalised by
`gr_canonical_associate`. `flags` is a combination of
`GR_FRACTION_NO_REDUCTION` and `GR_FRACTION_STRONGLY_CANONICAL`, or 0.

::: details
The constructor does not verify that `domain` is an integral domain.
Over a ring with zero divisors the behaviour is undefined: nothing
detects a product of two nonzero denominators becoming zero.
:::
:::

::: macro
GR_FRACTION_NO_REDUCTION
GR_FRACTION_STRONGLY_CANONICAL

Flags for `gr_ctx_init_gr_fraction`. `GR_FRACTION_NO_REDUCTION`
disables the reduction of fractions, which allows fraction field
arithmetic over domains without a GCD or canonical associates, and can
be faster when GCDs are more expensive than letting the entries grow.
`GR_FRACTION_STRONGLY_CANONICAL` asserts that reduced fractions are
strongly canonical ($a/b = c/d$ if and only if $a = c$ and $b = d$),
which allows faster equality tests; this holds in unique factorisation
domains with a correct GCD and canonical associates.
:::

::: function
void gr_ctx_init_fmpz_mpoly_q(gr_ctx_t ctx, slong nvars, const ordering_t ord)

Initialises `ctx` to the field $\mathbb{Q}(x_1, \ldots, x_n)$ of multivariate
rational functions with $n$ = `nvars`, with elements of type
`fmpz_mpoly_q_struct` (see [`fmpz_mpoly_q`](fmpz_mpoly_q)). Declared
only if `fmpz_mpoly_q.h` is included before `gr.h`.
:::

::: function
void gr_ctx_init_fmpz_mpoly_q_factored(gr_ctx_t ctx, slong nvars, const ordering_t ord)

Initialises `ctx` to the same field $\mathbb{Q}(x_1, \ldots, x_n)$, with
elements of type `fmpz_mpoly_q_factored_struct` in which denominators
are kept fully factored: addition and multiplication then avoid
polynomial GCDs, while division requires factorisation. Elements can be
converted to and from the contexts of `gr_ctx_init_fmpz_mpoly_q` and
`gr_ctx_init_fmpz_mpoly` with the same number of variables and ordering
using `gr_set_other`.
:::

::: function
void gr_ctx_init_fmpz_mod_mpoly_q(gr_ctx_t ctx, slong nvars, const ordering_t ord, const fmpz_t mod)

Initialises `ctx` to the field $\mathbb{F}_p(x_1, \ldots, x_n)$ of
rational functions over $\mathbb{F}_p$ with $p$ = `mod`, with elements
of type `fmpz_mod_mpoly_q_struct`. The modulus must be prime; this is
not checked and a composite modulus gives undefined behaviour. Declared
only if `fmpz_mod_mpoly_q.h` is included before `gr.h`.
:::

## Groups {#groups}

These structures are multiplicative groups (`gr_ctx_is_multiplicative_group`
returns `T_TRUE`): `gr_mul`, `gr_inv`, `gr_div`, `gr_pow_si` and
`gr_one` are defined, while addition is not.

::: function
void gr_ctx_init_perm(gr_ctx_t ctx, ulong n)

Initialises `ctx` to the symmetric group $S_n$ of permutations of
$\{0, 1, \ldots, n - 1\}$. Elements are currently pointers to arrays (the
representation may change).
:::

::: function
void gr_ctx_init_psl2z(gr_ctx_t ctx)

Initialises `ctx` to the modular group $\mathrm{PSL}(2, \Z)$ with
elements of type `psl2z_struct`.
:::

::: function
int gr_ctx_init_dirichlet_group(gr_ctx_t ctx, ulong q)

Initialises `ctx` to the group $G_q$ of Dirichlet characters modulo
$q$, with elements of type `dirichlet_char_struct`. Returns `GR_DOMAIN`
if $q = 0$ and `GR_UNABLE` if $q$ has a prime factor larger than
$10^{16}$, which `dirichlet_group_init` does not support; in both cases
`ctx` is not initialised.
:::

## Symbolic expressions {#fexpr}

::: function
void gr_ctx_init_fexpr(gr_ctx_t ctx)

Initialises `ctx` to the set of symbolic expressions of type
`fexpr_struct` (see [`fexpr`](fexpr)). Arithmetic builds expression
trees without simplification. Elements of any ring can be converted to
expressions with `gr_get_fexpr`, and expressions can be evaluated in a
ring with `gr_set_fexpr` or `gr_set_other`.
:::

## Random rings {#random}

These constructors are used by the test suite to run generic algorithms
over many rings.

::: function
void gr_ctx_init_random(gr_ctx_t ctx, flint_rand_t state)
void gr_ctx_init_random_commutative_ring(gr_ctx_t ctx, flint_rand_t state)
void gr_ctx_init_random_field(gr_ctx_t ctx, flint_rand_t state)
void gr_ctx_init_random_finite_field(gr_ctx_t ctx, flint_rand_t state)

Initialise `ctx` to a random ring, commutative ring, field or finite
field. The rings are drawn from the built-in base rings ($\Z$, $\mathbb{Q}$,
$\Z[i]$, $\Z/n\Z$ in its various representations, finite fields, number
fields, balls, Calcium fields) and from rings built over some of them
(polynomials, multivariate polynomials, series, vectors, matrices and
debug wrappers). The base rings of the composite rings are kept in
static storage, so only `ctx` is cleared by the caller.

::: details
Currently `gr_ctx_init_random` returns, with probability 1/2, the ring of
power series of precision 1 over the zero ring $\Z/1\Z$, and only
otherwise draws from the list above.
:::
:::

::: function
void gr_ctx_init_random_poly(gr_ctx_t ctx, flint_rand_t state)
void gr_ctx_init_random_mpoly(gr_ctx_t ctx, flint_rand_t state)
void gr_ctx_init_random_series(gr_ctx_t ctx, flint_rand_t state)

Initialise `ctx` to a random ring of univariate polynomials, multivariate
polynomials or power series over a random base ring. The random series
rings (also those among the composite rings of `gr_ctx_init_random`)
have precision 1 to 5, and only with probability 1/20 precision 0, in
which every element is $O(x^0)$.
:::

## Debugging {#debug}

::: function
void gr_ctx_init_debug(gr_ctx_t ctx, gr_ctx_t elem_ctx, int flags, double unable_probability)

Initialises `ctx` to a wrapper around the ring `elem_ctx` that computes
as `elem_ctx` but adds checks. Only a few methods are wrapped (memory
management, printing, `gr_zero`, `gr_one`, `gr_is_zero`, `gr_is_one`,
`gr_equal`, `gr_set`, `gr_neg`, `gr_add`, `gr_sub`, `gr_mul`, `gr_div`,
`gr_inv`, `gr_is_invertible`, `gr_sqrt`, `gr_rsqrt`, `gr_pow`); the
others use the generic defaults. If `unable_probability` is positive,
the wrapped methods randomly return `GR_UNABLE` and predicates randomly
return `T_UNKNOWN` with that probability, which tests the error handling
of algorithms. `flags` is a combination of the `GR_DEBUG_*` flags.
:::

::: macro
GR_DEBUG_WRAP
GR_DEBUG_VERBOSE
GR_DEBUG_CHECK_ALWAYS_ABLE
GR_DEBUG_TIMING

Flags for `gr_ctx_init_debug`. `GR_DEBUG_WRAP` stores each element in a
box with a pointer to its context, which detects elements used with the
wrong context and can reveal bugs hidden when elements are stored
shallowly. `GR_DEBUG_VERBOSE` prints every operation with its inputs,
output and status. `GR_DEBUG_CHECK_ALWAYS_ABLE` throws an exception
when a wrapped method returns `GR_UNABLE` or `T_UNKNOWN`.
`GR_DEBUG_TIMING` prints the time taken by each operation.
:::

## Coercions {#coercions}

::: function
int gr_ctx_cmp_coercion(gr_ctx_t ctx1, gr_ctx_t ctx2)

Returns 1 if elements should be coerced into `ctx1` and $-1$ if into
`ctx2`, when an operation combines elements of both. The rule compares
the implementation codes (`gr_which_structure`, in which structures
built on other structures come later), recursing into the base rings of
two `gr_poly` or two matrix contexts; for equal codes it returns 1.
:::
