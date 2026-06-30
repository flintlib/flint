.. _gr-ore-poly:

**gr_ore_poly.h** -- dense univariate Ore polynomials over generic rings
===============================================================================

.. note::

    This module is under construction.

A :type:`gr_ore_poly_t` represents a univariate Ore polynomial `L \in R[D]`
implemented as a dense array of coefficients in a generic ring *R*.
The choice of Ore algebra structure (e.g. with `D` being the standard
derivative or Euler derivative) is stored in the context object
:type:`gr_ore_poly_ctx_t`.

Most functions are provided in two versions: an underscore method which
operates directly on pre-allocated arrays of coefficients and generally
has some restrictions (often requiring the lengths to be nonzero
and not supporting aliasing of the input and output arrays),
and a non-underscore method which performs automatic memory
management and handles degenerate cases.

Ore algebra types
--------------------------------------------------------------------------------

.. type:: ore_algebra_t

    Represents one of the following supported Ore algebra types:

    .. macro:: ORE_ALGEBRA_CUSTOM

        Custom Ore polynomials.

    .. macro:: ORE_ALGEBRA_COMMUTATIVE

        Standard polynomials.

    .. macro:: ORE_ALGEBRA_DERIVATIVE

        Linear differential operators in the standard derivative.

        The endomorphism `\sigma` is the identity, and the `\sigma`-derivation
        `\delta` is the derivative `\frac{d}{dx}` with respect to a generator
        `x` of the base ring.

    .. macro:: ORE_ALGEBRA_EULER_DERIVATIVE

        Linear differential operators in the Euler derivative.

        The endomorphism `\sigma` is the identity, and the `\sigma`-derivation
        `\delta` is the Euler derivative `x\cdot\frac{d}{dx}` with respect to a
        generator `x` of the base ring.

    .. macro:: ORE_ALGEBRA_FORWARD_SHIFT

        Linear difference operators in the standard forward shift.

        The endomorphism `\sigma` is the shift `x \mapsto x + 1` with respect
        to a generator `x` of the base ring, and the `\sigma`-derivation
        `\delta` is the zero map.

    .. macro:: ORE_ALGEBRA_FORWARD_DIFFERENCE

        Linear difference operator in the forward finite difference operator.

        The endomorphism `\sigma` is the shift `x \mapsto x + 1` with respect
        to a generator `x` of the base ring, and the `\sigma`-derivation
        `\delta` maps `x \mapsto 1`.

    .. macro:: ORE_ALGEBRA_BACKWARD_SHIFT

        Linear difference operators in the standard backward shift.

    .. macro:: ORE_ALGEBRA_BACKWARD_DIFFERENCE

        Linear difference operator in the backward finite difference operator.

    .. macro:: ORE_ALGEBRA_Q_SHIFT

        Linear q-difference operators.

    .. macro:: ORE_ALGEBRA_MAHLER

        Linear Mahler operators.

    .. macro:: ORE_ALGEBRA_FROBENIUS

        Ore polynomials over a field twisted by the Frobenius endomorphism.

.. function:: ore_algebra_t ore_algebra_randtest(flint_rand_t state)

    Return a random Ore algebra type.

Type compatibility
-------------------------------------------------------------------------------

The ``gr_ore_poly`` type has the same data layout as ``gr_poly``.
Methods of ``gr_poly`` can therefore be used for linear and container
operations on a ``gr_ore_poly``, given that one is careful about providing
the right context object.

Weak normalization
-------------------------------------------------------------------------------

A :type:`gr_ore_poly_t` is always normalised by removing leading zeros.
For rings without decidable equality (e.g. rings with inexact
representation), only coefficients that are provably zero will be
removed, and there can thus be spurious leading zeros in the
internal representation.
Methods that depend on knowing the exact degree of a polynomial
will act appropriately, typically by returning ``GR_UNABLE``
when it is unknown whether the leading stored coefficient is nonzero.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: gr_ore_poly_struct

.. type:: gr_ore_poly_t

    Contains a pointer to an array of coefficients (``coeffs``), the used
    length (``length``), and the allocated size of the array (``alloc``).

    A ``gr_ore_poly_t`` is defined as an array of length one of type
    ``gr_ore_poly_struct``, permitting a ``gr_ore_poly_t`` to
    be passed by reference.

Context object methods
-------------------------------------------------------------------------------

.. function:: void gr_ore_poly_ctx_init(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, slong base_var, const ore_algebra_t which_algebra)

    Initializes *ctx* to a ring of densely represented Ore polynomials over
    the given *base_ring*, with the choice of Ore algebra structure given by
    *which_algebra*. The Ore algebra structure may refer to a distinguished
    generator of *base_ring*; this will be the generator of index
    *base_var*.

    This function can be used with all Ore algebra types for which no more
    specific initialization function is listed below.

.. function:: int gr_ore_poly_ctx_init_q_shift(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, slong base_var, gr_srcptr q)
              int gr_ore_poly_ctx_init_q_difference(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, slong base_var, gr_srcptr q)
              int gr_ore_poly_ctx_init_mahler(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, slong base_var, long mahler_base)

    Like :func:`gr_ore_poly_ctx_init` for predefined Ore polynomial types where
    `\sigma` and `\delta` depend on parameters.

.. function:: void * gr_ore_poly_ctx_data_ptr(gr_ore_poly_ctx_t ctx)

.. function:: void gr_ore_poly_ctx_init_custom(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, const gr_ore_poly_sigma_delta_t sigma_delta, void * ore_data)

    Initializes *ctx* to a ring of densely represented Ore polynomials over
    the given *base_ring*, with a custom Ore algebra structure specified by a
    pointer *sigma_delta* to an implementation of
    :func:`gr_ore_poly_sigma_delta`.
    The *ore_data* argument is accessible to *sigma_delta* as
    ``gr_ore_poly_ctx_data_ptr(ctx)``.

.. function:: void gr_ore_poly_ctx_init_randtest(gr_ore_poly_ctx_t ctx, flint_rand_t state, gr_ctx_t base_ring)

    Initializes *ctx* with a random Ore algebra structure.

.. function:: void gr_ore_poly_ctx_init_randtest2(gr_ctx_t base_ring, gr_ore_poly_ctx_t ctx, flint_rand_t state)

    Initializes *ctx* with a random Ore algebra structure over a random base
    ring.

.. function:: void gr_ore_poly_ctx_clear(gr_ore_poly_ctx_t ctx)

    Clears the context object *ctx*.

The following methods implement parts of the standard interface
for *gr* context objects.

.. function:: int _gr_ore_poly_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
              int _gr_ore_poly_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)

    Sets the name of the generator to the string in *s*, respectively the
    first string in *s*.

.. function:: int gr_ore_poly_ctx_write(gr_stream_t out, gr_ore_poly_ctx_t ctx)
              truth_t gr_ore_poly_ctx_is_ring(gr_ore_poly_ctx_t ctx)
              truth_t gr_ore_poly_ctx_is_zero_ring(gr_ore_poly_ctx_t ctx)
              truth_t gr_ore_poly_ctx_is_commutative_ring(gr_ore_poly_ctx_t ctx)
              truth_t gr_ore_poly_ctx_is_integral_domain(gr_ore_poly_ctx_t ctx)
              truth_t gr_ore_poly_ctx_is_threadsafe(gr_ore_poly_ctx_t ctx)

Memory management
-------------------------------------------------------------------------------

.. function:: void gr_ore_poly_init(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

.. function:: void gr_ore_poly_init2(gr_ore_poly_t poly, slong len, gr_ore_poly_ctx_t ctx)

.. function:: void gr_ore_poly_clear(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

.. function:: gr_ptr gr_ore_poly_coeff_ptr(gr_ore_poly_t poly, slong i, gr_ore_poly_ctx_t ctx)
              gr_srcptr gr_ore_poly_coeff_srcptr(const gr_ore_poly_t poly, slong i, gr_ore_poly_ctx_t ctx)

.. function:: slong gr_ore_poly_length(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

.. function:: void gr_ore_poly_swap(gr_ore_poly_t poly1, gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx)

.. function:: void gr_ore_poly_fit_length(gr_ore_poly_t poly, slong len, gr_ore_poly_ctx_t ctx)

.. function:: void _gr_ore_poly_set_length(gr_ore_poly_t poly, slong len, gr_ore_poly_ctx_t ctx)

Basic manipulation
-------------------------------------------------------------------------------

.. function:: void _gr_ore_poly_normalise(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

.. function:: int gr_ore_poly_set(gr_ore_poly_t res, const gr_ore_poly_t src, gr_ore_poly_ctx_t ctx)

.. function:: int gr_ore_poly_truncate(gr_ore_poly_t res, const gr_ore_poly_t poly, slong newlen, gr_ore_poly_ctx_t ctx)

.. function:: int gr_ore_poly_zero(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_one(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_neg_one(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_gen(gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

.. function:: int gr_ore_poly_write(gr_stream_t out, const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
              int _gr_ore_poly_write(gr_stream_t out, gr_srcptr poly, slong n, gr_ore_poly_ctx_t ctx)
              int _gr_ore_poly_get_str(char ** res, const gr_ore_poly_t f, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_get_str(char ** res, const gr_ore_poly_t f, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_print(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

.. function:: int _gr_ore_poly_set_str(gr_ptr res, const char * s, slong len, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_set_str(gr_ore_poly_t res, const char * s, gr_ore_poly_ctx_t ctx)

    Parse Ore polynomial from an expression string, assuming the name of the
    generator is the one set in *ctx*. The underscore method zero-pads the
    result if the length of the parsed polynomial is shorter than *len*,
    and returns ``GR_UNABLE`` if the length of the parsed polynomial exceeds
    *len*. Intermediate terms are allowed to be longer than *len*.

.. function:: int gr_ore_poly_randtest(gr_ore_poly_t poly, flint_rand_t state, slong len, gr_ore_poly_ctx_t ctx)

.. function:: truth_t _gr_ore_poly_equal(gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx)
              truth_t gr_ore_poly_equal(const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx)

.. function:: truth_t gr_ore_poly_is_zero(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
              truth_t gr_ore_poly_is_one(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)
              truth_t gr_ore_poly_is_gen(const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

.. function:: int gr_ore_poly_set_si(gr_ore_poly_t poly, slong c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_set_ui(gr_ore_poly_t poly, ulong c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_set_fmpz(gr_ore_poly_t poly, const fmpz_t c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_set_fmpq(gr_ore_poly_t poly, const fmpq_t c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_set_other(gr_ore_poly_t poly, gr_srcptr x, gr_ctx_t x_ctx, gr_ore_poly_ctx_t ctx)

Action
-------------------------------------------------------------------------------

.. function:: int gr_ore_poly_sigma(gr_ptr res, gr_srcptr a, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_delta(gr_ptr res, gr_srcptr a, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_sigma_delta(gr_ptr sigma, gr_ptr delta, gr_srcptr a, gr_ore_poly_ctx_t ctx)

    Compute `\sigma(a)`, `\delta(a)`, or both, where *a* is an element of the base ring.
    In the *sigma_delta* variant, the output variables *sigma* or *delta* can be
    ``NULL``.

.. type:: gr_ore_poly_sigma_delta_t

    A pointer to a function with the same specification as
    :func:`gr_ore_poly_sigma_delta`.

.. function:: int gr_ore_poly_apply_custom(gr_ptr res, const gr_ore_poly_t P, gr_srcptr f, gr_srcptr d1, gr_ore_poly_ctx_t ctx)

    Sets *res* to the result of applying *P* to the base ring element *f*, where
    the generator `D` acts by `g \mapsto \sigma(g) \cdot d1 + \delta(g)` for the
    given value ``d1`` of `D(1)`. Any ``d1`` defines a valid action. Aliasing of
    *res* with *f* or *d1* is allowed.

.. function:: int gr_ore_poly_apply(gr_ptr res, const gr_ore_poly_t P, gr_srcptr f, gr_ore_poly_ctx_t ctx)

    As :func:`gr_ore_poly_apply_custom`, but with the standard `D(1)` for the
    algebra type (so derivative operators differentiate, shift operators shift,
    etc.). Returns ``GR_UNABLE`` when no standard action exists
    (`ORE_ALGEBRA_COMMUTATIVE`, `ORE_ALGEBRA_CUSTOM`) or it cannot be computed.

Conversions
-------------------------------------------------------------------------------

These functions rewrite an operator from one Ore algebra into another. The
underscore versions act on raw coefficient arrays over the *coefficient* ring
*ctx* (the base ring of the Ore ring, not the Ore ring itself); *x* denotes its
generator of index *var*, that is, entry *var* of :func:`gr_gens`. Unless stated
otherwise, *res* has the same length *len* as *op* and must not alias it.

A conversion is an exact rewriting up to a unit factor. Since *x* and the
forward shift `S` need not act invertibly in the chosen representation, the
result may be the input premultiplied **on the left** by a power of *x* or `S`,
returned in an output parameter; left multiplication does not change the
solution space of the operator. Conversions that substitute into or
Taylor-shift the variable require a :macro:`GR_CTX_GR_POLY` base ring and return
``GR_UNABLE`` otherwise.

The first pair rewrites between the derivative `D = d/dx` and the Euler
derivative `\theta = x D`.

.. function:: int _gr_ore_poly_ddx_to_euler(gr_ptr res, gr_srcptr op, slong len, slong var, gr_ctx_t ctx)
              int _gr_ore_poly_euler_to_ddx(gr_ptr res, gr_srcptr op, slong len, slong var, gr_ctx_t ctx)

    Rewrite *op* between `D` and `\theta`. As `x` may not be invertible,
    :func:`_gr_ore_poly_ddx_to_euler` returns `x^{len-1} \cdot \text{op}`; the
    reverse rewriting needs no factor of `x`.

The next pair converts among the four shift/difference algebras: the forward and
backward shifts `S`, `S^{-1}` and the forward and backward differences `S-1`,
`1-S^{-1}`.

.. function:: int _gr_ore_poly_shift_convert(gr_ptr res, slong * epow, gr_srcptr op, slong len, ore_algebra_t src_alg, ore_algebra_t dst_alg, slong var, gr_ctx_t ctx)

    Convert *op* from *src_alg* to *dst_alg*. Each generator is a degree-one
    Laurent polynomial in `S`, so on output ``*epow`` `= p` and *res* satisfy
    `\text{res}_{dst} = S^{p} \cdot \text{op}_{src}`. A nonzero `p` (hence a
    polynomial base ring) is needed only when crossing between the forward side
    `S`, `S-1` and the backward side `S^{-1}`, `1-S^{-1}`, where `p = \mp(len-1)`.
    Converting between the two difference types uses a single binomial transform.
    The variable is the base ring generator of index *var*; only ``var == 0``
    over a univariate polynomial ring is implemented, and other indices return
    ``GR_UNABLE``. Returns ``GR_DOMAIN`` if either algebra is not a
    shift/difference type.

The remaining functions bridge differential and shift operators through the
generating-series isomorphism `\theta \mapsto n`, `x \mapsto S^{-1}`: for
`f = \sum_n a_n x^n`, the Euler derivative `\theta = x d/dx` acts on `(a_n)` as
multiplication by `n` and `x` as the backward shift `S^{-1}`. All require a
`GR_CTX_GR_POLY` base ring and exchange the operator order with the coefficient
degree.

.. function:: int _gr_ore_poly_euler_to_backshift_univar(gr_ptr res, slong reslen, gr_srcptr op, slong len, gr_ctx_t ctx)
              int _gr_ore_poly_backshift_to_euler_univar(gr_ptr res, slong reslen, gr_srcptr op, slong len, gr_ctx_t ctx)

    The two inverse rewritings of the isomorphism between the Euler operator and
    the backward shift `S^{-1}`. The caller allocates *res* to *reslen*, one more
    than the largest coefficient degree of *op*.

.. function:: int gr_ore_poly_differential_to_shift(gr_ore_poly_t res, slong * power, const gr_ore_poly_t op, gr_ore_poly_ctx_t res_ctx, gr_ore_poly_ctx_t op_ctx)
              int gr_ore_poly_shift_to_differential(gr_ore_poly_t res, slong * power, const gr_ore_poly_t op, gr_ore_poly_ctx_t res_ctx, gr_ore_poly_ctx_t op_ctx)

    Convert a whole operator between a differential algebra
    (`ORE_ALGEBRA_DERIVATIVE` or `ORE_ALGEBRA_EULER_DERIVATIVE`) and a
    shift/difference algebra, routing through the Euler and backward-shift pivots
    above. On output ``*power`` `= p` is the left power absorbed (`S^p` for
    :func:`gr_ore_poly_differential_to_shift`, `x^p` for
    :func:`gr_ore_poly_shift_to_differential`), so `\text{op}(f) = 0` and
    `\text{res} = 0` have the same solutions up to the index shift by `p`. Both
    contexts must share a `GR_CTX_GR_POLY` base ring (else ``GR_UNABLE``); returns
    ``GR_DOMAIN`` if a context algebra is outside its family.

.. function:: int gr_ore_poly_convert(gr_ore_poly_t res, slong * power, const gr_ore_poly_t op, gr_ore_poly_ctx_t res_ctx, gr_ore_poly_ctx_t op_ctx)

    Convert *op* from the algebra of *op_ctx* to that of *res_ctx*, provided both
    belong to the same family: the differential family
    (`ORE_ALGEBRA_DERIVATIVE`, `ORE_ALGEBRA_EULER_DERIVATIVE`) or the
    shift/difference family (the four shift and difference types). This
    dispatches to the appropriate conversion above and handles memory
    management. On output ``*power`` `= p` is the left power absorbed, **but its
    meaning depends on the pair of algebras**: for a conversion within the
    differential family `\text{res} = x^{p} \cdot \text{op}` (a power of the base
    generator `x`), whereas for one within the shift/difference family
    `\text{res} = S^{p} \cdot \text{op}` (a power of the forward shift `S`). A
    conversion to the same algebra is the identity, with `p = 0`. Returns
    ``GR_UNABLE`` for an unsupported pair: one that crosses the
    differential/shift boundary (use :func:`gr_ore_poly_differential_to_shift` or
    :func:`gr_ore_poly_shift_to_differential` for those) or that involves an
    algebra outside the two families, such as `ORE_ALGEBRA_Q_SHIFT`,
    `ORE_ALGEBRA_MAHLER`, `ORE_ALGEBRA_FROBENIUS` or `ORE_ALGEBRA_COMMUTATIVE`.

Arithmetic
-------------------------------------------------------------------------------

.. function:: int gr_ore_poly_neg(gr_ore_poly_t res, const gr_ore_poly_t src, gr_ore_poly_ctx_t ctx)

.. function:: int _gr_ore_poly_add(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_add(gr_ore_poly_t res, const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx)

.. function:: int _gr_ore_poly_sub(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_sub(gr_ore_poly_t res, const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx)

.. function:: int gr_ore_poly_add_ui(gr_ore_poly_t res, const gr_ore_poly_t poly, ulong c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_add_si(gr_ore_poly_t res, const gr_ore_poly_t poly, slong c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_add_fmpz(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpz c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_add_fmpq(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpq c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_add_other(gr_ore_poly_t res, const gr_ore_poly_t poly, gr_srcptr x, gr_ctx_t x_ctx, gr_ore_poly_ctx_t ctx)

    Sets *res* to *poly* plus the scalar *c* which must be
    an element of or coercible to the coefficient ring.

.. function:: int gr_ore_poly_sub_ui(gr_ore_poly_t res, const gr_ore_poly_t poly, ulong c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_sub_si(gr_ore_poly_t res, const gr_ore_poly_t poly, slong c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_sub_fmpz(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpz c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_sub_fmpq(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpq c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_sub_other(gr_ore_poly_t res, const gr_ore_poly_t poly, gr_srcptr x, gr_ctx_t x_ctx, gr_ore_poly_ctx_t ctx)

    Sets *res* to *poly* minus *c* which must be
    an element of or coercible to the coefficient ring.

.. function:: int gr_ore_poly_mul_ui(gr_ore_poly_t res, const gr_ore_poly_t poly, ulong c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_mul_si(gr_ore_poly_t res, const gr_ore_poly_t poly, slong c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_mul_fmpz(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpz c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_mul_fmpq(gr_ore_poly_t res, const gr_ore_poly_t poly, const fmpq c, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_other_mul(gr_ore_poly_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

    Sets *res* to *poly* multiplied by *c* (or *x* multiplied by *poly*)
    which must be an element of or coercible to the coefficient ring.

.. function:: int _gr_ore_poly_lmul_gen(gr_ptr res, gr_srcptr poly, slong len, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_lmul_gen(gr_ore_poly_t res, const gr_ore_poly_t poly, gr_ore_poly_ctx_t ctx)

    Sets *res* to the result of the left multiplication:
    `D \cdot poly = \sigma(poly) \cdot D + \delta(poly)`.
    The underscore method assumes *len != 0*.

.. function:: int _gr_ore_poly_mul(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_mul(gr_ore_poly_t res, const gr_ore_poly_t poly1, const gr_ore_poly_t poly2, gr_ore_poly_ctx_t ctx)

    Sets *res* to *poly1* multiplied by *poly2*
    which must be two Ore polynomials in the Ore algebra *ctx*.
    The underscore method assumes *res != poly1*, *res != poly2* (no aliasing) and *len1* != 0, *len2* != 0.

.. function:: int _gr_ore_poly_divrem(gr_ptr Q, gr_ptr R, gr_srcptr U, slong lenU, gr_srcptr V, slong lenV, gr_ore_poly_ctx_t ctx)
              int gr_ore_poly_divrem(gr_ore_poly_t Q, gr_ore_poly_t R, const gr_ore_poly_t U, gr_ore_poly_t V, gr_ore_poly_ctx_t ctx)

    Sets *(Q, R)* to the unique pair such that `U = QV + R` and `ord(R) < ord(V)`.

.. function:: int gr_ore_poly_div(gr_ore_poly_t Q, const gr_ore_poly_t U, gr_ore_poly_t V, gr_ore_poly_ctx_t ctx)

    Version of the divrem function which outputs only the quotient.

.. function:: int gr_ore_poly_rem(gr_ore_poly_t R, const gr_ore_poly_t U, gr_ore_poly_t V, gr_ore_poly_ctx_t ctx)

    Version of the divrem function which outputs only the remainder.

.. raw:: latex

    \newpage

