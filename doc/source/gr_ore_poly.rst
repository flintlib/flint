.. _gr-ore-poly:

**gr_ore_poly.h** -- dense univariate Ore polynomials over generic rings
===============================================================================

.. note::

    This module is under construction. Functionality is currently limited to
    memory management, additive arithmetic, and multiplication from the left
    by an element of the base ring.

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

    Initializes ``ctx`` to a ring of densely represented Ore polynomials over
    the given ``base_ring``, with the choice of Ore algebra structure given by
    ``which_algebra``. The Ore algebra structure may refer to a distinguished
    generator of ``base_ring``; this will be the generator of index
    ``base_var``.

    This function can be used with all Ore algebra types for which no more
    specific initialization function is listed below.

.. function:: int gr_ore_poly_ctx_init_q_shift(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, slong base_var, gr_srcptr q)
              int gr_ore_poly_ctx_init_q_difference(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, slong base_var, gr_srcptr q)
              int gr_ore_poly_ctx_init_mahler(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, slong base_var, long mahler_base)

    Like :func:`gr_ore_poly_ctx_init` for predefined Ore polynomial types where
    `\sigma` and `\delta` depend on parameters.

.. function:: void * gr_ore_poly_ctx_data_ptr(gr_ore_poly_ctx_t ctx)

.. function:: void gr_ore_poly_ctx_init_custom(gr_ore_poly_ctx_t ctx, gr_ctx_t base_ring, const gr_ore_poly_sigma_delta_t sigma_delta, void * ore_data)

    Initializes ``ctx`` to a ring of densely represented Ore polynomials over
    the given ``base_ring``, with a custom Ore algebra structure specified by a
    pointer ``sigma_delta`` to an implementation of
    :func:`gr_ore_poly_sigma_delta`.
    The ``ore_data`` argument is accessible to ``sigma_delta`` as
    ``gr_ore_poly_ctx_data_ptr(ctx)``.

.. function:: void gr_ore_poly_ctx_init_randtest(gr_ore_poly_ctx_t ctx, flint_rand_t state, gr_ctx_t base_ring)

    Initializes ``ctx`` with a random Ore algebra structure.

.. function:: void gr_ore_poly_ctx_init_randtest2(gr_ctx_t base_ring, gr_ore_poly_ctx_t ctx, flint_rand_t state)

    Initializes ``ctx`` with a random Ore algebra structure over a random base
    ring.

.. function:: void gr_ore_poly_ctx_clear(gr_ore_poly_ctx_t ctx)

    Clears the context object ``ctx``.

The following methods implement parts of the standard interface
for ``gr`` context objects.

.. function:: int _gr_ore_poly_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
              int _gr_ore_poly_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)

    Sets the name of the generator to the string in ``s``, respectively the
    first string in ``s``.

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

    Compute *σ(a)*,  *δ(a)*, or both, where *a* is an element of the base ring.
    In the *sigma_delta* variant, the output variables *sigma* or *delta* can be
    `NULL`.

.. type:: gr_ore_poly_sigma_delta_t

    A pointer to a function with the same specification as
    :func:`gr_ore_poly_sigma_delta`.

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

.. raw:: latex

    \newpage

