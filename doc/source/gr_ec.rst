.. _gr-ec:

**gr_ec.h** -- elliptic curves over generic rings
===============================================================================

.. note::

    This module is under construction. Functionality is currently limited to
    curve invariants and basic point arithmetic, and the arithmetic is only
    implemented for curves given in short Weierstrass form
    (see :ref:`the section on curve models <gr-ec-models>`).

A :type:`gr_ec_ctx_t` represents an elliptic curve `E` over a generic
commutative ring *R*, given by a general (long) Weierstrass equation

.. math::

    E : y^2 + a_1 x y + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6

with coefficients `a_1, a_2, a_3, a_4, a_6 \in R`.
A :type:`gr_ec_point_t` represents a point of `E`, stored in homogeneous
projective coordinates over the same ring.

The general Weierstrass equation is used as the interface for all curves
since it is the form needed over rings of residue characteristic 2 and 3,
where the short Weierstrass form `y^2 = x^3 + a_4 x + a_6` is not general
enough.

Since the base ring is not assumed to be a field, no division is performed
in the point arithmetic: the projective coordinates of a point are only
inverted on explicit request, for example by :func:`gr_ec_point_get_affine`
or :func:`gr_ec_point_normalize`.

Like the rest of the *gr* interface, functions return a status flag which
is ``GR_SUCCESS`` on success, ``GR_DOMAIN`` if the result does not exist,
and ``GR_UNABLE`` if the implementation is unable to perform the operation.
The ``GR_UNABLE`` flag is in particular returned when a decision about a
coordinate (typically, whether it is zero or invertible) cannot be made in
the base ring; this happens for rings with inexact representation, and for
rings where such predicates are not implemented.

Unless otherwise stated, aliasing between input and output objects is
allowed.

.. _gr-ec-models:

Curve models
-------------------------------------------------------------------------------

The point arithmetic is dispatched on the *model* of the curve, which is
determined when the context object is created and which is a property of the
stored coefficients, not of the isomorphism class of the curve.

.. type:: gr_ec_model_t

    Represents one of the following supported curve models:

    .. macro:: GR_EC_LONG_WEIERSTRASS

        A curve `y^2 + a_1 x y + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6`
        with at least one of `a_1, a_2, a_3` not known to be zero.
        This model is valid over any base ring, in particular over rings of
        residue characteristic 2 and 3.

        The group law in this model is not yet implemented:
        :func:`gr_ec_point_add`, :func:`gr_ec_point_sub`,
        :func:`gr_ec_point_dbl` and the scalar multiplication functions
        currently return ``GR_UNABLE``. Everything that does not require
        the addition formulas is supported, including the curve invariants,
        :func:`gr_ec_point_neg`, :func:`gr_ec_point_is_on_curve` and,
        when 2 is invertible in the base ring,
        :func:`gr_ec_point_lift_x`.

    .. macro:: GR_EC_SHORT_WEIERSTRASS

        A curve `y^2 = x^3 + a_4 x + a_6`, that is, a curve with
        `a_1 = a_2 = a_3 = 0`. This is the model for which arithmetic is
        currently implemented. Every elliptic curve over a ring in which 6
        is invertible -- in particular over a field of characteristic 0 or
        of characteristic `p > 3` -- is isomorphic to a curve in this model.

    .. macro:: GR_EC_NUM_MODELS

        The number of models, for iteration over the models in test code.

A curve given by a long Weierstrass equation over a ring in which 6 is
invertible can be put in short Weierstrass form by hand: with `c_4` and
`c_6` as computed by :func:`gr_ec_ctx_c_invariants`, the curve is isomorphic
to `y^2 = x^3 - 27 c_4 x - 54 c_6`. An interface for changes of variables
`(x, y) \mapsto (u^2 x + r, u^3 y + u^2 s x + t)`, which would also transport
points between the two models, is planned but not currently provided.

Point representation
-------------------------------------------------------------------------------

Points are represented by a triple `(X : Y : Z)` of elements of the base
ring, satisfying the homogenized curve equation

.. math::

    Y^2 Z + a_1 X Y Z + a_3 Y Z^2 = X^3 + a_2 X^2 Z + a_4 X Z^2 + a_6 Z^3.

The point at infinity `\mathcal{O}`, which is the identity element of the
group law, is `(0 : 1 : 0)`; a point with `Z` invertible corresponds to the
affine point `(X/Z, Y/Z)`. The group law is written additively, so the
identity element is created by :func:`gr_ec_point_zero` and tested for by
:func:`gr_ec_point_is_zero`.

Triples are not normalized automatically: the same point has many
representations, and :func:`gr_ec_point_equal` compares points by testing
whether the `2 \times 2` minors

.. math::

    X_1 Y_2 - X_2 Y_1, \quad Y_1 Z_2 - Y_2 Z_1, \quad X_1 Z_2 - X_2 Z_1

all vanish. This decides equality of points when the base ring is an
integral domain. Over a general commutative ring, the vanishing of the
minors is only a necessary condition for the two triples to define the same
point, and the result should be interpreted accordingly.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: gr_ec_ctx_struct

.. type:: gr_ec_ctx_t

    Contains a pointer to the base ring (``base_ring``), a pointer to an
    array of :macro:`GR_EC_CTX_NUM_COEFFS` elements of the base ring
    (``coeffs``) holding the invariants
    `a_1, a_2, a_3, a_4, a_6, b_2, b_4, b_6, b_8, \Delta` of the curve,
    and the model of the curve (``which_model``).

    A ``gr_ec_ctx_t`` is defined as an array of length one of type
    ``gr_ec_ctx_struct``, permitting a ``gr_ec_ctx_t`` to be passed by
    reference.

.. type:: gr_ec_point_struct

.. type:: gr_ec_point_t

    Contains a pointer (``coords``) to an array of three elements of the
    base ring, holding the projective coordinates `X`, `Y`, `Z` of the
    point.

    A ``gr_ec_point_t`` is defined as an array of length one of type
    ``gr_ec_point_struct``, permitting a ``gr_ec_point_t`` to be passed by
    reference.

.. macro:: GR_EC_CTX_NUM_COEFFS

    The number of base ring elements cached in a context object, currently
    10.

.. macro:: GR_EC_ELEM_CTX(ctx)
           GR_EC_SIZEOF_ELEM(ctx)

    The base ring of *ctx*, respectively the size in bytes of an element of
    that ring.

.. macro:: GR_EC_COEFF(ctx, i)

    Pointer to the cached coefficient of index *i* in *ctx*, where *i* is
    an integer with `0 \le i < 10`. The named macros below should normally
    be preferred.

.. macro:: GR_EC_A1(ctx)
           GR_EC_A2(ctx)
           GR_EC_A3(ctx)
           GR_EC_A4(ctx)
           GR_EC_A6(ctx)
           GR_EC_B2(ctx)
           GR_EC_B4(ctx)
           GR_EC_B6(ctx)
           GR_EC_B8(ctx)
           GR_EC_DISC(ctx)

    Pointers to the individual invariants `a_1, a_2, a_3, a_4, a_6`,
    `b_2, b_4, b_6, b_8` and `\Delta` cached in *ctx*.

.. macro:: GR_EC_POINT_X(P, ctx)
           GR_EC_POINT_Y(P, ctx)
           GR_EC_POINT_Z(P, ctx)

    Pointers to the projective coordinates of the point *P*.

Context object methods
-------------------------------------------------------------------------------

.. function:: int gr_ec_ctx_init(gr_ec_ctx_t ctx, gr_ctx_t base_ring, gr_srcptr a1, gr_srcptr a2, gr_srcptr a3, gr_srcptr a4, gr_srcptr a6)

    Initializes *ctx* to the elliptic curve

    .. math::

        y^2 + a_1 x y + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6

    over *base_ring*, where the coefficients are elements of *base_ring*.
    The `b`-invariants and the discriminant of the curve are computed and
    cached in the context object.

    The model of the curve is set to :macro:`GR_EC_SHORT_WEIERSTRASS` if
    `a_1`, `a_2` and `a_3` are all provably zero, and to
    :macro:`GR_EC_LONG_WEIERSTRASS` otherwise.

    Returns ``GR_DOMAIN`` if the discriminant of the curve is provably zero,
    that is, if the given equation defines a singular curve.
    Returns ``GR_UNABLE`` if the base ring is unable to compute the
    invariants. In both cases, *ctx* is left uninitialized and must not be
    cleared. If the discriminant can be computed but not proved to be
    nonzero, initialization succeeds; :func:`gr_ec_ctx_is_smooth` can then be
    used to inspect what is known about the curve.

    A *base_ring* which is not a commutative ring is not supported;
    ``GR_DOMAIN`` is returned in that case.

.. function:: int gr_ec_ctx_init_si(gr_ec_ctx_t ctx, gr_ctx_t base_ring, slong a1, slong a2, slong a3, slong a4, slong a6)

    As :func:`gr_ec_ctx_init`, with the `a`-invariants given as machine
    integers, mapped into *base_ring* by :func:`gr_set_si`.

.. function:: int gr_ec_ctx_init_short_weierstrass(gr_ec_ctx_t ctx, gr_ctx_t base_ring, gr_srcptr a4, gr_srcptr a6)
              int gr_ec_ctx_init_short_weierstrass_si(gr_ec_ctx_t ctx, gr_ctx_t base_ring, slong a4, slong a6)

    Initializes *ctx* to the curve `y^2 = x^3 + a_4 x + a_6` over
    *base_ring*, that is, to the curve with `a`-invariants
    `(0, 0, 0, a_4, a_6)`. The model is set to
    :macro:`GR_EC_SHORT_WEIERSTRASS`. The return value is as for
    :func:`gr_ec_ctx_init`; in particular, ``GR_DOMAIN`` is returned when
    `-16 (4 a_4^3 + 27 a_6^2)` is provably zero.

.. function:: int gr_ec_ctx_init_randtest(gr_ec_ctx_t ctx, flint_rand_t state, gr_ctx_t base_ring)

    Initializes *ctx* to a random nonsingular curve over *base_ring*,
    for use in test code. Returns ``GR_UNABLE`` if no such curve could be
    generated after a few attempts, which can happen for example over very
    small finite rings; in that case *ctx* is left uninitialized.

.. function:: void gr_ec_ctx_clear(gr_ec_ctx_t ctx)

    Clears the context object *ctx*, freeing the cached invariants.
    The base ring is not cleared.

.. function:: gr_ctx_struct * gr_ec_ctx_base_ring(gr_ec_ctx_t ctx)

    Returns a pointer to the base ring of *ctx*. The result is a borrowed
    reference which must not be cleared.

.. function:: gr_ec_model_t gr_ec_ctx_model(gr_ec_ctx_t ctx)

    Returns the model of the curve *ctx*.

.. function:: int gr_ec_ctx_write(gr_stream_t out, gr_ec_ctx_t ctx)
              int gr_ec_ctx_get_str(char ** res, gr_ec_ctx_t ctx)
              int gr_ec_ctx_print(gr_ec_ctx_t ctx)

    Writes a description of the curve *ctx* to the stream *out*, to a
    string, or to standard output. The description lists the `a`-invariants
    of the curve and the base ring.

Curve invariants
-------------------------------------------------------------------------------

.. function:: gr_srcptr gr_ec_ctx_a_invariants_srcptr(gr_ec_ctx_t ctx)
              gr_srcptr gr_ec_ctx_b_invariants_srcptr(gr_ec_ctx_t ctx)

    Returns a pointer to the array of the five `a`-invariants
    `a_1, a_2, a_3, a_4, a_6`, respectively of the four `b`-invariants
    `b_2, b_4, b_6, b_8`, cached in *ctx*. These are borrowed references
    valid until *ctx* is cleared, and they must not be modified.

.. function:: int gr_ec_ctx_a_invariants(gr_ptr a1, gr_ptr a2, gr_ptr a3, gr_ptr a4, gr_ptr a6, gr_ec_ctx_t ctx)

    Sets *a1*, *a2*, *a3*, *a4*, *a6* to copies of the `a`-invariants of
    the curve *ctx*.

.. function:: int gr_ec_ctx_b_invariants(gr_ptr b2, gr_ptr b4, gr_ptr b6, gr_ptr b8, gr_ec_ctx_t ctx)

    Sets *b2*, *b4*, *b6*, *b8* to the `b`-invariants of the curve *ctx*,
    given by

    .. math::

        b_2 &= a_1^2 + 4 a_2

        b_4 &= 2 a_4 + a_1 a_3

        b_6 &= a_3^2 + 4 a_6

        b_8 &= a_1^2 a_6 + 4 a_2 a_6 - a_1 a_3 a_4 + a_2 a_3^2 - a_4^2

.. function:: int gr_ec_ctx_c_invariants(gr_ptr c4, gr_ptr c6, gr_ec_ctx_t ctx)

    Sets *c4* and *c6* to the `c`-invariants of the curve *ctx*, given by

    .. math::

        c_4 &= b_2^2 - 24 b_4

        c_6 &= -b_2^3 + 36 b_2 b_4 - 216 b_6

.. function:: int gr_ec_ctx_discriminant(gr_ptr res, gr_ec_ctx_t ctx)

    Sets *res* to the discriminant

    .. math::

        \Delta = -b_2^2 b_8 - 8 b_4^3 - 27 b_6^2 + 9 b_2 b_4 b_6

    of the curve *ctx*. The discriminant is cached in the context object, so
    this only copies the stored value.

.. function:: int gr_ec_ctx_j_invariant(gr_ptr res, gr_ec_ctx_t ctx)

    Sets *res* to the `j`-invariant `j = c_4^3 / \Delta` of the curve *ctx*.
    Returns ``GR_DOMAIN`` if the discriminant is not invertible in the base
    ring, and ``GR_UNABLE`` if invertibility cannot be decided or the
    division cannot be performed.

.. function:: truth_t gr_ec_ctx_is_smooth(gr_ec_ctx_t ctx)

    Returns ``T_TRUE`` if the discriminant of *ctx* is provably nonzero, and
    ``T_UNKNOWN`` if this cannot be decided in the base ring. The value
    ``T_FALSE`` is never returned, since a context object cannot be
    initialized with a provably zero discriminant.

    Note that over a base ring which is not a field, a nonzero discriminant
    does not imply that the equation defines a smooth curve over the whole
    base: for that, `\Delta` must be invertible.

Memory management
-------------------------------------------------------------------------------

.. function:: void gr_ec_point_init(gr_ec_point_t P, gr_ec_ctx_t ctx)

    Initializes *P* for use as a point of the curve *ctx*, and sets it to
    the point at infinity `(0 : 1 : 0)`.

.. function:: void gr_ec_point_clear(gr_ec_point_t P, gr_ec_ctx_t ctx)

    Clears the point *P*.

.. function:: void gr_ec_point_swap(gr_ec_point_t P, gr_ec_point_t Q, gr_ec_ctx_t ctx)

    Swaps *P* and *Q* efficiently.

.. function:: gr_ptr gr_ec_point_x_ptr(gr_ec_point_t P, gr_ec_ctx_t ctx)
              gr_ptr gr_ec_point_y_ptr(gr_ec_point_t P, gr_ec_ctx_t ctx)
              gr_ptr gr_ec_point_z_ptr(gr_ec_point_t P, gr_ec_ctx_t ctx)
              gr_srcptr gr_ec_point_x_srcptr(const gr_ec_point_t P, gr_ec_ctx_t ctx)
              gr_srcptr gr_ec_point_y_srcptr(const gr_ec_point_t P, gr_ec_ctx_t ctx)
              gr_srcptr gr_ec_point_z_srcptr(const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Returns a pointer to the projective coordinate `X`, `Y` or `Z` of *P*.
    Writing to a coordinate through the non-const version can leave *P* off
    the curve; it is the responsibility of the caller to restore a valid
    representation.

Basic manipulation
-------------------------------------------------------------------------------

.. function:: int gr_ec_point_set(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to a copy of *P*.

.. function:: int gr_ec_point_zero(gr_ec_point_t res, gr_ec_ctx_t ctx)

    Sets *res* to the point at infinity `\mathcal{O} = (0 : 1 : 0)`, the
    identity element of the group law.

.. function:: int gr_ec_point_set_affine(gr_ec_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx)
              int _gr_ec_point_set_affine(gr_ec_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx)

    Sets *res* to the point with affine coordinates `(x, y)`, that is, to
    `(x : y : 1)`.

    The non-underscore version verifies that the point lies on the curve,
    returning ``GR_DOMAIN`` if it provably does not and ``GR_UNABLE`` if
    this cannot be decided. The underscore version performs no check and
    always succeeds unless the base ring fails to copy the coordinates.

.. function:: int gr_ec_point_set_projective(gr_ec_point_t res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ec_ctx_t ctx)
              int _gr_ec_point_set_projective(gr_ec_point_t res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ec_ctx_t ctx)

    Sets *res* to the point with projective coordinates `(x : y : z)`.
    Checking is as for :func:`gr_ec_point_set_affine`; in addition, the
    non-underscore version returns ``GR_DOMAIN`` if `x`, `y` and `z` are all
    zero, which is not a valid representation of a point.

.. function:: int gr_ec_point_get_affine(gr_ptr x, gr_ptr y, const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Sets *x* and *y* to the affine coordinates `X/Z` and `Y/Z` of *P*.
    Returns ``GR_DOMAIN`` if *P* is the point at infinity or, more
    generally, if `Z` is not invertible in the base ring, and ``GR_UNABLE``
    if this cannot be decided or the division cannot be performed.

.. function:: int gr_ec_point_normalize(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to the normalized representation of *P*: this is
    `(X/Z : Y/Z : 1)` if `Z` is invertible, and `(0 : 1 : 0)` if *P* is the
    point at infinity. Returns ``GR_DOMAIN`` if `Z` is neither zero nor
    invertible, which can happen over a base ring that is not a field, and
    ``GR_UNABLE`` if the base ring cannot decide whether `Z` is zero or
    invertible, or cannot perform the division. In both cases *res* is left
    in an unspecified state.

    Normalization is never required for correctness, but it makes
    subsequent operations on the same point cheaper and gives a canonical
    representation for printing.

.. function:: int gr_ec_point_lift_x(gr_ec_point_t res, gr_srcptr x, gr_ec_ctx_t ctx)

    Sets *res* to a point of the curve with affine `x`-coordinate *x*.
    When there are two such points, which of the two is returned is
    unspecified but deterministic; the other one is its negative, obtained
    with :func:`gr_ec_point_neg`.

    Returns ``GR_DOMAIN`` if no such point exists over the base ring, and
    ``GR_UNABLE`` if this cannot be decided, which includes the case of a
    curve in the model :macro:`GR_EC_LONG_WEIERSTRASS` over a base ring in
    which 2 is not invertible.

Comparisons and properties
-------------------------------------------------------------------------------

.. function:: truth_t gr_ec_point_is_zero(const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Returns whether *P* is the point at infinity, that is, whether its
    coordinate `Z` is zero.

.. function:: truth_t gr_ec_point_equal(const gr_ec_point_t P, const gr_ec_point_t Q, gr_ec_ctx_t ctx)

    Returns whether *P* and *Q* are the same point, by testing the
    vanishing of the `2 \times 2` minors of the matrix formed by their
    coordinates, as described above. Returns ``T_UNKNOWN`` if the base ring
    cannot decide whether the minors are zero.

.. function:: truth_t gr_ec_point_is_on_curve(const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Returns whether the coordinates of *P* satisfy the homogenized
    Weierstrass equation of *ctx*. This is checked in the general
    (long) form, and is therefore supported for every model.

Input and output
-------------------------------------------------------------------------------

.. function:: int gr_ec_point_write(gr_stream_t out, const gr_ec_point_t P, gr_ec_ctx_t ctx)
              int gr_ec_point_get_str(char ** res, const gr_ec_point_t P, gr_ec_ctx_t ctx)
              int gr_ec_point_print(const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Writes *P* to the stream *out*, to a string, or to standard output.
    The point is printed as the triple of its projective coordinates
    ``(X : Y : Z)``, without normalizing it first.

Random generation
-------------------------------------------------------------------------------

.. function:: int gr_ec_point_randtest(gr_ec_point_t res, flint_rand_t state, gr_ec_ctx_t ctx)

    Sets *res* to a random point of the curve *ctx*, for use in test code.
    The point at infinity is generated with some positive probability, and
    the projective representation of the result is not necessarily
    normalized.

    Returns ``GR_UNABLE`` if no point could be generated, which is the case
    whenever :func:`gr_ec_point_lift_x` is unable to lift the `x`-coordinates
    that were tried.

Arithmetic
-------------------------------------------------------------------------------

The functions in this section implement the group law of the curve, written
additively with the point at infinity as the identity element. They
dispatch on the model of *ctx*; see :ref:`the section on curve models
<gr-ec-models>` for which models are currently implemented.

.. function:: int gr_ec_point_neg(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to `-P`, which in projective coordinates is
    `(X : -Y - a_1 X - a_3 Z : Z)`. This is supported for every model.

.. function:: int gr_ec_point_add(gr_ec_point_t res, const gr_ec_point_t P, const gr_ec_point_t Q, gr_ec_ctx_t ctx)
              int gr_ec_point_sub(gr_ec_point_t res, const gr_ec_point_t P, const gr_ec_point_t Q, gr_ec_ctx_t ctx)

    Sets *res* to `P + Q`, respectively `P - Q`.

    The addition formulas are not uniform: the cases where one of the
    points is the point at infinity, and the case `P = \pm Q`, are detected
    and handled separately. Consequently, these functions return
    ``GR_UNABLE`` when the base ring cannot decide the corresponding
    equalities, even when the generic formulas would apply.

.. function:: int gr_ec_point_dbl(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to `2 P`. This is more efficient than calling
    :func:`gr_ec_point_add` with two aliased arguments.

.. function:: int gr_ec_point_mul_ui(gr_ec_point_t res, const gr_ec_point_t P, ulong n, gr_ec_ctx_t ctx)
              int gr_ec_point_mul_si(gr_ec_point_t res, const gr_ec_point_t P, slong n, gr_ec_ctx_t ctx)
              int gr_ec_point_mul_fmpz(gr_ec_point_t res, const gr_ec_point_t P, const fmpz_t n, gr_ec_ctx_t ctx)

    Sets *res* to `n P`, the sum of *n* copies of *P*, with the convention
    that `n P = (-n)(-P)` for negative *n* and that `0 P` is the point at
    infinity.

Model-specific implementations
-------------------------------------------------------------------------------

The following low-level functions implement the group law for a specific
model. They assume that the curve *ctx* is given in that model, in
particular that `a_1 = a_2 = a_3 = 0` for the short Weierstrass versions,
and they do not check this; the results are undefined otherwise. Use the
functions of the previous section unless the model is known in advance.

.. function:: int _gr_ec_point_neg_long_weierstrass(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx)
              int _gr_ec_point_neg_short_weierstrass(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Implementations of :func:`gr_ec_point_neg`.

.. function:: int _gr_ec_point_add_long_weierstrass(gr_ec_point_t res, const gr_ec_point_t P, const gr_ec_point_t Q, gr_ec_ctx_t ctx)
              int _gr_ec_point_add_short_weierstrass(gr_ec_point_t res, const gr_ec_point_t P, const gr_ec_point_t Q, gr_ec_ctx_t ctx)

    Implementations of :func:`gr_ec_point_add`, using the standard
    homogeneous projective chord-and-tangent formulas for the respective
    model. The long Weierstrass version is not yet implemented and returns
    ``GR_UNABLE``.

.. function:: int _gr_ec_point_dbl_long_weierstrass(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx)
              int _gr_ec_point_dbl_short_weierstrass(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Implementations of :func:`gr_ec_point_dbl`. The long Weierstrass
    version is not yet implemented and returns ``GR_UNABLE``.

.. function:: int _gr_ec_point_mul_fmpz_binary(gr_ec_point_t res, const gr_ec_point_t P, const fmpz_t n, gr_ec_ctx_t ctx)

    Implementation of :func:`gr_ec_point_mul_fmpz` using left-to-right
    binary double-and-add, calling :func:`gr_ec_point_dbl` and
    :func:`gr_ec_point_add` and therefore working for any model for which
    those are implemented. Aliasing of *res* and *P* is allowed.
