.. _gr-ec:

**gr_ec.h** -- elliptic curves over generic rings
===============================================================================

.. note::

    This module is under construction.

A :type:`gr_ec_ctx_t` represents an elliptic curve `E` over a generic
commutative ring *R*, given (by default) by a general (long) Weierstrass equation

.. math::

    E : y^2 + a_1 x y + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6

with coefficients `a_1, a_2, a_3, a_4, a_6 \in R`.

Points of `E` are represented by one of three types --
:type:`gr_ec_point_t`, :type:`gr_ec_aff_point_t` and
:type:`gr_ec_jac_point_t` -- which differ in the coordinate system; see
:ref:`the section on point representations <gr-ec-representations>`.
Each type has its own set of arithmetic functions, and conversions between
the three are provided.

.. _gr-ec-models:

Curve models
-------------------------------------------------------------------------------

The point arithmetic is dispatched on the *model* of the curve; multiple models
are available all with their pros and cons in performance/versatility.

.. type:: gr_ec_model_t

    Represents one of the following supported curve models:

    .. macro:: GR_EC_LONG_WEIERSTRASS

        A curve `y^2 + a_1 x y + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6`
        with at least one of `a_1, a_2, a_3` none zero.
        This model is valid over any base ring, and is the one to use in
        residue characteristic 2 and 3, where the short Weierstrass
        equation is always singular.


    .. macro:: GR_EC_SHORT_WEIERSTRASS

        A curve `y^2 = x^3 + a_4 x + a_6`, that is, a curve with
        `a_1 = a_2 = a_3 = 0`. Every elliptic curve over a ring in which 6
        is invertible -- in particular over a field of characteristic 0 or
        of characteristic `p > 3` -- is isomorphic to a curve in this model.

.. _gr-ec-representations:

Point representations
-------------------------------------------------------------------------------

Three representations of points are available. They all represent the
same point and the choice between them is a matter of efficiency:

.. list-table::
   :header-rows: 1
   :widths: 22 22 20 36

   * - Type
     - Stored data
     - Base ring
     - Group law
   * - :type:`gr_ec_point_t`
     - `X, Y, Z`
     - any commutative ring
     - inversion-free; the interchange representation
   * - :type:`gr_ec_aff_point_t`
     - `x, y` and a flag
     - field only
     - one inversion per operation
   * - :type:`gr_ec_jac_point_t`
     - `X, Y, Z` and a flag
     - any commutative ring
     - inversion-free; fastest for scalar multiplication

:type:`gr_ec_point_t` uses homogeneous projective coordinates
`(X : Y : Z)`, satisfying the homogenized curve equation

.. math::

    Y^2 Z + a_1 X Y Z + a_3 Y Z^2 = X^3 + a_2 X^2 Z + a_4 X Z^2 + a_6 Z^3.

The point at infinity `\mathcal{O}` is `(0 : 1 : 0)` by convention, and a point with `Z`
invertible corresponds to the affine point `(X/Z, Y/Z)`. This is the default
representation.

:type:`gr_ec_aff_point_t` stores the two affine coordinates `x` and `y`
satisfying the Weierstrass equation, together with a :type:`truth_t` flag
recording whether the point is `\mathcal{O}`. This representation requires 
the base ring to be a field; see the section on affine points below. This
representation is the most compact and the most convenient for input,
output and testing, but it is the slowest for repeated arithmetic.

:type:`gr_ec_jac_point_t` uses Jacobian coordinates `(X, Y, Z)`, in which
the affine point represented is `(X/Z^2, Y/Z^3)`, satisfying

.. math::

    Y^2 + a_1 X Y Z + a_3 Y Z^3 = X^3 + a_2 X^2 Z^2 + a_4 X Z^4 + a_6 Z^6.

As for affine points, the point at infinity is recorded by a separate
:type:`truth_t` flag rather than by a zero coordinate. Storing the flag
means that `Z` is nonzero in every valid finite point, so the doubling and
addition formulas do not have to test `Z` against zero; the cost is one
machine word per point. Jacobian coordinates give the cheapest doubling of
the three representations, which is what makes them the representation of
choice for scalar multiplication.

The infinity flag
...............................................................................

For :type:`gr_ec_aff_point_t` and :type:`gr_ec_jac_point_t`, the field
``is_infinity`` has the following meaning:

* ``T_TRUE`` -- the point is `\mathcal{O}`. The stored coordinates are
  undefined behavior.
* ``T_FALSE`` -- the point is the finite point given by the stored
  coordinates.
* ``T_UNKNOWN`` -- it is not known whether the point is `\mathcal{O}`.
  Such a point is *invalid*: the stored coordinates are unspecified,
  predicates applied to it return ``T_UNKNOWN``, and functions taking it
  as an input return ``GR_UNABLE``.


Equality of points
...............................................................................

Points are not normalized automatically, so the same point of `E` has many
representations for `gr_ec_point_t`:

:func:`gr_ec_point_equal` tests equality by cross-multiplying, so equality is 
guaranteed only when the base ring is an integral domain. Over a general 
commutative ring they are only necessary conditions for the two representations 
to define the same point, and the result should be interpreted accordingly. 

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: gr_ec_ctx_struct

.. type:: gr_ec_ctx_t

    Contains a pointer to the base ring (``base_ring``), a pointer to an
    array of :macro:`GR_EC_CTX_NUM_COEFFS` elements of the base ring
    (``coeffs``) holding the invariants
    `a_1, a_2, a_3, a_4, a_6, b_2, b_4, b_6, b_8, \Delta` of the curve,
    and the model of the curve (``model``).

.. type:: gr_ec_point_struct

.. type:: gr_ec_point_t

    Contains a pointer (``coords``) to an array of three elements of the
    base ring, holding the homogeneous projective coordinates `X`, `Y`, `Z`
    of the point.

.. type:: gr_ec_aff_point_struct

.. type:: gr_ec_aff_point_t

    Contains a pointer (``coords``) to an array of two elements of the base
    ring, holding the affine coordinates `x`, `y` of the point, and a
    :type:`truth_t` (``is_infinity``) recording whether the point is the
    point at infinity.

.. type:: gr_ec_jac_point_struct

.. type:: gr_ec_jac_point_t

    Contains a pointer (``coords``) to an array of three elements of the
    base ring, holding the Jacobian coordinates `X`, `Y`, `Z` of the point,
    and a :type:`truth_t` (``is_infinity``) recording whether the point is
    the point at infinity.

    Each of the three point types is defined as an array of length one of
    the corresponding struct type, permitting it to be passed by reference.

.. macro:: GR_EC_CTX_NUM_COEFFS

    The maximum number of base ring elements cached in a context object, currently
    10.

.. macro:: GR_EC_ELEM_CTX(ctx)
           GR_EC_SIZEOF_ELEM(ctx)

    The base ring of *ctx*, respectively the size in bytes of an element of
    that ring.

.. macro:: GR_EC_COEFF(ctx, i)

    Pointer to the cached coefficient of index *i* in *ctx*, where *i* is
    an integer with `0 \le i < GR_EC_CTX_NUM_COEFFS`.

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
           GR_EC_AFF_POINT_X(P, ctx)
           GR_EC_AFF_POINT_Y(P, ctx)
           GR_EC_JAC_POINT_X(P, ctx)
           GR_EC_JAC_POINT_Y(P, ctx)
           GR_EC_JAC_POINT_Z(P, ctx)

    Pointers to the stored coordinates of the point *P*, for each of the
    three point types.

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

    Returns a pointer to the base ring of *ctx*.

.. function:: gr_ec_model_t gr_ec_ctx_model(gr_ec_ctx_t ctx)

    Returns the model of the curve *ctx*.

.. function:: truth_t gr_ec_ctx_is_over_field(gr_ec_ctx_t ctx)

    Returns whether the base ring of *ctx* is a field. This is the
    condition under which the affine point functions are usable.

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

Projective points: memory management
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

Projective points: basic manipulation
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

    `gr_ec_point_set_affine` verifies that the point lies on the curve,
    returning ``GR_DOMAIN`` if it provably does not and ``GR_UNABLE`` if
    this cannot be decided. `_gr_ec_point_set_affine` performs no check and
    always succeeds unless the base ring fails to copy the coordinates.

.. function:: int gr_ec_point_set_projective(gr_ec_point_t res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ec_ctx_t ctx)
              int _gr_ec_point_set_projective(gr_ec_point_t res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ec_ctx_t ctx)

    Sets *res* to the point with projective coordinates `(x : y : z)`.
    Checking is as for :func:`gr_ec_point_set_affine`; in addition,
    `gr_ec_point_set_projective` returns ``GR_DOMAIN`` if `x`, `y` and `z` are all
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

.. function:: int gr_ec_point_lift_x(gr_ec_point_t res, gr_srcptr x, gr_ec_ctx_t ctx)

    Sets *res* to a point of the curve with affine `x`-coordinate *x*.
    When there are two such points, which of the two is returned is
    unspecified;

    Returns ``GR_DOMAIN`` if no such point exists over the base ring, and
    ``GR_UNABLE`` if this cannot be decided, which includes the case of a
    curve in the model :macro:`GR_EC_LONG_WEIERSTRASS` over a base ring in
    which 2 is not invertible.

Projective points: comparisons and properties
-------------------------------------------------------------------------------

.. function:: truth_t gr_ec_point_is_inf(const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Returns whether *P* is the point at infinity.

.. function:: truth_t gr_ec_point_equal(const gr_ec_point_t P, const gr_ec_point_t Q, gr_ec_ctx_t ctx)

    Returns whether *P* and *Q* are the same point. Returns ``T_UNKNOWN`` if the base ring
    cannot allow to decide equality.

.. function:: truth_t gr_ec_point_is_on_curve(const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Returns whether the coordinates of *P* satisfy the homogenized
    Weierstrass equation of *ctx*. This is checked in the general
    (long) form, and is therefore supported for every model.

Projective points: input, output and random generation
-------------------------------------------------------------------------------

.. function:: int gr_ec_point_write(gr_stream_t out, const gr_ec_point_t P, gr_ec_ctx_t ctx)
              int gr_ec_point_get_str(char ** res, const gr_ec_point_t P, gr_ec_ctx_t ctx)
              int gr_ec_point_print(const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Writes *P* to the stream *out*, to a string, or to standard output.
    The point is printed as the triple of its projective coordinates
    ``(X : Y : Z)``.

.. function:: int gr_ec_point_randtest(gr_ec_point_t res, flint_rand_t state, gr_ec_ctx_t ctx)

    Sets *res* to a random point of the curve *ctx*, for use in test code.
    The point at infinity is generated with some positive probability, and
    the projective representation of the result is not necessarily
    normalized.

    Returns ``GR_UNABLE`` if no point could be generated, which is the case
    whenever :func:`gr_ec_point_lift_x` is unable to lift the `x`-coordinates
    that were tried.

Projective points: arithmetic
-------------------------------------------------------------------------------

The functions in this section implement the group law of the curve.

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

    The work is done in Jacobian coordinates, where the ladder needs no
    inversion and the conversions in and out are inversion-free.

Affine points
-------------------------------------------------------------------------------

The affine group law divides, so it is only defined over a field:
:func:`gr_ec_aff_point_add`, :func:`gr_ec_aff_point_sub`,
:func:`gr_ec_aff_point_dbl` and the scalar multiplication functions return
``GR_DOMAIN`` when :func:`gr_ec_ctx_is_over_field` returns ``T_FALSE``.
A ring that merely does not advertise itself as a field is still accepted,
and the division then reports ``GR_DOMAIN`` itself if it meets a non-unit.

The rest of this section -- memory management, the predicates, and the
functions that only move coordinates around -- works over any base ring;
so do the conversions, which fail through their inversion exactly when the
affine point does not exist.

.. function:: void gr_ec_aff_point_init(gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
              void gr_ec_aff_point_clear(gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
              void gr_ec_aff_point_swap(gr_ec_aff_point_t P, gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)

    Initializes *P* to the point at infinity, clears *P*, respectively
    swaps *P* and *Q* efficiently.

.. function:: gr_ptr gr_ec_aff_point_x_ptr(gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
              gr_ptr gr_ec_aff_point_y_ptr(gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
              gr_srcptr gr_ec_aff_point_x_srcptr(const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
              gr_srcptr gr_ec_aff_point_y_srcptr(const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Returns a pointer to the affine coordinate `x` or `y` of *P*.

.. function:: int gr_ec_aff_point_set(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to a copy of *P*.

.. function:: int gr_ec_aff_point_zero(gr_ec_aff_point_t res, gr_ec_ctx_t ctx)

    Sets *res* to the point at infinity, that is, sets ``is_infinity`` to
    ``T_TRUE``.

.. function:: int gr_ec_aff_point_set_affine(gr_ec_aff_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx)
              int _gr_ec_aff_point_set_affine(gr_ec_aff_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx)

    Sets *res* to the finite point with affine coordinates `(x, y)`.
    `gr_ec_aff_point_set_affine` verifies that the point lies on the curve; 
    `_gr_ec_aff_point_set_affine` performs no check.

.. function:: int gr_ec_aff_point_get_affine(gr_ptr x, gr_ptr y, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Sets *x* and *y* to the affine coordinates of *P*, which for this
    representation is just a copy. Returns ``GR_DOMAIN`` if *P* is the
    point at infinity.

.. function:: int gr_ec_aff_point_lift_x(gr_ec_aff_point_t res, gr_srcptr x, gr_ec_ctx_t ctx)

    Sets *res* to a point of the curve with `x`-coordinate *x*, as
    :func:`gr_ec_point_lift_x`.

.. function:: truth_t gr_ec_aff_point_is_inf(const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Returns ``P->is_infinity``, that is, whether *P* is the point at
    infinity.

.. function:: truth_t gr_ec_aff_point_equal(const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)

    Returns whether *P* and *Q* are the same point. Since the affine
    representation is unique, this compares the infinity flags and then the
    coordinates.

.. function:: truth_t gr_ec_aff_point_is_on_curve(const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Returns whether *P* is the point at infinity or its coordinates satisfy
    the Weierstrass equation of *ctx*, in the general (long) form.

.. function:: int gr_ec_aff_point_write(gr_stream_t out, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
              int gr_ec_aff_point_get_str(char ** res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
              int gr_ec_aff_point_print(const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Writes *P* to the stream *out*, to a string, or to standard output.
    A finite point is printed as ``(x, y)`` and the point at infinity as
    ``O``.

.. function:: int gr_ec_aff_point_randtest(gr_ec_aff_point_t res, flint_rand_t state, gr_ec_ctx_t ctx)

    Sets *res* to a random point of the curve *ctx*, for use in test code,
    as :func:`gr_ec_point_randtest`.

.. function:: int gr_ec_aff_point_neg(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to `-P`, which in affine coordinates is
    `(x, -y - a_1 x - a_3)`. This is supported for every model and requires
    no inversion.

.. function:: int gr_ec_aff_point_add(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)
              int gr_ec_aff_point_sub(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)
              int gr_ec_aff_point_dbl(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to `P + Q`, `P - Q`, respectively `2 P`, using the
    chord-and-tangent formulas. Each of these performs one inversion in the
    base ring.

    The exceptional cases are decided from the infinity flags and from the
    equality of the `x`-coordinates, so no comparison of full points is
    needed.

.. function:: int gr_ec_aff_point_mul_ui(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, ulong n, gr_ec_ctx_t ctx)
              int gr_ec_aff_point_mul_si(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, slong n, gr_ec_ctx_t ctx)
              int gr_ec_aff_point_mul_fmpz(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const fmpz_t n, gr_ec_ctx_t ctx)

    Sets *res* to `n P`.

    These convert *P* to Jacobian coordinates, run the scalar
    multiplication there, and convert the result back, so that a single
    inversion is performed instead of one per bit of *n*.

Jacobian points
-------------------------------------------------------------------------------

.. function:: void gr_ec_jac_point_init(gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
              void gr_ec_jac_point_clear(gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
              void gr_ec_jac_point_swap(gr_ec_jac_point_t P, gr_ec_jac_point_t Q, gr_ec_ctx_t ctx)

    Initializes *P* to the point at infinity, clears *P*, respectively
    swaps *P* and *Q* efficiently.

.. function:: gr_ptr gr_ec_jac_point_x_ptr(gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
              gr_ptr gr_ec_jac_point_y_ptr(gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
              gr_ptr gr_ec_jac_point_z_ptr(gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
              gr_srcptr gr_ec_jac_point_x_srcptr(const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
              gr_srcptr gr_ec_jac_point_y_srcptr(const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
              gr_srcptr gr_ec_jac_point_z_srcptr(const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Returns a pointer to the Jacobian coordinate `X`, `Y` or `Z` of *P*.
    These are the weighted coordinates, not the affine ones: the point
    represented is `(X/Z^2, Y/Z^3)`. They are only meaningful when
    ``P->is_infinity`` is ``T_FALSE``.

.. function:: int gr_ec_jac_point_set(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to a copy of *P*.

.. function:: int gr_ec_jac_point_zero(gr_ec_jac_point_t res, gr_ec_ctx_t ctx)

    Sets *res* to the point at infinity, that is, sets ``is_infinity`` to
    ``T_TRUE``.

.. function:: int gr_ec_jac_point_set_affine(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx)
              int _gr_ec_jac_point_set_affine(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx)

    Sets *res* to the finite point with affine coordinates `(x, y)`, that
    is, to the Jacobian triple `(x, y, 1)`. Checking is as for
    :func:`gr_ec_point_set_affine`.

.. function:: int gr_ec_jac_point_set_jacobian(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ec_ctx_t ctx)
              int _gr_ec_jac_point_set_jacobian(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ec_ctx_t ctx)

    Sets *res* to the finite point with Jacobian coordinates `(x, y, z)`,
    representing the affine point `(x/z^2, y/z^3)`.

    The non-underscore version verifies that the point lies on the curve
    and returns ``GR_DOMAIN`` if *z* is provably zero; the point at
    infinity must be created with :func:`gr_ec_jac_point_zero` instead. The
    underscore version performs no check.

.. function:: int gr_ec_jac_point_get_affine(gr_ptr x, gr_ptr y, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Sets *x* and *y* to the affine coordinates `X/Z^2` and `Y/Z^3` of *P*.
    Returns ``GR_DOMAIN`` if *P* is the point at infinity or if `Z` is not
    invertible in the base ring, and ``GR_UNABLE`` if this cannot be
    decided or the division cannot be performed. One inversion is
    performed.

.. function:: int gr_ec_jac_point_normalize(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to the normalized representation `(X/Z^2, Y/Z^3, 1)` of *P*,
    or to the point at infinity if *P* is. Failure is as for
    :func:`gr_ec_jac_point_get_affine`, and on failure *res* is set to an
    invalid point.

    Normalizing before a run of mixed additions is worthwhile, since a
    normalized point can be converted to an affine point for free.

.. function:: truth_t gr_ec_jac_point_is_inf(const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Returns ``P->is_infinity``, that is, whether *P* is the point at
    infinity.

.. function:: truth_t gr_ec_jac_point_equal(const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx)

    Returns whether *P* and *Q* are the same point, by testing
    `X_1 Z_2^2 = X_2 Z_1^2` and `Y_1 Z_2^3 = Y_2 Z_1^3` for finite points.

.. function:: truth_t gr_ec_jac_point_is_on_curve(const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Returns whether *P* is the point at infinity or its coordinates satisfy
    the weighted-homogeneous curve equation

    .. math::

        Y^2 + a_1 X Y Z + a_3 Y Z^3 = X^3 + a_2 X^2 Z^2 + a_4 X Z^4 + a_6 Z^6.

.. function:: int gr_ec_jac_point_write(gr_stream_t out, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_get_str(char ** res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_print(const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Writes *P* to the stream *out*, to a string, or to standard output.
    A finite point is printed as its Jacobian triple ``(X, Y, Z)``, without
    normalizing it first, and the point at infinity as ``O``.

.. function:: int gr_ec_jac_point_randtest(gr_ec_jac_point_t res, flint_rand_t state, gr_ec_ctx_t ctx)

    Sets *res* to a random point of the curve *ctx*, for use in test code,
    as :func:`gr_ec_point_randtest`. The `Z` coordinate of the result is
    not necessarily 1.

.. function:: int gr_ec_jac_point_neg(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to `-P`, which in Jacobian coordinates is
    `(X, -Y - a_1 X Z - a_3 Z^3, Z)`. This is supported for every model.

.. function:: int gr_ec_jac_point_add(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_sub(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx)

    Sets *res* to `P + Q`, respectively `P - Q`, using the standard
    Jacobian addition formulas. As in projective coordinates, the formulas
    are not uniform: the case `P = \pm Q` is detected from the intermediate
    quantities and handled separately, and ``GR_UNABLE`` is returned when
    the base ring cannot decide the corresponding equalities.

.. function:: int gr_ec_jac_point_dbl(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to `2 P`. This is the cheapest doubling of the three
    representations and is what makes Jacobian coordinates the default for
    scalar multiplication.

.. function:: int gr_ec_jac_point_add_aff_point(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_sub_aff_point(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)

    Sets *res* to `P + Q`, respectively `P - Q`, where *Q* is an affine
    point. These are the mixed-addition formulas: knowing that the second
    operand has `Z = 1` saves several multiplications compared to
    :func:`gr_ec_jac_point_add`, which is the reason to keep precomputed
    tables of affine points.

    Unlike the other affine functions, these do not require the base ring
    to be a field, since no inversion is performed.

.. function:: int gr_ec_jac_point_mul_ui(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, ulong n, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_mul_si(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, slong n, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_mul_fmpz(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpz_t n, gr_ec_ctx_t ctx)

    Sets *res* to `n P`, with the same conventions as
    :func:`gr_ec_point_mul_fmpz`, using
    :func:`_gr_ec_jac_point_mul_fmpz_naf` and falling back to
    :func:`_gr_ec_jac_point_mul_fmpz_binary` over a base ring where the
    window table cannot be normalized.

Conversions between representations
-------------------------------------------------------------------------------

The conversions to :type:`gr_ec_point_t` and :type:`gr_ec_jac_point_t` are
inversion-free and work over any base ring; the conversions to
:type:`gr_ec_aff_point_t` perform one inversion and therefore require the
base ring to be a field:

.. list-table::
   :header-rows: 1
   :widths: 30 46 24

   * - Conversion
     - Function
     - Cost
   * - affine to projective
     - :func:`gr_ec_point_set_aff_point`
     - copy
   * - affine to Jacobian
     - :func:`gr_ec_jac_point_set_aff_point`
     - copy
   * - Jacobian to projective
     - :func:`gr_ec_point_set_jac_point`
     - 2M + 1S
   * - projective to Jacobian
     - :func:`gr_ec_jac_point_set_point`
     - 2M + 1S
   * - projective to affine
     - :func:`gr_ec_aff_point_set_point`
     - 1I + 2M
   * - Jacobian to affine
     - :func:`gr_ec_aff_point_set_jac_point`
     - 1I + 3M + 1S

Here M, S and I denote a multiplication, a squaring and an inversion in the
base ring.

.. function:: int gr_ec_point_set_aff_point(gr_ec_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to the projective point `(x : y : 1)` corresponding to the
    affine point *P*, or to `(0 : 1 : 0)` if *P* is the point at infinity.

.. function:: int gr_ec_point_set_jac_point(gr_ec_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to the projective point `(X Z : Y : Z^3)` corresponding to
    the Jacobian point *P*, or to `(0 : 1 : 0)` if *P* is the point at
    infinity.

.. function:: int gr_ec_jac_point_set_point(gr_ec_jac_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to the Jacobian point `(X Z, Y Z^2, Z)` corresponding to the
    projective point *P*, or to the point at infinity if `Z` is zero.

    Returns ``GR_UNABLE``, and sets *res* to an invalid point, if the base
    ring cannot decide whether `Z` is zero.

.. function:: int gr_ec_jac_point_set_aff_point(gr_ec_jac_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to the Jacobian point `(x, y, 1)` corresponding to the
    affine point *P*, or to the point at infinity if *P* is. This is a copy
    of the coordinates and the flag, with no arithmetic.

.. function:: int gr_ec_aff_point_set_point(gr_ec_aff_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to the affine point `(X/Z, Y/Z)` corresponding to the
    projective point *P*, or to the point at infinity if `Z` is zero.

.. function:: int gr_ec_aff_point_set_jac_point(gr_ec_aff_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to the affine point `(X/Z^2, Y/Z^3)` corresponding to the
    Jacobian point *P*, or to the point at infinity if *P* is.

Model-specific implementations
-------------------------------------------------------------------------------

The following low-level functions implement the group law of one
representation for one specific model. They assume that the curve *ctx* is
given in that model, in particular that `a_1 = a_2 = a_3 = 0` for the short
Weierstrass versions, and they do not check this; the results are undefined
otherwise. Use the functions of the sections above unless the model is
known in advance.

Negation is not split by model, since the general formula specializes
correctly when `a_1 = a_3 = 0`.

.. function:: int _gr_ec_point_add_long_weierstrass(gr_ec_point_t res, const gr_ec_point_t P, const gr_ec_point_t Q, gr_ec_ctx_t ctx)
              int _gr_ec_point_add_short_weierstrass(gr_ec_point_t res, const gr_ec_point_t P, const gr_ec_point_t Q, gr_ec_ctx_t ctx)
              int _gr_ec_point_dbl_long_weierstrass(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx)
              int _gr_ec_point_dbl_short_weierstrass(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx)

    Implementations of :func:`gr_ec_point_add` and :func:`gr_ec_point_dbl`
    in homogeneous projective coordinates. The short Weierstrass versions
    are add-1998-cmo-2 and dbl-1998-cmo-2; the long Weierstrass ones are
    obtained by homogenizing the affine chord-and-tangent formulas, and
    reduce to the former when `a_1 = a_2 = a_3 = 0`.

.. function:: int _gr_ec_aff_point_add_long_weierstrass(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)
              int _gr_ec_aff_point_add_short_weierstrass(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)
              int _gr_ec_aff_point_dbl_long_weierstrass(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
              int _gr_ec_aff_point_dbl_short_weierstrass(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Implementations of :func:`gr_ec_aff_point_add` and
    :func:`gr_ec_aff_point_dbl`, computing the chord, respectively tangent,
    slope and one inversion. The long Weierstrass versions use

    .. math::

        \lambda_{\mathrm{chord}} = \frac{y_2 - y_1}{x_2 - x_1}

        \lambda_{\mathrm{tangent}} = \frac{3x^2 + 2 a_2 x + a_4 - a_1 y}{2y + a_1 x + a_3}

    and then `x_3 = \lambda^2 + a_1 \lambda - a_2 - x_1 - x_2`,
    `y_3 = -(\lambda + a_1) x_3 - \nu - a_3` with `\nu = y_1 - \lambda x_1`.

.. function:: int _gr_ec_jac_point_add_long_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx)
              int _gr_ec_jac_point_add_short_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx)
              int _gr_ec_jac_point_dbl_long_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
              int _gr_ec_jac_point_dbl_short_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
              int _gr_ec_jac_point_add_aff_point_long_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)
              int _gr_ec_jac_point_add_aff_point_short_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)

    Implementations of :func:`gr_ec_jac_point_add`,
    :func:`gr_ec_jac_point_dbl` and
    :func:`gr_ec_jac_point_add_aff_point`. The short Weierstrass versions
    are the classical add-1998-cmo, the `M = 3X^2 + a_4 Z^4` doubling and
    the corresponding mixed addition; the long Weierstrass ones take
    `Z_3 = H Z_1 Z_2`, respectively `Z_3 = Z D` with
    `D = 2Y + a_1 X Z + a_3 Z^3`, and reduce to the former when
    `a_1 = a_2 = a_3 = 0`.

    These write their result directly into *res* rather than through
    temporaries, so *res* is clobbered even when they fail.

.. function:: int _gr_ec_point_mul_fmpz_binary(gr_ec_point_t res, const gr_ec_point_t P, const fmpz_t n, gr_ec_ctx_t ctx)
              int _gr_ec_aff_point_mul_fmpz_binary(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const fmpz_t n, gr_ec_ctx_t ctx)
              int _gr_ec_jac_point_mul_fmpz_binary(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpz_t n, gr_ec_ctx_t ctx)

    Scalar multiplication by left-to-right binary double-and-add, staying
    in the given representation. Each of these calls the ``dbl`` and
    ``add`` functions of its representation and therefore works for any
    model for which those are implemented. Aliasing of *res* and *P* is
    allowed.

    The Jacobian one is the fallback used by
    :func:`gr_ec_jac_point_mul_fmpz` over base rings where the NAF window
    table cannot be normalized; the affine one is mostly useful for
    testing, since it performs an inversion at every step.

.. function:: int _gr_ec_jac_point_mul_fmpz_naf(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpz_t n, gr_ec_ctx_t ctx)

    Scalar multiplication by a width-`w` non-adjacent form ladder. The odd
    multiples `P, 3P, 5P, \ldots, (2^{w-1} - 1)P` are computed in Jacobian
    coordinates and then normalized to affine coordinates by a single call
    to :func:`_gr_ec_jac_point_vec_get_aff_point_vec`, so that every
    addition in the main loop is a mixed addition. The window width is
    chosen from the size of *n*.

    Returns ``GR_UNABLE`` if the table cannot be normalized, which is the
    case over a base ring that is not a field; the caller should fall back
    to :func:`_gr_ec_jac_point_mul_fmpz_binary` there.

.. function:: int _gr_ec_jac_point_vec_get_aff_point_vec(gr_ec_aff_point_struct * res, const gr_ec_jac_point_struct * P, slong len, gr_ec_ctx_t ctx)

    Converts a vector of *len* Jacobian points to affine coordinates using
    Montgomery's trick: one inversion for the whole vector, plus three
    multiplications per point. Points at infinity in the input contribute
    nothing to the product and come out as points at infinity.
