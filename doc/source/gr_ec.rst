.. _gr-ec:

**gr_ec.h** -- elliptic curves over generic rings
===============================================================================

.. note::

    This module is under construction.

A :type:`gr_ec_ctx_t` represents an elliptic curve `E` over a commutative
ring *R*, given by a Weierstrass equation

.. math::

    y^2 + a_1 x y + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6.

Points of `E` can be stored in two ways: in Jacobian coordinates
(:type:`gr_ec_jac_point_t`), which is the representation to use for
computations and works over any ring, or in affine coordinates
(:type:`gr_ec_aff_point_t`), which is convenient for input and output but
whose arithmetic needs a field. The module also provides division
polynomials and several algorithms for counting the points of a curve over
a finite field.

A curve is also a domain of the generic interface (see :ref:`gr-ec-generic`).

Over a ring that is not an integral domain, such as `\mathbb{Z}/n\mathbb{Z}`
with `n` composite, the points of `E` do not in general form a group. The
functions remain usable there, but their output only has a meaning when the
quantities they treat as nonzero are units; see
:func:`gr_ec_jac_point_mul_fmpz_witness`.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: gr_ec_ctx_struct

.. type:: gr_ec_ctx_t

    A curve. This is the same type as :type:`gr_ctx_t`; the curve data
    (a pointer to the base ring, the invariants of the curve and, once
    known, its number of points) is stored inside the context object.

.. type:: gr_ec_model_t

    The shape of the equation, which determines the formulas used for the
    group law. It takes one of the two following values.

.. macro:: GR_EC_SHORT_WEIERSTRASS

    The equation `y^2 = x^3 + a_4 x + a_6`, that is, `a_1 = a_2 = a_3 = 0`.

.. macro:: GR_EC_LONG_WEIERSTRASS

    The general equation. This is needed in characteristic 2 and 3,
    where the short equation is always singular.

.. type:: gr_ec_aff_point_struct

.. type:: gr_ec_aff_point_t

    A point in affine coordinates `(x, y)`, or the point at infinity.

.. type:: gr_ec_jac_point_struct

.. type:: gr_ec_jac_point_t

    A point in Jacobian coordinates `(X, Y, Z)`, standing for the affine
    point `(X/Z^2, Y/Z^3)`, or the point at infinity.

Both point types store their coordinates together with a :type:`truth_t`
field ``is_infinity``, which is ``T_TRUE`` for the point at infinity
`\mathcal{O}` and ``T_FALSE`` for a finite point. A function that fails may
set it to ``T_UNKNOWN``; such a point should not be used as input.

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

    Pointers to the invariants `a_1, \ldots, a_6`, `b_2, \ldots, b_8` and
    the discriminant `\Delta` stored in *ctx*. They must not be modified.

.. macro:: GR_EC_AFF_POINT_X(P, ctx)
           GR_EC_AFF_POINT_Y(P, ctx)
           GR_EC_JAC_POINT_X(P, ctx)
           GR_EC_JAC_POINT_Y(P, ctx)
           GR_EC_JAC_POINT_Z(P, ctx)

    Pointers to the coordinates of the point *P*.

Context objects
-------------------------------------------------------------------------------

.. function:: int gr_ec_ctx_init(gr_ec_ctx_t ctx, gr_ctx_t base_ring, gr_srcptr a1, gr_srcptr a2, gr_srcptr a3, gr_srcptr a4, gr_srcptr a6)
              int gr_ec_ctx_init_si(gr_ec_ctx_t ctx, gr_ctx_t base_ring, slong a1, slong a2, slong a3, slong a4, slong a6)

    Initializes *ctx* to the curve with coefficients `a_1, a_2, a_3, a_4,
    a_6` over *base_ring*. The model is set to
    :macro:`GR_EC_SHORT_WEIERSTRASS` if `a_1`, `a_2` and `a_3` are all
    zero, and to :macro:`GR_EC_LONG_WEIERSTRASS` otherwise.

    The base ring is not copied, and must remain valid until *ctx* is
    cleared.

    Returns ``GR_DOMAIN`` if the discriminant is zero (the curve is
    singular) or if *base_ring* is not a commutative ring. If the
    initialization fails, *ctx* must not be cleared.

.. function:: int gr_ec_ctx_init_short_weierstrass(gr_ec_ctx_t ctx, gr_ctx_t base_ring, gr_srcptr a4, gr_srcptr a6)
              int gr_ec_ctx_init_short_weierstrass_si(gr_ec_ctx_t ctx, gr_ctx_t base_ring, slong a4, slong a6)

    Initializes *ctx* to the curve `y^2 = x^3 + a_4 x + a_6`.

.. function:: int gr_ec_ctx_init_randtest(gr_ec_ctx_t ctx, flint_rand_t state, gr_ctx_t base_ring)

    Initializes *ctx* to a random nonsingular curve over *base_ring*, for
    test purposes. Returns ``GR_UNABLE`` if none was found.

.. function:: void gr_ec_ctx_clear(gr_ec_ctx_t ctx)

    Clears *ctx*. This is the same as :func:`gr_ctx_clear`.

.. function:: gr_ctx_struct * gr_ec_ctx_base_ring(gr_ec_ctx_t ctx)

    Returns the base ring of *ctx*.

.. function:: gr_ec_model_t gr_ec_ctx_model(gr_ec_ctx_t ctx)

    Returns the model of *ctx*.

Invariants
-------------------------------------------------------------------------------

The notation is that of Silverman, *The Arithmetic of Elliptic Curves*,
Section III.1.

.. function:: int gr_ec_ctx_a_invariants(gr_ptr a1, gr_ptr a2, gr_ptr a3, gr_ptr a4, gr_ptr a6, gr_ec_ctx_t ctx)

    Sets the outputs to the coefficients of the equation of *ctx*.

.. function:: int gr_ec_ctx_b_invariants(gr_ptr b2, gr_ptr b4, gr_ptr b6, gr_ptr b8, gr_ec_ctx_t ctx)

    Sets the outputs to `b_2 = a_1^2 + 4 a_2`, `b_4 = 2 a_4 + a_1 a_3`,
    `b_6 = a_3^2 + 4 a_6` and
    `b_8 = a_1^2 a_6 + 4 a_2 a_6 - a_1 a_3 a_4 + a_2 a_3^2 - a_4^2`.

.. function:: int gr_ec_ctx_c_invariants(gr_ptr c4, gr_ptr c6, gr_ec_ctx_t ctx)

    Sets the outputs to `c_4 = b_2^2 - 24 b_4` and
    `c_6 = -b_2^3 + 36 b_2 b_4 - 216 b_6`.

.. function:: int gr_ec_ctx_discriminant(gr_ptr res, gr_ec_ctx_t ctx)

    Sets *res* to the discriminant
    `\Delta = -b_2^2 b_8 - 8 b_4^3 - 27 b_6^2 + 9 b_2 b_4 b_6`.

.. function:: int gr_ec_ctx_j_invariant(gr_ptr res, gr_ec_ctx_t ctx)

    Sets *res* to the `j`-invariant `c_4^3 / \Delta`. Returns
    ``GR_DOMAIN`` if `\Delta` is not invertible.

.. function:: truth_t gr_ec_ctx_is_smooth(gr_ec_ctx_t ctx)

    Returns whether the discriminant of *ctx* is nonzero. Since a
    singular curve cannot be created, this never returns ``T_FALSE``, but
    it returns ``T_UNKNOWN`` when the base ring cannot decide.

.. _gr-ec-generic:

Generic interface
-------------------------------------------------------------------------------

A :type:`gr_ec_ctx_t` is also a :type:`gr_ctx_t`, whose elements are the
points of the curve, stored as :type:`gr_ec_jac_point_t`. The domain is the
abelian group `E(R)`, written additively: :func:`gr_zero` is the point at
infinity, :func:`gr_add`, :func:`gr_sub` and :func:`gr_neg` are the group
law, and :func:`gr_mul_ui`, :func:`gr_mul_si`, :func:`gr_mul_fmpz`,
:func:`gr_mul_two` and :func:`gr_mul_2exp_si` compute integer multiples of a
point. :func:`gr_ctx_cardinality_fmpz` returns the number of points.

Operations that make sense only in a ring, such as :func:`gr_mul`,
:func:`gr_one` or :func:`gr_set_si`, return ``GR_DOMAIN``. Points must be
created with the functions of this module, for example
:func:`gr_ec_jac_point_set_affine` or :func:`gr_ec_jac_point_lift_x`.

Affine points
-------------------------------------------------------------------------------

The functions computing the group law in affine coordinates need to divide,
and return ``GR_DOMAIN`` when the base ring is not a field.

.. function:: void gr_ec_aff_point_init(gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Initializes *P* to the point at infinity.

.. function:: void gr_ec_aff_point_clear(gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

.. function:: void gr_ec_aff_point_swap(gr_ec_aff_point_t P, gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)

.. function:: int gr_ec_aff_point_set(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

.. function:: int gr_ec_aff_point_zero(gr_ec_aff_point_t res, gr_ec_ctx_t ctx)

    Sets *res* to the point at infinity.

.. function:: int _gr_ec_aff_point_set_affine(gr_ec_aff_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx)
              int gr_ec_aff_point_set_affine(gr_ec_aff_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx)

    Sets *res* to the point `(x, y)`. The non-underscore version returns
    ``GR_DOMAIN`` if the point is not on the curve.

.. function:: int gr_ec_aff_point_get_affine(gr_ptr x, gr_ptr y, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Sets *x* and *y* to the coordinates of *P*. Returns ``GR_DOMAIN`` if
    *P* is the point at infinity.

.. function:: int gr_ec_aff_point_lift_x(gr_ec_aff_point_t res, gr_srcptr x, gr_ec_ctx_t ctx)

    Sets *res* to one of the points with `x`-coordinate *x*. Returns
    ``GR_DOMAIN`` if there is no such point. For a curve which is not in
    short Weierstrass form, this requires 2 to be invertible.

.. function:: truth_t gr_ec_aff_point_is_inf(const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

.. function:: truth_t gr_ec_aff_point_equal(const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)

.. function:: truth_t gr_ec_aff_point_is_on_curve(const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

.. function:: int gr_ec_aff_point_write(gr_stream_t out, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
              int gr_ec_aff_point_print(const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Writes *P* as ``(x, y)``, or ``O`` for the point at infinity.

.. function:: int gr_ec_aff_point_randtest(gr_ec_aff_point_t res, flint_rand_t state, gr_ec_ctx_t ctx)

    Sets *res* to a random point of the curve, for test purposes. Returns
    ``GR_UNABLE`` if no point was found.

.. function:: int gr_ec_aff_point_neg(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
              int gr_ec_aff_point_add(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)
              int gr_ec_aff_point_sub(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)
              int gr_ec_aff_point_dbl(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to `-P`, `P + Q`, `P - Q` and `2P` respectively.

.. function:: int gr_ec_aff_point_mul_ui(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, ulong n, gr_ec_ctx_t ctx)
              int gr_ec_aff_point_mul_si(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, slong n, gr_ec_ctx_t ctx)
              int gr_ec_aff_point_mul_fmpz(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const fmpz_t n, gr_ec_ctx_t ctx)

    Sets *res* to `nP`. The computation is done in Jacobian coordinates.

Jacobian points
-------------------------------------------------------------------------------

The functions in this section work over any commutative ring.

.. function:: void gr_ec_jac_point_init(gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Initializes *P* to the point at infinity.

.. function:: void gr_ec_jac_point_clear(gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

.. function:: void gr_ec_jac_point_swap(gr_ec_jac_point_t P, gr_ec_jac_point_t Q, gr_ec_ctx_t ctx)

.. function:: int gr_ec_jac_point_set(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

.. function:: int gr_ec_jac_point_zero(gr_ec_jac_point_t res, gr_ec_ctx_t ctx)

    Sets *res* to the point at infinity.

.. function:: int _gr_ec_jac_point_set_affine(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_set_affine(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y, gr_ec_ctx_t ctx)

    Sets *res* to the point with affine coordinates `(x, y)`. The
    non-underscore version returns ``GR_DOMAIN`` if the point is not on
    the curve.

.. function:: int _gr_ec_jac_point_set_jacobian(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_set_jacobian(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ec_ctx_t ctx)

    Sets *res* to the point with Jacobian coordinates `(x, y, z)`. The
    non-underscore version returns ``GR_DOMAIN`` if `z` is zero or if the
    point is not on the curve.

.. function:: int gr_ec_jac_point_get_affine(gr_ptr x, gr_ptr y, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Sets *x* and *y* to the affine coordinates of *P*. Returns
    ``GR_DOMAIN`` if *P* is the point at infinity or if its `Z`-coordinate
    is not invertible.

.. function:: int gr_ec_jac_point_normalize(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to the representative of *P* with `Z = 1`. Fails like
    :func:`gr_ec_jac_point_get_affine`, except that the point at infinity
    is allowed.

.. function:: int gr_ec_jac_point_lift_x(gr_ec_jac_point_t res, gr_srcptr x, gr_ec_ctx_t ctx)

    Sets *res* to one of the points with affine `x`-coordinate *x*, as
    :func:`gr_ec_aff_point_lift_x`.

.. function:: truth_t gr_ec_jac_point_is_inf(const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

.. function:: truth_t gr_ec_jac_point_equal(const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx)

    Returns whether *P* and *Q* represent the same point. The answer is only
    reliable over an integral domain.

.. function:: truth_t gr_ec_jac_point_is_on_curve(const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

.. function:: int gr_ec_jac_point_write(gr_stream_t out, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_print(const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Writes *P* as ``(X, Y, Z)``, or ``O`` for the point at infinity.

.. function:: int gr_ec_jac_point_randtest(gr_ec_jac_point_t res, flint_rand_t state, gr_ec_ctx_t ctx)

    Sets *res* to a random point of the curve, for test purposes. Returns
    ``GR_UNABLE`` if no point was found.

.. function:: int gr_ec_jac_point_neg(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_add(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_sub(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_dbl(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to `-P`, `P + Q`, `P - Q` and `2P` respectively.

.. function:: int gr_ec_jac_point_add_aff_point(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_sub_aff_point(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx)

    Sets *res* to `P + Q` or `P - Q` where *Q* is an affine point. This is
    faster than first converting *Q* to Jacobian coordinates.

.. function:: int gr_ec_jac_point_mul_ui(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, ulong n, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_mul_si(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, slong n, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_mul_fmpz(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpz_t n, gr_ec_ctx_t ctx)

    Sets *res* to `nP`. If the number of points of the curve is stored in
    *ctx* (see :func:`gr_ec_ctx_cardinality`), *n* is first reduced modulo
    it.

.. function:: int gr_ec_jac_point_mul_2exp_si(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, slong k, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_mul_2exp_fmpz(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpz_t k, gr_ec_ctx_t ctx)

    Sets *res* to `2^k P`. Returns ``GR_DOMAIN`` if `k` is negative, and
    ``GR_UNABLE`` if `k` does not fit in a *slong*.

.. function:: int gr_ec_jac_point_mul_fmpz_witness(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpz_t n, gr_ptr w, gr_ec_ctx_t ctx)

    Sets *res* to `nP`, and multiplies *w* by every value that the group
    law needed to be invertible along the way.

    This is meant for rings such as `\mathbb{Z}/m\mathbb{Z}` with `m` not
    known to be prime. The result is correct when *w*, initially set to
    one, is a unit at the end; otherwise `\gcd(w, m)` is a nontrivial
    factor of `m`. Elliptic curve primality proving (:ref:`ecpp`) relies on
    this.

Conversions
-------------------------------------------------------------------------------

.. function:: int gr_ec_jac_point_set_aff_point(gr_ec_jac_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to *P*, with `Z = 1`.

.. function:: int gr_ec_aff_point_set_jac_point(gr_ec_aff_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)

    Sets *res* to *P*. Fails like :func:`gr_ec_jac_point_get_affine`,
    except that the point at infinity is allowed.

Division by an integer
-------------------------------------------------------------------------------

Given a point `P` and an integer `n`, the points `Q` with `nQ = P` are
either absent or as many as the `n`-torsion points of the curve. Two
variants are provided: one that wants a unique answer and one that accepts
any of them.

When `n` is coprime to the number of points of a curve over a finite field,
the answer is found with a single multiplication, after counting the points
if needed. Otherwise the solutions are found among the roots of a
polynomial of degree about `n^2`; this is only done for the short
Weierstrass model and for `n` up to 64, and ``GR_UNABLE`` is returned in
other cases.

.. function:: int gr_ec_jac_point_div_fmpz(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpz_t n, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_div_ui(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, ulong n, gr_ec_ctx_t ctx)

    Sets *res* to the unique point `Q` with `nQ = P`. Returns ``GR_DOMAIN``
    if there is no such point or if there are several.

.. function:: int gr_ec_jac_point_div_fmpz_nonunique(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpz_t n, gr_ec_ctx_t ctx)

    Sets *res* to some point `Q` with `nQ = P`. Returns ``GR_DOMAIN`` if
    there is none.

.. function:: int gr_ec_jac_point_mul_fmpq(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpq_t c, gr_ec_ctx_t ctx)
              int gr_ec_jac_point_mul_fmpq_nonunique(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpq_t c, gr_ec_ctx_t ctx)

    Writing `c = a/b`, sets *res* to a point `Q` with `bQ = aP`, with the
    same conventions as above.

Division polynomials
-------------------------------------------------------------------------------

The division polynomials `\psi_n` describe the `n`-torsion of `E`: a finite
point `P` satisfies `nP = \mathcal{O}` exactly when `\psi_n(P) = 0`. When `n`
is odd, `\psi_n` is a polynomial in `x`; when `n` is even, it is `\psi_2 =
2y + a_1 x + a_3` times a polynomial in `x`. The functions below return the
polynomial in `x`, that is

.. math::

    \Psi_n = \psi_n \quad (n \text{ odd}), \qquad
    \Psi_n = \psi_n / \psi_2 \quad (n \text{ even}).

In particular `\Psi_0 = 0` and `\Psi_1 = \Psi_2 = 1`. Note that PARI/GP's
``elldivpol`` returns `\Psi_n \psi_2^2` for even `n`.

.. function:: int gr_ec_ctx_psi2_sqr(gr_poly_t res, gr_ec_ctx_t ctx)

    Sets *res* to `\psi_2^2 = 4x^3 + b_2 x^2 + 2 b_4 x + b_6`.

.. function:: int gr_ec_ctx_division_poly(gr_poly_t res, ulong n, gr_ec_ctx_t ctx)

    Sets *res* to `\Psi_n`.

.. function:: int gr_ec_ctx_division_poly_vec(gr_poly_struct * res, slong len, gr_ec_ctx_t ctx)

    Sets the entries of *res* to `\Psi_0, \ldots, \Psi_{len-1}`. This is
    faster than computing them one by one.

Point counting
-------------------------------------------------------------------------------

The functions in this section compute the number of points
`\#E(\mathbb{F}_q) = q + 1 - t` of a curve over a finite field, where
`|t| \le 2 \sqrt{q}`. They return ``GR_DOMAIN`` if the base ring is not
known to be a field (see :func:`gr_ctx_set_is_field`), except
:func:`gr_ec_ctx_cardinality` which also handles
`\mathbb{Z}/n\mathbb{Z}` for composite `n`.

The number of points is stored in the context by
:func:`gr_ec_ctx_cardinality` and :func:`gr_ec_ctx_set_cardinality`, and is
then used to speed up scalar multiplication and division. These two
functions modify *ctx*, so they should not be called while other threads
are using the same curve.

.. function:: int gr_ec_ctx_cardinality(fmpz_t res, gr_ec_ctx_t ctx)

    Sets *res* to the number of points of the curve, choosing an algorithm
    automatically, and stores the result in *ctx*. If the number is
    already stored, it is returned directly.

    Over `\mathbb{Z}/n\mathbb{Z}` with `n` not known to be prime, this
    uses :func:`gr_ec_ctx_cardinality_crt`.

.. function:: int gr_ec_ctx_set_cardinality(gr_ec_ctx_t ctx, const fmpz_t N)

    Stores *N* as the number of points of the curve. Returns ``GR_DOMAIN``
    if *N* is found to be wrong, which is checked against the Hasse bound,
    by counting over small fields, and by multiplying random points by *N*.
    Since these checks cannot prove that *N* is correct, it is the
    caller's responsibility to provide the right value.

.. function:: int gr_ec_ctx_cardinality_naive(fmpz_t res, gr_ec_ctx_t ctx)

    Counts the points by going through all `x \in \mathbb{F}_q`. Returns
    ``GR_UNABLE`` if `q` is larger than `10^8`.

.. function:: int gr_ec_ctx_cardinality_bsgs(fmpz_t res, gr_ec_ctx_t ctx)

    Uses the baby-step giant-step algorithm of Shanks and Mestre
    ([Coh1996]_, Section 7.4), which takes `O(q^{1/4})` group operations.
    Returns ``GR_UNABLE`` if the order could not be determined this way,
    which can happen over very small fields.

.. function:: int gr_ec_ctx_cardinality_schoof(fmpz_t res, gr_ec_ctx_t ctx)
              int gr_ec_ctx_cardinality_sea(fmpz_t res, gr_ec_ctx_t ctx)

    Uses Schoof's algorithm, respectively its improvement by Elkies and
    Atkin, whose running time is polynomial in `\log q`. Only curves in
    short Weierstrass form in characteristic greater than 3 are supported;
    ``GR_DOMAIN`` is returned otherwise.

.. function:: int gr_ec_ctx_cardinality_cm(fmpz_t res, gr_ec_ctx_t ctx)

    Counts the points of a curve in short Weierstrass form over a prime
    field `\mathbb{F}_p`, `p > 3`, whose `j`-invariant is that of a curve
    with complex multiplication by an imaginary quadratic order of class
    number at most 4. Returns ``GR_UNABLE`` for other curves.

.. function:: int gr_ec_ctx_cardinality_subfield(fmpz_t res, gr_ec_ctx_t ctx)

    Counts the points of a curve over `\mathbb{F}_{p^k}`, `k > 1`, whose
    coefficients lie in `\mathbb{F}_p`, by counting over `\mathbb{F}_p`.
    Returns ``GR_UNABLE`` for other curves.

.. function:: int gr_ec_ctx_cardinality_crt(fmpz_t res, gr_ec_ctx_t ctx)

    Counts the points of a curve over `\mathbb{Z}/n\mathbb{Z}` (the base
    ring must be of type *nmod*, *fmpz_mod* or *mpn_mod*) by counting
    modulo each prime factor of `n`. Returns ``GR_DOMAIN`` if the
    discriminant is not a unit modulo `n`, and ``GR_UNABLE`` if `n` could
    not be factored.
