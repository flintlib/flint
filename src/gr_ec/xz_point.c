/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Montgomery curves B y^2 = x^3 + A x^2 + x, with points carried as
    (X : Z) and no y at all.

    Dropping y costs the ability to add two arbitrary points -- the sum of
    P and Q depends on y, and only the pair {P+Q, P-Q} is determined by
    x(P) and x(Q) -- and buys three things. Doubling and *differential*
    addition, where x(P-Q) is supplied, need no inversion and fewer
    multiplications than any Weierstrass formula; the ladder built from
    them does the same work whatever the bits of the scalar are, so it
    leaks nothing through timing; and nothing in it ever branches on
    whether a quantity vanishes, which is what makes it the right
    arithmetic over Z/n for a composite n. Elliptic curve factorization
    wants all three.

    These functions take the base ring and the curve constant directly
    rather than a gr_ec_ctx_t. An (X : Z) pair is not an element of a
    group, so it is not a gr domain, and the only curve datum the formulas
    read is a24 = (A + 2)/4.

    The point at infinity is Z = 0. In x-only coordinates it is not
    distinguishable from anything by its x, so is_zero is the only
    predicate that means what it says.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_ec.h"

#define XX(P) ((P)->coords)
#define XZ(P) GR_ENTRY((P)->coords, 1, R->sizeof_elem)

void
gr_ec_xz_point_init(gr_ec_xz_point_t P, gr_ctx_t R)
{
    P->coords = gr_heap_init_vec(2, R);
}

void
gr_ec_xz_point_clear(gr_ec_xz_point_t P, gr_ctx_t R)
{
    gr_heap_clear_vec(P->coords, 2, R);
    P->coords = NULL;
}

int
gr_ec_xz_point_set(gr_ec_xz_point_t res, const gr_ec_xz_point_t P, gr_ctx_t R)
{
    int status = GR_SUCCESS;

    status |= gr_set(XX(res), XX(P), R);
    status |= gr_set(XZ(res), XZ(P), R);

    return status;
}

void
gr_ec_xz_point_swap(gr_ec_xz_point_t P, gr_ec_xz_point_t Q, gr_ctx_t FLINT_UNUSED(R))
{
    FLINT_SWAP(gr_ec_xz_point_struct, *P, *Q);
}

/* the point at infinity, (1 : 0) */
int
gr_ec_xz_point_zero(gr_ec_xz_point_t res, gr_ctx_t R)
{
    int status = GR_SUCCESS;

    status |= gr_one(XX(res), R);
    status |= gr_zero(XZ(res), R);

    return status;
}

truth_t
gr_ec_xz_point_is_zero(const gr_ec_xz_point_t P, gr_ctx_t R)
{
    return gr_is_zero(XZ(P), R);
}

int
gr_ec_xz_point_set_x(gr_ec_xz_point_t res, gr_srcptr x, gr_ctx_t R)
{
    int status = GR_SUCCESS;

    status |= gr_set(XX(res), x, R);
    status |= gr_one(XZ(res), R);

    return status;
}

int
gr_ec_xz_point_get_x(gr_ptr x, const gr_ec_xz_point_t P, gr_ctx_t R)
{
    return gr_div(x, XX(P), XZ(P), R);
}

/* x(P) = x(Q) as projective pairs, so X1 Z2 = X2 Z1 */
truth_t
gr_ec_xz_point_equal(const gr_ec_xz_point_t P, const gr_ec_xz_point_t Q,
        gr_ctx_t R)
{
    gr_ptr u, v;
    truth_t res;

    GR_TMP_INIT2(u, v, R);

    if (gr_mul(u, XX(P), XZ(Q), R) != GR_SUCCESS
            || gr_mul(v, XX(Q), XZ(P), R) != GR_SUCCESS)
        res = T_UNKNOWN;
    else
        res = gr_equal(u, v, R);

    GR_TMP_CLEAR2(u, v, R);

    return res;
}

/*
    A Montgomery curve is a Weierstrass curve in disguise:

        B v^2 = u^3 + A u^2 + u

    becomes y^2 = x^3 + a x + b under x = (u + A/3)/B, v = y/B, with

        a = (3 - A^2)/(3 B^2),   b = (2 A^3 - 9 A)/(27 B^3).

    That map is why the Montgomery form does not need to be a model of its
    own in this module: everything else here -- point counting, division
    polynomials, the generic interface, the group law itself -- is reached
    by building the Weierstrass curve and carrying x-coordinates across,
    while the x-only arithmetic above is used where its speed and its
    indifference to zero divisors are what matter.
*/
int
gr_ec_montgomery_x_to_weierstrass(gr_ptr x, gr_srcptr u, gr_srcptr A,
        gr_srcptr B, gr_ctx_t R)
{
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, R);

    status |= gr_set_ui(t, 3, R);
    status |= gr_div(t, A, t, R);           /* A/3 */
    status |= gr_add(x, u, t, R);

    if (status == GR_SUCCESS)
        status = gr_div(x, x, B, R);

    GR_TMP_CLEAR(t, R);

    return status;
}

int
gr_ec_weierstrass_x_to_montgomery(gr_ptr u, gr_srcptr x, gr_srcptr A,
        gr_srcptr B, gr_ctx_t R)
{
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, R);

    status |= gr_set_ui(t, 3, R);
    status |= gr_div(t, A, t, R);           /* A/3 */
    status |= gr_mul(u, x, B, R);
    status |= gr_sub(u, u, t, R);

    GR_TMP_CLEAR(t, R);

    return status;
}

int
gr_ec_ctx_init_from_montgomery(gr_ec_ctx_t ctx, gr_ctx_t R, gr_srcptr A,
        gr_srcptr B)
{
    gr_ptr a, b, t, u, three;
    int status = GR_SUCCESS;

    GR_TMP_INIT5(a, b, t, u, three, R);

    status |= gr_set_ui(three, 3, R);

    /* a = (3 - A^2) / (3 B^2) */
    status |= gr_sqr(t, A, R);
    status |= gr_sub(a, three, t, R);
    status |= gr_sqr(u, B, R);
    status |= gr_mul(u, u, three, R);

    if (status == GR_SUCCESS)
        status = gr_div(a, a, u, R);

    /* b = (2 A^3 - 9 A) / (27 B^3) */
    if (status == GR_SUCCESS)
    {
        status |= gr_mul(t, t, A, R);           /* A^3 */
        status |= gr_mul_two(b, t, R);
        status |= gr_mul_ui(t, A, 9, R);
        status |= gr_sub(b, b, t, R);

        status |= gr_sqr(u, B, R);
        status |= gr_mul(u, u, B, R);
        status |= gr_mul_ui(u, u, 27, R);

        if (status == GR_SUCCESS)
            status = gr_div(b, b, u, R);
    }

    if (status == GR_SUCCESS)
        status = gr_ec_ctx_init_short_weierstrass(ctx, R, a, b);

    GR_TMP_CLEAR5(a, b, t, u, three, R);

    return status;
}

int
gr_ec_montgomery_a24(gr_ptr a24, gr_srcptr A, gr_ctx_t R)
{
    gr_ptr four;
    int status = GR_SUCCESS;

    GR_TMP_INIT(four, R);

    status |= gr_set_ui(four, 4, R);
    status |= gr_set_ui(a24, 2, R);
    status |= gr_add(a24, a24, A, R);

    if (status == GR_SUCCESS)
        status = gr_div(a24, a24, four, R);

    GR_TMP_CLEAR(four, R);

    return status;
}

/*
    2P from P.

      X2 = (X + Z)^2 (X - Z)^2
      Z2 = 4 X Z ((X - Z)^2 + a24 * 4 X Z)

    with 4 X Z = (X + Z)^2 - (X - Z)^2. Two squarings, three
    multiplications, and no test of any kind.
*/
int
gr_ec_xz_point_dbl(gr_ec_xz_point_t res, const gr_ec_xz_point_t P,
        gr_srcptr a24, gr_ctx_t R)
{
    gr_ptr t0, t1, t2;
    int status = GR_SUCCESS;

    GR_TMP_INIT3(t0, t1, t2, R);

    status |= gr_add(t0, XX(P), XZ(P), R);
    status |= gr_sqr(t0, t0, R);
    status |= gr_sub(t1, XX(P), XZ(P), R);
    status |= gr_sqr(t1, t1, R);
    status |= gr_sub(t2, t0, t1, R);             /* 4 X Z */

    status |= gr_mul(XX(res), t0, t1, R);

    status |= gr_mul(t0, t2, a24, R);
    status |= gr_add(t0, t0, t1, R);
    status |= gr_mul(XZ(res), t2, t0, R);

    GR_TMP_CLEAR3(t0, t1, t2, R);

    return status;
}

/*
    P + Q from P, Q and P - Q.

      X3 = Z_ ((X1 - Z1)(X2 + Z2) + (X1 + Z1)(X2 - Z2))^2
      Z3 = X_ ((X1 - Z1)(X2 + Z2) - (X1 + Z1)(X2 - Z2))^2

    where (X_ : Z_) is P - Q. Four multiplications and two squarings, and
    again nothing is tested.

    The formula needs x(P - Q) to be neither zero nor infinity, which is
    visible in it: X_ = 0 forces Z3 = 0 and Z_ = 0 forces X3 = 0, so the
    answer degenerates to (0 : 0) and says nothing. The caller is
    responsible for that; the ladder below keeps P - Q equal to its own
    argument throughout and deals with the two bad values of it up front.
*/
/*
    The same, for a difference already scaled to Z = 1, which saves the
    multiplication by it. The ladder differentially adds against the same
    point at every step, so when that point arrives as an affine
    x-coordinate -- which is how gr_ec_xz_point_set_x leaves it, and how
    callers normally have it -- this is one multiplication in seven saved
    on the whole scalar multiplication.
*/
static int
_gr_ec_xz_point_dadd_z1(gr_ec_xz_point_t res, const gr_ec_xz_point_t P,
        const gr_ec_xz_point_t Q, gr_srcptr xPmQ, gr_ctx_t R)
{
    gr_ptr t0, t1, t2, t3;
    int status = GR_SUCCESS;

    GR_TMP_INIT4(t0, t1, t2, t3, R);

    status |= gr_sub(t0, XX(P), XZ(P), R);
    status |= gr_add(t1, XX(Q), XZ(Q), R);
    status |= gr_mul(t0, t0, t1, R);

    status |= gr_add(t1, XX(P), XZ(P), R);
    status |= gr_sub(t2, XX(Q), XZ(Q), R);
    status |= gr_mul(t1, t1, t2, R);

    status |= gr_add(t2, t0, t1, R);
    status |= gr_sqr(t2, t2, R);
    status |= gr_sub(t3, t0, t1, R);
    status |= gr_sqr(t3, t3, R);

    status |= gr_mul(t3, t3, xPmQ, R);

    status |= gr_set(XX(res), t2, R);
    status |= gr_set(XZ(res), t3, R);

    GR_TMP_CLEAR4(t0, t1, t2, t3, R);

    return status;
}

int
gr_ec_xz_point_dadd(gr_ec_xz_point_t res, const gr_ec_xz_point_t P,
        const gr_ec_xz_point_t Q, const gr_ec_xz_point_t PmQ, gr_ctx_t R)
{
    gr_ptr t0, t1, t2, t3;
    int status = GR_SUCCESS;

    GR_TMP_INIT4(t0, t1, t2, t3, R);

    status |= gr_sub(t0, XX(P), XZ(P), R);
    status |= gr_add(t1, XX(Q), XZ(Q), R);
    status |= gr_mul(t0, t0, t1, R);            /* (X1-Z1)(X2+Z2) */

    status |= gr_add(t1, XX(P), XZ(P), R);
    status |= gr_sub(t2, XX(Q), XZ(Q), R);
    status |= gr_mul(t1, t1, t2, R);            /* (X1+Z1)(X2-Z2) */

    status |= gr_add(t2, t0, t1, R);
    status |= gr_sqr(t2, t2, R);
    status |= gr_sub(t3, t0, t1, R);
    status |= gr_sqr(t3, t3, R);

    /* res may alias P, Q or PmQ, so both products are formed first */
    status |= gr_mul(t2, t2, XZ(PmQ), R);
    status |= gr_mul(t3, t3, XX(PmQ), R);

    status |= gr_set(XX(res), t2, R);
    status |= gr_set(XZ(res), t3, R);

    GR_TMP_CLEAR4(t0, t1, t2, t3, R);

    return status;
}

/*
    k P by the Montgomery ladder.

    Two points are carried with the invariant that their difference is
    always P, which is what lets the differential addition be used: at
    every bit one of them is doubled and the other becomes their sum, so
    the work per bit does not depend on the bit.
*/
int
gr_ec_xz_point_mul_fmpz(gr_ec_xz_point_t res, const gr_ec_xz_point_t P,
        const fmpz_t k, gr_srcptr a24, gr_ctx_t R)
{
    gr_ec_xz_point_t R0, R1, base;
    fmpz_t n;
    slong i, bits;
    int z1, status = GR_SUCCESS;

    if (fmpz_is_zero(k))
        return gr_ec_xz_point_zero(res, R);

    /*
        The ladder differentially adds against P at every step, so it needs
        x(P) to be neither infinity nor zero. Both exceptions have an
        answer that needs no ladder at all: the first is fixed by
        everything, and the second is the point of order two at the origin,
        which k either fixes or kills according to its parity.
    */
    {
        truth_t inf = gr_ec_xz_point_is_zero(P, R);
        truth_t two;

        if (inf == T_UNKNOWN)
            return GR_UNABLE;

        if (inf == T_TRUE)
            return gr_ec_xz_point_zero(res, R);

        two = gr_is_zero(XX(P), R);

        if (two == T_UNKNOWN)
            return GR_UNABLE;

        if (two == T_TRUE)
            return fmpz_is_even(k) ? gr_ec_xz_point_zero(res, R)
                                   : gr_ec_xz_point_set(res, P, R);
    }

    fmpz_init(n);
    fmpz_abs(n, k);                 /* x(-P) = x(P), so the sign is nothing */

    bits = fmpz_bits(n);

    if (bits == 1)
    {
        status = gr_ec_xz_point_set(res, P, R);
        fmpz_clear(n);
        return status;
    }

    gr_ec_xz_point_init(R0, R);
    gr_ec_xz_point_init(R1, R);
    gr_ec_xz_point_init(base, R);

    status |= gr_ec_xz_point_set(base, P, R);
    status |= gr_ec_xz_point_set(R0, P, R);
    status |= gr_ec_xz_point_dbl(R1, P, a24, R);

    /* the difference is the same point at every step, so whether it is
       already scaled is worth asking once rather than per step */
    z1 = (gr_is_one(XZ(base), R) == T_TRUE);

    for (i = bits - 2; i >= 0 && status == GR_SUCCESS; i--)
    {
        if (fmpz_tstbit(n, i))
        {
            status |= z1 ? _gr_ec_xz_point_dadd_z1(R0, R0, R1, XX(base), R)
                         : gr_ec_xz_point_dadd(R0, R0, R1, base, R);
            status |= gr_ec_xz_point_dbl(R1, R1, a24, R);
        }
        else
        {
            status |= z1 ? _gr_ec_xz_point_dadd_z1(R1, R0, R1, XX(base), R)
                         : gr_ec_xz_point_dadd(R1, R0, R1, base, R);
            status |= gr_ec_xz_point_dbl(R0, R0, a24, R);
        }
    }

    if (status == GR_SUCCESS)
        status = gr_ec_xz_point_set(res, R0, R);

    gr_ec_xz_point_clear(base, R);
    gr_ec_xz_point_clear(R1, R);
    gr_ec_xz_point_clear(R0, R);
    fmpz_clear(n);

    return status;
}

int
gr_ec_xz_point_mul_ui(gr_ec_xz_point_t res, const gr_ec_xz_point_t P,
        ulong k, gr_srcptr a24, gr_ctx_t R)
{
    fmpz_t t;
    int status;

    fmpz_init_set_ui(t, k);
    status = gr_ec_xz_point_mul_fmpz(res, P, t, a24, R);
    fmpz_clear(t);

    return status;
}
