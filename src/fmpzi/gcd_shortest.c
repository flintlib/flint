/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"
#include "fmpq.h"

/*
return with (u1,u2) = t1*(c,b) + t2*(0,a)
                    = the l_infty shortest vector in ZZ*(c, b) + ZZ*(0, a)

!!! Note !!!
must return (u1,u2) = (1,0) or (0,1) in the case [c b; 0 a] = [1 0; 0 1]

The shortest is u = y*(c,b) - x*(0,a) for some t = (y,-x) in ZZ^2
(x/y) is either the last convergent of b/a that lies outside [(b-c)/a, (b+c)/a]
             or the first convergent that lies inside

Simple fact: If I is a finite interval containing an integer with I > 1, then
             for any z in I, either the first or second convergent of z is in I.

Outputs shouldn't alias the inputs.
*/
void _fmpz_mat22_shortest_l_infinity(
    fmpz_t u1, fmpz_t u2,
    fmpz_t t1, fmpz_t t2,
    const fmpz_t c, const fmpz_t b, const fmpz_t a)
{
    int vcmp, tries_left;
    _fmpq_cfrac_list_t s;
    _fmpz_mat22_t m;
    _fmpq_ball_t x;
    fmpz_t Q;

    FLINT_ASSERT(fmpz_cmp_si(c, 0) > 0);
    FLINT_ASSERT(fmpz_cmp(a, b) > 0);
    FLINT_ASSERT(fmpz_cmp_si(b, 0) >= 0);

    fmpz_add(u1, b, c);
    fmpz_sub(u2, b, c);

    if (fmpz_cmp(a, c) <= 0)
    {
        fmpz_zero(u1);
        fmpz_set(u2, a);
        fmpz_zero(t1);
        fmpz_one(t2);
        return;
    }
    else if (fmpz_sgn(u2) <= 0)
    {
        fmpz_set(u1, c);
        fmpz_set(u2, b);
        fmpz_one(t1);
        fmpz_zero(t2);
        return;
    }
    else if (fmpz_cmp(a, u1) <= 0)
    {
        fmpz_set(u1, c);
        fmpz_sub(u2, b, a);
        fmpz_one(t1);
        fmpz_set_si(t2, -1);
        return;
    }

    fmpz_init(Q);

    _fmpq_cfrac_list_init(s);
    s->length = -1;  /* don't need partial quotients */

    _fmpz_mat22_init(m);
    _fmpz_mat22_one(m);

    /* initialize x manually */
    fmpz_init_set(x->left_num, a);
    fmpz_init(x->left_den); fmpz_swap(x->left_den, u1);
    fmpz_init_set(x->right_num, a);
    fmpz_init(x->right_den); fmpz_swap(x->right_den, u2);
    x->exact = 0;

    _fmpq_ball_get_cfrac(s, m, 1, x);

#define v12 x->left_den
#define v22 x->left_num
#define v11 x->right_den
#define v21 x->right_num

    fmpz_add(v12, v12, v11);
    fmpz_fdiv_q_2exp(v12, v12, 1);
    fmpz_add(v22, v22, v21);
    fmpz_fdiv_q_2exp(v22, v22, 1);
    if (m->det < 0)
        fmpz_neg(v12, v12);
    else
        fmpz_neg(v22, v22);
    fmpz_mul(v11, m->_11, c);
    fmpz_mul(v21, m->_12, c);

#ifdef FLINT_WANT_ASSERT
    {
        fmpz_t tt;
        fmpz_init(tt);
        /*
            The remainder ball provides v12 and v22 easily in:
                v11 == m11*c
                v12 == m11*b - m21*a
                v21 == m12*c
                v22 == m12*b - m22*a
        */
        fmpz_mul(tt, m->_11, c);
        FLINT_ASSERT(fmpz_equal(tt, v11));
        fmpz_fmms(tt, m->_11, b, m->_21, a);
        FLINT_ASSERT(fmpz_equal(tt, v12));
        fmpz_mul(tt, m->_12, c);
        FLINT_ASSERT(fmpz_equal(tt, v21));
        fmpz_fmms(tt, m->_12, b, m->_22, a);
        FLINT_ASSERT(fmpz_equal(tt, v22));
        fmpz_clear(tt);
    }
#endif

    vcmp = fmpz_cmpabs(v11, v12);

    /*
        get_cfrac ensures I = M^-1([a/(b+c) a/(b-c)]) satisfies the simple fact.
        We have |m11*c| >= |m11*b - m21*a| iff a/(b-c) <= m11/m21 <= a/(b+c).
    */
    FLINT_ASSERT(vcmp < 0); /* since infty = M^-1(m11/m21) is outside of I. */

    /* u is best, t is transformation to best */
    fmpz_set(u1, v11);
    fmpz_set(u2, v12);
    fmpz_set(t1, m->_11);
    fmpz_neg(t2, m->_21);

    /*
        The simple fact is satisfied with z = M^-1(a/b): generate at most two
        more convergents q1 and q2. |v12| is decreasing and |v11| is increasing.
        We are done as soon as |v11| >= |v12|, since then
            a/(b+c) <= m11/21 <= a/(b-c),
        which is the same thing as q1 or q1+1/q2 in I.
    */
    tries_left = 2;
    while (--tries_left >= 0 && vcmp < 0)
    {
        fmpz_tdiv_q(Q, v22, v12);
        FLINT_ASSERT(fmpz_cmp_si(Q, 0) < 0);
        FLINT_ASSERT(fmpz_cmpabs(u1, u2) < 0);
        fmpz_submul(m->_12, m->_11, Q); fmpz_swap(m->_12, m->_11);
        fmpz_submul(m->_22, m->_21, Q); fmpz_swap(m->_22, m->_21);
        fmpz_submul(v21, v11, Q); fmpz_swap(v21, v11);
        fmpz_submul(v22, v12, Q); fmpz_swap(v22, v12);
        vcmp = fmpz_cmpabs(v11, v12);
        if (fmpz_cmpabs(vcmp < 0 ? v12 : v11, u2) < 0)
        {
            /* these could be swaps on last interation */
            fmpz_set(u1, v11);
            fmpz_set(u2, v12);
            fmpz_set(t1, m->_11);
            fmpz_neg(t2, m->_21);
        }
    }

#undef v11
#undef v12
#undef v21
#undef v22

    fmpz_clear(Q);
    _fmpq_cfrac_list_clear(s);
    _fmpz_mat22_clear(m);
    _fmpq_ball_clear(x);
}


void _fmpzi_gcd_fmpz_shortest(
    fmpz_t gx, fmpz_t gy,
    const fmpz_t ax_, const fmpz_t ay_,
    const fmpz_t b)
{
    fmpz_t A, B, C, ga, ua, va, g, u, v, axog, ayog, m1, m2, m3, m4;
    fmpz_t t, ax_copy, ay_copy;
    const fmpz* ax = ax_;
    const fmpz* ay = ay_;

    if (fmpz_is_zero(b))
    {
        fmpz_set(gx, ax);
        fmpz_set(gy, ay);
        return;
    }

    fmpz_init(t);
    fmpz_init(ax_copy);
    fmpz_init(ay_copy);

    /* reduce ax+i*ay by b */
    if (fmpz_cmpabs(ax, b) > 0)
    {
        fmpz_tdiv_qr(t, ax_copy, ax, b);
        ax = ax_copy;
    }

    if (fmpz_cmpabs(ay, b) > 0)
    {
        fmpz_tdiv_qr(t, ay_copy, ay, b);
        ay = ay_copy;
    }

    if (fmpz_is_zero(ax))
    {
        fmpz_gcd(gx, ay, b);
        fmpz_zero(gy);
        goto cleanup_stage1;
    }
    else if (fmpz_is_zero(ay))
    {
        fmpz_gcd(gx, ax, b);
        fmpz_zero(gy);
        goto cleanup_stage1;
    }

    fmpz_init(A);
    fmpz_init(B);
    fmpz_init(C);
    fmpz_init(ga);
    fmpz_init(ua);
    fmpz_init(va);
    fmpz_init(g);
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(axog);
    fmpz_init(ayog);
    fmpz_init(m1);
    fmpz_init(m2);
    fmpz_init(m3);
    fmpz_init(m4);

    fmpz_xgcd(ga, ua, va, ax, ay);
    fmpz_xgcd(g, u, v, ga, b);  /* v not needed */
    fmpz_divexact(axog, ax, g);
    fmpz_divexact(ayog, ay, g);
    fmpz_fmms(m1, ayog, ua, axog, va);
    fmpz_fmma(m2, ax, axog, ay, ayog);
    fmpz_divexact(m2, m2, ga);
    fmpz_divexact(B, b, g);
    fmpz_gcd(A, m2, B);
    fmpz_one(C);
    fmpz_mul(B, m1, u);
    fmpz_mod(B, B, A);

    _fmpz_mat22_shortest_l_infinity(gx, gy, u, v, C, B, A);
    fmpz_mul(gx, gx, g);
    fmpz_mul(gy, gy, g);

    fmpz_clear(A);
    fmpz_clear(B);
    fmpz_clear(C);
    fmpz_clear(ga);
    fmpz_clear(ua);
    fmpz_clear(va);
    fmpz_clear(g);
    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(axog);
    fmpz_clear(ayog);
    fmpz_clear(m1);
    fmpz_clear(m2);
    fmpz_clear(m3);
    fmpz_clear(m4);

cleanup_stage1:
    fmpz_clear(t);
    fmpz_clear(ax_copy);
    fmpz_clear(ay_copy);
}

/*
     ax ay           g*1 g*B            g = gcd(gx, gy)
    -ay ax  row ops   0  g*A        g^2*A = norm(gx + i*gy)
     bx by  ------>   0   0
    -by bx            0   0

    Except when the gcd is purely real or purely imaginary, it is any shortest
    vector in the above lattice in the l_infinity norm. In the exceptional case,
    it turns out that A = 1 and B = 0, and we rely on _shortest_l_infinity
    returning (1,0) or (0,1) instead of (1,1) in this case.
*/
void _fmpzi_gcd_shortest(
    fmpz_t gx, fmpz_t gy,
    const fmpz_t ax, const fmpz_t ay,
    const fmpz_t bx_, const fmpz_t by_)
{
    fmpz_t A, B, C, ag, t1, t2, bg, bu, bv, g, u, v, m1, m2, m3, m4;
    fmpz_t bx_copy, by_copy;
    const fmpz* bx = bx_;
    const fmpz* by = by_;

    /* ensure norm(a) <= norm(b) approximately */
    {
        slong ax_size = fmpz_size(ax);
        slong ay_size = fmpz_size(ay);
        slong bx_size = fmpz_size(bx);
        slong by_size = fmpz_size(by);
        if (FLINT_MAX(ax_size, ay_size) > FLINT_MAX(bx_size, by_size))
        {
            FLINT_SWAP(const fmpz *, ax, bx);
            FLINT_SWAP(const fmpz *, ay, by);
        }
    }

    if (fmpz_is_zero(ax))
    {
        _fmpzi_gcd_fmpz_shortest(gx, gy, bx, by, ay);
        return;
    }
    else if (fmpz_is_zero(ay))
    {
        _fmpzi_gcd_fmpz_shortest(gx, gy, bx, by, ax);
        return;
    }

    fmpz_init(A);
    fmpz_init(B);
    fmpz_init(C);
    fmpz_init(ag);
    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_init(bg);
    fmpz_init(bu);
    fmpz_init(bv);
    fmpz_init(g);
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(m1);
    fmpz_init(m2);
    fmpz_init(m3);
    fmpz_init(m4);
    fmpz_init(bx_copy);
    fmpz_init(by_copy);

    /*
        first find t1*ax + t2*ay = ag
         [ ax, ay]      [ag, t1*ay - t2*ax ]
         [-ay, ax] ---> [ 0, (ax^2+ay^2)/ag]
    */
    fmpz_xgcd(ag, t1, t2, ax, ay);
    fmpz_fmms(m1, t1, ay, t2, ax);
    fmpz_fmma(m2, ax, ax, ay, ay);
    fmpz_divexact(m2, m2, ag);

    /* next reduce bx+i*by by (ax^2+ay^2)/ag */
    if (fmpz_cmpabs(bx, m2) > 0)
    {
        fmpz_tdiv_qr(t1, bx_copy, bx, m2);
        bx = bx_copy;
    }

    if (fmpz_cmpabs(by, m2) > 0)
    {
        fmpz_tdiv_qr(t1, by_copy, by, m2);
        by = by_copy;
    }

    /* now reduce [bx by; -by bx] */
    if (fmpz_is_zero(bx) && fmpz_is_zero(by))
    {
        /* avoid division by bg below */
        fmpz_set(gx, ax);
        fmpz_set(gy, ay);
        goto cleanup;
    }

    fmpz_xgcd(bg, bu, bv, bx, by);
    fmpz_xgcd(g, u, v, ag, bg);

    if (!fmpz_is_one(g))
    {
        fmpz_divexact(m1, m1, g);
        fmpz_divexact(m2, m2, g);
        fmpz_divexact(t1, bx, g);
        fmpz_divexact(t2, by, g);
        fmpz_fmms(m3, bu, t2, bv, t1);
        fmpz_fmma(m4, bx, t1, by, t2);
    }
    else
    {
        fmpz_fmms(m3, bu, by, bv, bx);
        fmpz_fmma(m4, bx, bx, by, by);
    }

    fmpz_divexact(m4, m4, bg);

    /* now reduce the one remaining 2x2 */
    fmpz_fmms(t1, m3, ag, m1, bg);
    fmpz_fmma(m1, m1, u, m3, v);
    if (fmpz_is_one(g))
        fmpz_swap(m3, t1);
    else
        fmpz_divexact(m3, t1, g);

    /* reduce last column */
    fmpz_gcd3(A, m2, m3, m4);
    fmpz_fdiv_qr(t1, B, m1, A);
    fmpz_one(C);

    _fmpz_mat22_shortest_l_infinity(gx, gy, u, v, C, B, A);
    fmpz_mul(gx, gx, g);
    fmpz_mul(gy, gy, g);

cleanup:

    fmpz_clear(A);
    fmpz_clear(B);
    fmpz_clear(C);
    fmpz_clear(ag);
    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_clear(bg);
    fmpz_clear(bu);
    fmpz_clear(bv);
    fmpz_clear(g);
    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(m1);
    fmpz_clear(m2);
    fmpz_clear(m3);
    fmpz_clear(m4);
    fmpz_clear(bx_copy);
    fmpz_clear(by_copy);
}


void fmpzi_gcd_shortest(fmpzi_t g, const fmpzi_t a, const fmpzi_t b)
{
    _fmpzi_gcd_shortest(fmpzi_realref(g), fmpzi_imagref(g),
                        fmpzi_realref(a), fmpzi_imagref(a),
                        fmpzi_realref(b), fmpzi_imagref(b));
    fmpzi_canonicalise_unit(g, g);
}

