/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/*
    Thin mpz-like interface to the flint_mpn division and square root
    routines: the same semantics as the corresponding GMP mpz functions
    (signs, rounding modes, aliasing), with the limb-level work done by
    flint_mpn_tdiv_qr & co, which dispatch to GMP's mpn layer for small
    operands and to Newton iteration with FLINT's multiplication for large
    ones.
*/

#define ROUND_T 0
#define ROUND_F 1
#define ROUND_C 2

/* q = round(a / b), r = a - q b (either may be NULL), q, r not aliased
   with a, b */
static void
_flint_mpz_div_qr_noalias(mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b, int mode)
{
    mp_size_t an = FLINT_ABS(a->_mp_size), bn = FLINT_ABS(b->_mp_size);
    mp_size_t qn = an - bn + 1, rn;
    int asgn = (a->_mp_size < 0) ? -1 : 1;
    int bsgn = (b->_mp_size < 0) ? -1 : 1;
    int qsgn = asgn * bsgn;
    mp_ptr qd = NULL, rd;
    mpz_t rtmp;
    int rtmp_used = 0;

    FLINT_ASSERT(bn >= 1);

    if (an < bn)
    {
        /* truncated quotient 0, remainder a */
        if (q != NULL)
            q->_mp_size = 0;

        if (r != NULL)
        {
            mpz_set(r, a);
        }
        else if (mode != ROUND_T && an != 0)
        {
            /* need the sign of the remainder only: it is that of a */
            if ((mode == ROUND_F && asgn != bsgn) || (mode == ROUND_C && asgn == bsgn))
                mpz_set_si(q, (mode == ROUND_F) ? -1 : 1);
            return;
        }
        else
        {
            return;
        }

        goto adjust;
    }

    if (q != NULL)
        qd = FLINT_MPZ_REALLOC(q, qn + 1);   /* +1 for the rounding carry */

    if (r == NULL)
    {
        if (mode == ROUND_T)
        {
            flint_mpn_tdiv_q(qd, a->_mp_d, an, b->_mp_d, bn);
            while (qn > 0 && qd[qn - 1] == 0)
                qn--;
            q->_mp_size = (qsgn < 0) ? -qn : qn;
            return;
        }

        /* the remainder is needed to decide the rounding */
        mpz_init2(rtmp, FLINT_BITS * bn);
        rtmp_used = 1;
        r = rtmp;
    }

    rd = FLINT_MPZ_REALLOC(r, bn);

    if (q != NULL)
        _flint_mpn_tdiv_qr(qd, rd, a->_mp_d, an, b->_mp_d, bn);
    else
        flint_mpn_tdiv_r(rd, a->_mp_d, an, b->_mp_d, bn);

    if (q != NULL)
    {
        while (qn > 0 && qd[qn - 1] == 0)
            qn--;
        q->_mp_size = (qsgn < 0) ? -qn : qn;
    }

    rn = bn;
    while (rn > 0 && rd[rn - 1] == 0)
        rn--;
    r->_mp_size = (asgn < 0) ? -rn : rn;

adjust:
    /* truncated division leaves r with the sign of a; floor wants the sign
       of b, ceiling the opposite sign */
    if (mode != ROUND_T && r->_mp_size != 0)
    {
        int rsgn = (r->_mp_size < 0) ? -1 : 1;

        if (mode == ROUND_F && rsgn != bsgn)
        {
            if (q != NULL)
                mpz_sub_ui(q, q, 1);
            if (!rtmp_used)
                mpz_add(r, r, b);
        }
        else if (mode == ROUND_C && rsgn == bsgn)
        {
            if (q != NULL)
                mpz_add_ui(q, q, 1);
            if (!rtmp_used)
                mpz_sub(r, r, b);
        }
    }

    if (rtmp_used)
        mpz_clear(rtmp);
}

static void
_flint_mpz_div_qr(mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b, int mode)
{
    mpz_t qt, rt;
    mpz_ptr qq = q, rr = r;
    int qalias = (q != NULL) && (q == a || q == b);
    int ralias = (r != NULL) && (r == a || r == b);

    if (b->_mp_size == 0)
        flint_throw(FLINT_DIVZERO, "flint_mpz division by zero\n");

    /* GMP's mpz layer has dedicated fast paths for one- and two-limb
       divisors which the generic mpn code cannot match */
    if (FLINT_ABS(b->_mp_size) <= 2)
    {
        if (q != NULL && r != NULL)
        {
            if (mode == ROUND_T) mpz_tdiv_qr(q, r, a, b);
            else if (mode == ROUND_F) mpz_fdiv_qr(q, r, a, b);
            else mpz_cdiv_qr(q, r, a, b);
        }
        else if (q != NULL)
        {
            if (mode == ROUND_T) mpz_tdiv_q(q, a, b);
            else if (mode == ROUND_F) mpz_fdiv_q(q, a, b);
            else mpz_cdiv_q(q, a, b);
        }
        else
        {
            if (mode == ROUND_T) mpz_tdiv_r(r, a, b);
            else if (mode == ROUND_F) mpz_fdiv_r(r, a, b);
            else mpz_cdiv_r(r, a, b);
        }
        return;
    }

    if (qalias)
    {
        mpz_init(qt);
        qq = qt;
    }
    if (ralias)
    {
        mpz_init(rt);
        rr = rt;
    }

    _flint_mpz_div_qr_noalias(qq, rr, a, b, mode);

    /* copy rather than swap: the outputs may be fmpz-pooled mpz's, which
       must keep their own allocation */
    if (qalias)
    {
        mpz_set(q, qt);
        mpz_clear(qt);
    }
    if (ralias)
    {
        mpz_set(r, rt);
        mpz_clear(rt);
    }
}

void _flint_mpz_tdiv_qr(mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b) { _flint_mpz_div_qr(q, r, a, b, ROUND_T); }
void _flint_mpz_tdiv_q(mpz_ptr q, mpz_srcptr a, mpz_srcptr b) { _flint_mpz_div_qr(q, NULL, a, b, ROUND_T); }
void _flint_mpz_tdiv_r(mpz_ptr r, mpz_srcptr a, mpz_srcptr b) { _flint_mpz_div_qr(NULL, r, a, b, ROUND_T); }
void _flint_mpz_fdiv_qr(mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b) { _flint_mpz_div_qr(q, r, a, b, ROUND_F); }
void _flint_mpz_fdiv_q(mpz_ptr q, mpz_srcptr a, mpz_srcptr b) { _flint_mpz_div_qr(q, NULL, a, b, ROUND_F); }
void _flint_mpz_fdiv_r(mpz_ptr r, mpz_srcptr a, mpz_srcptr b) { _flint_mpz_div_qr(NULL, r, a, b, ROUND_F); }
void _flint_mpz_cdiv_qr(mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b) { _flint_mpz_div_qr(q, r, a, b, ROUND_C); }
void _flint_mpz_cdiv_q(mpz_ptr q, mpz_srcptr a, mpz_srcptr b) { _flint_mpz_div_qr(q, NULL, a, b, ROUND_C); }
void _flint_mpz_cdiv_r(mpz_ptr r, mpz_srcptr a, mpz_srcptr b) { _flint_mpz_div_qr(NULL, r, a, b, ROUND_C); }

/* r = a mod b with 0 <= r < |b| */
void
_flint_mpz_mod(mpz_ptr r, mpz_srcptr a, mpz_srcptr b)
{
    _flint_mpz_div_qr(NULL, r, a, b, (b->_mp_size < 0) ? ROUND_C : ROUND_F);
}

void
_flint_mpz_divexact(mpz_ptr q, mpz_srcptr a, mpz_srcptr b)
{
    mp_size_t an = FLINT_ABS(a->_mp_size), bn = FLINT_ABS(b->_mp_size), qn;
    int qsgn = ((a->_mp_size < 0) ? -1 : 1) * ((b->_mp_size < 0) ? -1 : 1);
    mp_ptr qd;

    if (bn == 0)
        flint_throw(FLINT_DIVZERO, "flint_mpz_divexact: division by zero\n");

    if (bn <= 2)
    {
        mpz_divexact(q, a, b);
        return;
    }

    if (an < bn)
    {
        q->_mp_size = 0;
        return;
    }

    if (q == a || q == b)
    {
        mpz_t t;
        mpz_init(t);
        _flint_mpz_divexact(t, a, b);
        mpz_set(q, t);      /* not swap: q may be an fmpz-pooled mpz */
        mpz_clear(t);
        return;
    }

    qn = an - bn + 1;
    qd = FLINT_MPZ_REALLOC(q, qn);
    flint_mpn_divexact(qd, a->_mp_d, an, b->_mp_d, bn);
    while (qn > 0 && qd[qn - 1] == 0)
        qn--;
    q->_mp_size = (qsgn < 0) ? -qn : qn;
}

/* s = floor(sqrt(a)), r = a - s^2 (r may be NULL); a >= 0 */
void
_flint_mpz_sqrtrem(mpz_ptr s, mpz_ptr r, mpz_srcptr a)
{
    mp_size_t an = a->_mp_size, sn, rn;
    mp_ptr sd, rd;

    if (an < 0)
        flint_throw(FLINT_ERROR, "flint_mpz_sqrtrem: negative input\n");

    if (an == 0)
    {
        s->_mp_size = 0;
        if (r != NULL)
            r->_mp_size = 0;
        return;
    }

    if (s == a || r == a)
    {
        mpz_t t;
        mpz_init(t);
        mpz_set(t, a);
        _flint_mpz_sqrtrem(s, r, t);
        mpz_clear(t);
        return;
    }

    sn = (an + 1) / 2;
    sd = FLINT_MPZ_REALLOC(s, sn);

    if (r == NULL)
    {
        flint_mpn_sqrtrem(sd, NULL, a->_mp_d, an);
    }
    else if (an > 2 && an < FLINT_MPN_SQRTREM_NEWTON_CUTOFF)
    {
        /* GMP writes the remainder in place given an limbs of room, so no
           temporary and copy are needed here */
        rd = FLINT_MPZ_REALLOC(r, an);
        rn = mpn_sqrtrem(sd, rd, a->_mp_d, an);
        r->_mp_size = rn;
    }
    else
    {
        rd = FLINT_MPZ_REALLOC(r, FLINT_MAX(an, sn + 1));
        rn = flint_mpn_sqrtrem(sd, rd, a->_mp_d, an);
        r->_mp_size = rn;
    }

    s->_mp_size = sn;   /* the root has exactly sn limbs */
}
