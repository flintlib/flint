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
    mpz-like interface to the flint_mpn division and square root routines:
    the same semantics as the corresponding GMP mpz functions (signs,
    rounding modes, aliasing), with the limb-level work done by
    _flint_mpn_tdiv_qr, flint_mpn_divexact and flint_mpn_sqrtrem. Each
    function is a single out-of-line function with the rounding mode and
    the wanted outputs as compile-time constants, and the rounding is done
    at the limb level, so that the overhead is no larger than that of GMP's
    mpz layer for any operand size.
*/

#define ROUND_T 0   /* truncate */
#define ROUND_F 1   /* floor */
#define ROUND_C 2   /* ceiling */
#define ROUND_M 3   /* floor with |b|: 0 <= r < |b| (mpz_mod) */

/* The case an < bn (truncated quotient 0) and division by zero, handled
   in the exported functions before tail-calling the main code, whose stack
   frame is heavier. q and r may alias a and b (but not each other). */
FLINT_FORCE_INLINE void
_mpz_div_qr_short(mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b,
    const int MODE, const int WANT_Q, const int WANT_R)
{
    mp_size_t as = a->_mp_size, bs = b->_mp_size;
    mp_size_t an = FLINT_ABS(as), bn = FLINT_ABS(bs), rn;
    mp_ptr qd, rd;
    int adjust;

    if (FLINT_UNLIKELY(bn == 0))
        flint_throw(FLINT_DIVZERO, "flint_mpz division by zero\n");

    /* truncated quotient 0, remainder a */
    adjust = (an != 0) &&
        ((MODE == ROUND_F && (as ^ bs) < 0) ||
         (MODE == ROUND_C && (as ^ bs) >= 0) ||
         (MODE == ROUND_M && as < 0));

    if (WANT_R)
    {
        if (adjust)
        {
            /* |r| = |b| - |a|; r is written before q, which may alias
               a or b */
            rd = FLINT_MPZ_REALLOC(r, bn);
            mpn_sub(rd, b->_mp_d, bn, a->_mp_d, an);
            rn = bn;
            while (rd[rn - 1] == 0)
                rn--;
            if (MODE == ROUND_F)
                r->_mp_size = (bs < 0) ? -rn : rn;
            else if (MODE == ROUND_C)
                r->_mp_size = (bs < 0) ? rn : -rn;
            else
                r->_mp_size = rn;
        }
        else if (r != a)
        {
            rd = FLINT_MPZ_REALLOC(r, an);
            flint_mpn_copyi(rd, a->_mp_d, an);
            r->_mp_size = as;
        }
    }

    if (WANT_Q)
    {
        if (adjust)
        {
            qd = FLINT_MPZ_REALLOC(q, 1);
            qd[0] = 1;
            q->_mp_size = (MODE == ROUND_F) ? -1 : 1;
        }
        else
        {
            q->_mp_size = 0;
        }
    }
}

/* q = round(a / b), r = a - q b for an >= bn >= 1, computing q only if
   WANT_Q and r only if WANT_R; q and r may alias a and b (but not each
   other) */
FLINT_FORCE_INLINE void
_mpz_div_qr_main(mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b,
    const int MODE, const int WANT_Q, const int WANT_R)
{
    mp_size_t as = a->_mp_size, bs = b->_mp_size;
    mp_size_t an = FLINT_ABS(as), bn = FLINT_ABS(bs), qn, rn;
    mp_srcptr ad, bd;
    mp_ptr qd, rd;
    int adjust, need_r;
    TMP_INIT;

    qn = an - bn + 1;
    need_r = WANT_R || MODE != ROUND_T;
    ad = a->_mp_d;
    bd = b->_mp_d;

    TMP_START;

    /* inputs overwritten by an output are copied first */
    if ((WANT_Q && q == a) || (WANT_R && r == a))
    {
        mp_ptr t = TMP_ALLOC(an * sizeof(mp_limb_t));
        flint_mpn_copyi(t, ad, an);
        ad = t;
    }
    if ((WANT_Q && q == b) || (WANT_R && r == b))
    {
        mp_ptr t = TMP_ALLOC(bn * sizeof(mp_limb_t));
        flint_mpn_copyi(t, bd, bn);
        bd = t;
    }

    if (WANT_Q)
        qd = FLINT_MPZ_REALLOC(q, qn + (MODE != ROUND_T));
    else
        qd = TMP_ALLOC(qn * sizeof(mp_limb_t));

    if (WANT_R)
        rd = FLINT_MPZ_REALLOC(r, bn);
    else if (need_r)
        rd = TMP_ALLOC(bn * sizeof(mp_limb_t));
    else
        rd = NULL;

    _flint_mpn_tdiv_qr(qd, rd, ad, an, bd, bn);

    rn = 0;
    adjust = 0;
    if (need_r)
    {
        rn = bn;
        while (rn > 0 && rd[rn - 1] == 0)
            rn--;
        adjust = (rn != 0) &&
            ((MODE == ROUND_F && (as ^ bs) < 0) ||
             (MODE == ROUND_C && (as ^ bs) >= 0) ||
             (MODE == ROUND_M && as < 0));
    }

    if (WANT_Q)
    {
        /* the top quotient limb may be zero, and only that one */
        qn -= (qd[qn - 1] == 0);

        if (adjust)
        {
            if (qn == 0)
                qd[qn++] = 1;
            else if (mpn_add_1(qd, qd, qn, 1))
                qd[qn++] = 1;
        }

        q->_mp_size = ((as ^ bs) < 0) ? -qn : qn;
    }

    if (WANT_R)
    {
        if (adjust)
        {
            /* |r| = |b| - |r| */
            mpn_sub(rd, bd, bn, rd, rn);
            rn = bn;
            while (rd[rn - 1] == 0)
                rn--;

            if (MODE == ROUND_F)
                r->_mp_size = (bs < 0) ? -rn : rn;
            else if (MODE == ROUND_C)
                r->_mp_size = (bs < 0) ? rn : -rn;
            else
                r->_mp_size = rn;
        }
        else
        {
            r->_mp_size = (MODE != ROUND_M && as < 0) ? -rn : rn;
        }
    }

    TMP_END;
}

#define DEF_DIV(name, PARAMS, Q, R, MODE, WANT_Q, WANT_R) \
FLINT_STATIC_NOINLINE void \
_mpz_##name(mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b) \
{ \
    _mpz_div_qr_main(q, r, a, b, MODE, WANT_Q, WANT_R); \
} \
void flint_mpz_##name PARAMS \
{ \
    if (FLINT_ABS(a->_mp_size) < FLINT_ABS(b->_mp_size) || b->_mp_size == 0) \
        _mpz_div_qr_short(Q, R, a, b, MODE, WANT_Q, WANT_R); \
    else \
        _mpz_##name(Q, R, a, b); \
}

DEF_DIV(tdiv_qr, (mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b), q, r, ROUND_T, 1, 1)
DEF_DIV(tdiv_q, (mpz_ptr q, mpz_srcptr a, mpz_srcptr b), q, NULL, ROUND_T, 1, 0)
DEF_DIV(tdiv_r, (mpz_ptr r, mpz_srcptr a, mpz_srcptr b), NULL, r, ROUND_T, 0, 1)
DEF_DIV(fdiv_qr, (mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b), q, r, ROUND_F, 1, 1)
DEF_DIV(fdiv_q, (mpz_ptr q, mpz_srcptr a, mpz_srcptr b), q, NULL, ROUND_F, 1, 0)
DEF_DIV(fdiv_r, (mpz_ptr r, mpz_srcptr a, mpz_srcptr b), NULL, r, ROUND_F, 0, 1)
DEF_DIV(cdiv_qr, (mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b), q, r, ROUND_C, 1, 1)
DEF_DIV(cdiv_q, (mpz_ptr q, mpz_srcptr a, mpz_srcptr b), q, NULL, ROUND_C, 1, 0)
DEF_DIV(cdiv_r, (mpz_ptr r, mpz_srcptr a, mpz_srcptr b), NULL, r, ROUND_C, 0, 1)
DEF_DIV(mod, (mpz_ptr r, mpz_srcptr a, mpz_srcptr b), NULL, r, ROUND_M, 0, 1)

FLINT_STATIC_NOINLINE void
_mpz_divexact(mpz_ptr q, mpz_srcptr a, mpz_srcptr b)
{
    mp_size_t as = a->_mp_size, bs = b->_mp_size;
    mp_size_t an = FLINT_ABS(as), bn = FLINT_ABS(bs), qn;
    mp_srcptr ad, bd;
    mp_ptr qd;
    TMP_INIT;

    qn = an - bn + 1;
    ad = a->_mp_d;
    bd = b->_mp_d;

    TMP_START;

    if (q == a)
    {
        mp_ptr t = TMP_ALLOC(an * sizeof(mp_limb_t));
        flint_mpn_copyi(t, ad, an);
        ad = t;
    }
    else if (q == b)
    {
        mp_ptr t = TMP_ALLOC(bn * sizeof(mp_limb_t));
        flint_mpn_copyi(t, bd, bn);
        bd = t;
    }

    qd = FLINT_MPZ_REALLOC(q, qn);
    flint_mpn_divexact(qd, ad, an, bd, bn);
    qn -= (qd[qn - 1] == 0);
    q->_mp_size = ((as ^ bs) < 0) ? -qn : qn;

    TMP_END;
}

void
flint_mpz_divexact(mpz_ptr q, mpz_srcptr a, mpz_srcptr b)
{
    mp_size_t an = FLINT_ABS(a->_mp_size), bn = FLINT_ABS(b->_mp_size);

    if (FLINT_UNLIKELY(bn == 0))
        flint_throw(FLINT_DIVZERO, "flint_mpz_divexact: division by zero\n");

    if (an < bn)
        q->_mp_size = 0;
    else
        _mpz_divexact(q, a, b);
}

/* s = floor(sqrt(a)), r = a - s^2 (r may be NULL); a >= 0 */
FLINT_FORCE_INLINE void
_mpz_sqrtrem(mpz_ptr s, mpz_ptr r, mpz_srcptr a)
{
    mp_size_t an = a->_mp_size, sn, rn;
    mp_srcptr ad;
    mp_ptr sd, rd;
    TMP_INIT;

    if (FLINT_UNLIKELY(an <= 0))
    {
        if (an < 0)
            flint_throw(FLINT_ERROR, "flint_mpz_sqrtrem: negative input\n");
        s->_mp_size = 0;
        if (r != NULL)
            r->_mp_size = 0;
        return;
    }

    sn = (an + 1) / 2;
    ad = a->_mp_d;

    /* inputs of at most FLINT_MPN_SQRTREM_SMALL limbs: dedicated code
       needing sn + 1 limbs of remainder space */
    if (an <= FLINT_MPN_SQRTREM_SMALL && s != a && r != a)
    {
        sd = FLINT_MPZ_REALLOC(s, sn);
        if (r != NULL)
        {
            rd = FLINT_MPZ_REALLOC(r, sn + 1);
            r->_mp_size = _flint_mpn_sqrtrem(sd, rd, ad, an);
        }
        else
        {
            _flint_mpn_sqrtrem(sd, NULL, ad, an);
        }
        s->_mp_size = sn;
        return;
    }

    TMP_START;

    if (s == a || r == a)
    {
        mp_ptr t = TMP_ALLOC(an * sizeof(mp_limb_t));
        flint_mpn_copyi(t, ad, an);
        ad = t;
    }

    sd = FLINT_MPZ_REALLOC(s, sn);

    if (r == NULL)
    {
        flint_mpn_sqrtrem(sd, NULL, ad, an);
    }
    else if (FLINT_MPN_SQRTREM_USE_GMP(an))
    {
        /* GMP writes the remainder in place given an limbs of room */
        rd = FLINT_MPZ_REALLOC(r, an);
        rn = mpn_sqrtrem(sd, rd, ad, an);
        r->_mp_size = rn;
    }
    else
    {
        rd = FLINT_MPZ_REALLOC(r, FLINT_MAX(an, sn + 1));
        rn = flint_mpn_sqrtrem(sd, rd, ad, an);
        r->_mp_size = rn;
    }

    s->_mp_size = sn;   /* the root has exactly sn limbs */

    TMP_END;
}

void flint_mpz_sqrtrem(mpz_ptr s, mpz_ptr r, mpz_srcptr a) { _mpz_sqrtrem(s, r, a); }
void flint_mpz_sqrt(mpz_ptr s, mpz_srcptr a) { _mpz_sqrtrem(s, NULL, a); }
