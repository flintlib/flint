/*
    From GMP's mpz/gcd.c, mpz/gcdext.c, mpz/invert.c:
    Copyright 1991, 1993, 1994, 1996, 2000-2002, 2005, 2010, 2011, 2012,
    2015 Free Software Foundation, Inc.

    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "impl.h"

/*
    mpz-like interface to flint_mpn_gcd and _flint_mpn_gcdext_inplace, with
    the semantics of GMP's mpz_gcd, mpz_gcdext and mpz_invert (including
    the choice of cofactors, so that results agree with GMP).
*/

void
flint_mpz_gcd(mpz_ptr g, mpz_srcptr u, mpz_srcptr v)
{
    mp_size_t usize, vsize;
    mp_ptr gp;

    usize = FLINT_ABS(u->_mp_size);
    vsize = FLINT_ABS(v->_mp_size);

    /* gcd(0, v) = |v| */
    if (usize == 0)
    {
        if (g != v)
        {
            gp = FLINT_MPZ_REALLOC(g, vsize);
            flint_mpn_copyi(gp, v->_mp_d, vsize);
        }
        g->_mp_size = vsize;
        return;
    }

    /* gcd(u, 0) = |u| */
    if (vsize == 0)
    {
        if (g != u)
        {
            gp = FLINT_MPZ_REALLOC(g, usize);
            flint_mpn_copyi(gp, u->_mp_d, usize);
        }
        g->_mp_size = usize;
        return;
    }

    /* if g is u or v, this does not reallocate, and flint_mpn_gcd allows
       its output to be one of the inputs */
    gp = FLINT_MPZ_REALLOC(g, FLINT_MIN(usize, vsize));
    g->_mp_size = flint_mpn_gcd(gp, u->_mp_d, usize, v->_mp_d, vsize);
}

/*
    g = gcd(a, b) = s a + t b with GMP's choice of cofactors. Any of g, s
    and t may be NULL (not computed); the outputs must be distinct but may
    alias the inputs.
*/
void
flint_mpz_gcdext(mpz_ptr g, mpz_ptr s, mpz_ptr t, mpz_srcptr a, mpz_srcptr b)
{
    mp_size_t asize, bsize;
    mp_ptr tmp_ap, tmp_bp;
    mp_size_t gsize, ssize, tmp_ssize;
    mp_ptr gp, tmp_gp, tmp_sp;
    TMP_INIT;

    /* _flint_mpn_gcdext_inplace requires that asize >= bsize. Therefore, we
       often have to swap a and b. The computed cofactor will be the
       "smallest" one, which is faster to produce. The wanted one will be
       computed here; this is needed anyway when both are requested. */
    asize = FLINT_ABS(a->_mp_size);
    bsize = FLINT_ABS(b->_mp_size);

    if (asize < bsize)
    {
        FLINT_SWAP(mpz_srcptr, a, b);
        FLINT_SWAP(mp_size_t, asize, bsize);
        FLINT_SWAP(mpz_ptr, s, t);
    }

    if (bsize == 0)
    {
        /* g = |a|, s = sgn(a), t = 0. */
        ssize = a->_mp_size >= 0 ? (asize != 0) : -1;

        if (g != NULL)
        {
            gp = FLINT_MPZ_REALLOC(g, asize);
            if (g != a)
                flint_mpn_copyi(gp, a->_mp_d, asize);
            g->_mp_size = asize;
        }
        if (t != NULL)
            t->_mp_size = 0;
        if (s != NULL)
        {
            FLINT_MPZ_REALLOC(s, 1)[0] = 1;
            s->_mp_size = ssize;
        }
        return;
    }

    TMP_START;

    tmp_gp = TMP_ALLOC((bsize + asize + bsize + bsize + 1) * sizeof(mp_limb_t));
    tmp_sp = tmp_gp + bsize;
    tmp_bp = tmp_sp + bsize + 1;
    tmp_ap = tmp_bp + bsize;

    flint_mpn_copyi(tmp_ap, a->_mp_d, asize);
    flint_mpn_copyi(tmp_bp, b->_mp_d, bsize);

    gsize = _flint_mpn_gcdext_inplace(tmp_gp, tmp_sp, &tmp_ssize, tmp_ap, asize, tmp_bp, bsize);

    ssize = FLINT_ABS(tmp_ssize);
    tmp_ssize = a->_mp_size >= 0 ? tmp_ssize : -tmp_ssize;

    if (t != NULL)
    {
        /* t = (g - s a) / b = sgn(b) (g - s' |a|) / |b|, where s' is the
           cofactor of |a| computed above (t may alias a or b, so it is
           computed in scratch space first) */
        mp_ptr qp = TMP_ALLOC(asize * sizeof(mp_limb_t));
        mp_size_t qn, sn1 = (a->_mp_size >= 0) ? tmp_ssize : -tmp_ssize;
        mp_ptr tp;

        _flint_mpn_gcdext_cofactor(qp, &qn, tmp_gp, gsize, tmp_sp, sn1,
            a->_mp_d, asize, b->_mp_d, bsize);

        tp = FLINT_MPZ_REALLOC(t, FLINT_ABS(qn));
        flint_mpn_copyi(tp, qp, FLINT_ABS(qn));
        t->_mp_size = (b->_mp_size >= 0) ? qn : -qn;
    }

    if (s != NULL)
    {
        mp_ptr sp = FLINT_MPZ_REALLOC(s, ssize);
        flint_mpn_copyi(sp, tmp_sp, ssize);
        s->_mp_size = tmp_ssize;
    }

    if (g != NULL)
    {
        gp = FLINT_MPZ_REALLOC(g, gsize);
        flint_mpn_copyi(gp, tmp_gp, gsize);
        g->_mp_size = gsize;
    }

    TMP_END;
}

/*
    r = x^(-1) mod |n| in [0, |n|) if gcd(x, n) = 1 (returning 1), else
    returns 0 with r undefined. n must be nonzero; r may alias x or n.

    GMP's mpz_invert calls mpz_gcdext(x, n), which for the usual |x| < |n|
    computes the cofactor of n and recovers that of x by a multiplication
    and an exact division. Since the inverse is unique, we instead reduce
    x modulo n and pass it first, zero-padded to the length of n, so that
    _flint_mpn_gcdext_inplace returns its cofactor directly.
*/
int
flint_mpz_invert(mpz_ptr r, mpz_srcptr x, mpz_srcptr n)
{
    mp_size_t xn, nn, gn, sn;
    mp_ptr xp, np, gp, sp, rp;
    int result;
    TMP_INIT;

    xn = FLINT_ABS(x->_mp_size);
    nn = FLINT_ABS(n->_mp_size);

    if (FLINT_UNLIKELY(nn == 0))
        flint_throw(FLINT_DIVZERO, "flint_mpz_invert: division by zero\n");

    TMP_START;

    xp = TMP_ALLOC((4 * nn + 3) * sizeof(mp_limb_t));
    np = xp + nn + 1;
    gp = np + nn + 1;
    sp = gp + nn;

    /* xp = |x| mod |n|, padded to nn limbs */
    if (xn >= nn)
    {
        mp_ptr qp = TMP_ALLOC((xn - nn + 1) * sizeof(mp_limb_t));
        flint_mpn_tdiv_qr(qp, xp, x->_mp_d, xn, n->_mp_d, nn);
    }
    else
    {
        flint_mpn_copyi(xp, x->_mp_d, xn);
        flint_mpn_zero(xp + xn, nn - xn);
    }

    if (flint_mpn_zero_p(xp, nn))
    {
        /* invertible only modulo 1, with inverse 0 */
        result = (nn == 1 && n->_mp_d[0] == 1);
        if (result)
            r->_mp_size = 0;
        TMP_END;
        return result;
    }

    flint_mpn_copyi(np, n->_mp_d, nn);

    gn = _flint_mpn_gcdext_inplace(gp, sp, &sn, xp, nn, np, nn);

    result = (gn == 1 && gp[0] == 1);

    if (result)
    {
        /* s |x| = 1 mod |n|, |s| < |n|; the inverse of x is s sgn(x) */
        if (x->_mp_size < 0)
            sn = -sn;

        if (sn == 0)
        {
            r->_mp_size = 0;
        }
        else if (sn > 0)
        {
            rp = FLINT_MPZ_REALLOC(r, sn);
            flint_mpn_copyi(rp, sp, sn);
            r->_mp_size = sn;
        }
        else
        {
            mp_size_t rn = nn;
            /* |n| - |s|, read before r (which may alias n) is written */
            flint_mpn_copyi(gp, n->_mp_d, nn);
            mpn_sub(gp, gp, nn, sp, -sn);
            MPN_NORM(gp, rn);
            rp = FLINT_MPZ_REALLOC(r, rn);
            flint_mpn_copyi(rp, gp, rn);
            r->_mp_size = rn;
        }
    }

    TMP_END;
    return result;
}
