/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"
#include "fmpz.h"
#include "gr.h"

void
radix_integer_init(radix_integer_t res, const radix_t radix)
{
    res->d = NULL;
    res->size = 0;
    res->alloc = 0;
}

void
radix_integer_clear(radix_integer_t res, const radix_t radix)
{
    if (res->d != NULL)
        flint_free(res->d);
}

nn_ptr
radix_integer_fit_limbs(radix_integer_t res, slong nlimbs, const radix_t radix)
{
    if (res->alloc < nlimbs)
    {
        /* todo: should we overallocate? */
        res->d = flint_realloc(res->d, nlimbs * sizeof(ulong));
        res->alloc = nlimbs;
    }

    return res->d;
}

void
radix_integer_zero(radix_integer_t res, const radix_t radix)
{
    /* todo: do we want to deallocate? */
    res->size = 0;
}

void
radix_integer_rand_limbs(radix_integer_t res, flint_rand_t state, slong n, const radix_t radix)
{
    nn_ptr d = radix_integer_fit_limbs(res, n, radix);
    radix_rand_limbs(d, state, n, radix);
    MPN_NORM(d, n);
    res->size = (n_randlimb(state) & 1) ? n : -n;
}

void
radix_integer_randtest_limbs(radix_integer_t res, flint_rand_t state, slong max_limbs, const radix_t radix)
{
    slong n = n_randint(state, max_limbs + 1);
    nn_ptr d = radix_integer_fit_limbs(res, n, radix);
    radix_randtest_limbs(d, state, n, radix);
    MPN_NORM(d, n);
    res->size = (n_randlimb(state) & 1) ? n : -n;
}

void
radix_integer_one(radix_integer_t res, const radix_t radix)
{
    radix_integer_fit_limbs(res, 1, radix)[0] = 1;
    res->size = 1;
}

void
radix_integer_neg_one(radix_integer_t res, const radix_t radix)
{
    radix_integer_fit_limbs(res, 1, radix)[0] = 1;
    res->size = -1;
}

int
radix_integer_is_zero(const radix_integer_t x, const radix_t radix)
{
    return x->size == 0;
}

int
radix_integer_is_one(const radix_integer_t x, const radix_t radix)
{
    return x->size == 1 && x->d[0] == 1;
}

int
radix_integer_is_neg_one(const radix_integer_t x, const radix_t radix)
{
    return x->size == -1 && x->d[0] == 1;
}

int
radix_integer_equal(const radix_integer_t x, const radix_integer_t y, const radix_t radix)
{
    return (x->size == y->size) && flint_mpn_equal_p(x->d, y->d, FLINT_ABS(x->size));
}

int
radix_integer_cmp(const radix_integer_t x, const radix_integer_t y, const radix_t radix)
{
    slong xsize = x->size;
    slong ysize = y->size;
    if (xsize != ysize)
        return (xsize > ysize) ? 1 : -1;
    int cmp = mpn_cmp(x->d, y->d, FLINT_ABS(xsize));
    return (xsize >= 0) ? FLINT_SGN(cmp) : -FLINT_SGN(cmp);
}

int
radix_integer_cmpabs(const radix_integer_t x, const radix_integer_t y, const radix_t radix)
{
    slong xsize = FLINT_ABS(x->size);
    slong ysize = FLINT_ABS(y->size);
    if (xsize != ysize)
        return (xsize > ysize) ? 1 : -1;
    int cmp = mpn_cmp(x->d, y->d, xsize);
    return FLINT_SGN(cmp);
}

void
radix_integer_set(radix_integer_t res, const radix_integer_t x, const radix_t radix)
{
    slong size = x->size;
    slong n = FLINT_ABS(size);

    if (n != 0)
    {
        radix_integer_fit_limbs(res, n, radix);
        flint_mpn_copyi(res->d, x->d, n);
    }

    res->size = size;
}

void
radix_integer_neg(radix_integer_t res, const radix_integer_t x, const radix_t radix)
{
    if (res != x)
        radix_integer_set(res, x, radix);

    res->size = -res->size;
}


void
radix_integer_abs(radix_integer_t res, const radix_integer_t x, const radix_t radix)
{
    if (res != x)
        radix_integer_set(res, x, radix);

    res->size = FLINT_ABS(res->size);
}

int
radix_integer_sgn(const radix_integer_t x, const radix_t radix)
{
    return FLINT_SGN(x->size);
}

void
radix_integer_set_ui(radix_integer_t res, ulong x, const radix_t radix)
{
    if (x == 0)
    {
        radix_integer_zero(res, radix);
    }
    else if (x < LIMB_RADIX(radix))
    {
        radix_integer_fit_limbs(res, 1, radix)[0] = x;
        res->size = 1;
    }
    else
    {
        res->size = radix_set_mpn(radix_integer_fit_limbs(res,
            radix_set_mpn_need_alloc(1, radix), radix), &x, 1, radix);
    }
}

void
radix_integer_set_si(radix_integer_t res, slong x, const radix_t radix)
{
    ulong ux = FLINT_UABS(x);

    if (x == 0)
    {
        radix_integer_zero(res, radix);
    }
    else if (ux < LIMB_RADIX(radix))
    {
        radix_integer_fit_limbs(res, 1, radix)[0] = ux;
        res->size = (x < 0) ? -1 : 1;
    }
    else
    {
        nn_ptr rd = radix_integer_fit_limbs(res, radix_set_mpn_need_alloc(1, radix), radix);
        slong rn = radix_set_mpn(rd, &ux, 1, radix);
        res->size = (x < 0) ? -rn : rn;
    }
}

void
radix_integer_set_fmpz(radix_integer_t res, const fmpz_t x, const radix_t radix)
{
    if (!COEFF_IS_MPZ(*x))
    {
        radix_integer_set_si(res, *x, radix);
    }
    else
    {
        slong xsize = COEFF_TO_PTR(*x)->_mp_size;
        nn_srcptr xd = COEFF_TO_PTR(*x)->_mp_d;
        nn_ptr rd = radix_integer_fit_limbs(res, radix_set_mpn_need_alloc(FLINT_ABS(xsize), radix), radix);
        slong rn = radix_set_mpn(rd, xd, FLINT_ABS(xsize), radix);
        res->size = (xsize > 0) ? rn : -rn;
    }
}

void
radix_integer_get_fmpz(fmpz_t res, const radix_integer_t x, const radix_t radix)
{
    slong xsize = x->size;

    if (xsize == 0)
        fmpz_zero(res);
    else if (xsize == 1)
        fmpz_set_ui(res, x->d[0]);
    else if (xsize == -1)
        fmpz_neg_ui(res, x->d[0]);
    else
    {
        __mpz_struct * r = _fmpz_promote(res);
        slong xn = FLINT_ABS(xsize);
        nn_srcptr xd = x->d;
        nn_ptr rd = FLINT_MPZ_REALLOC(r, xn);
        slong rn = radix_get_mpn(rd, xd, xn, radix);
        /* todo: check for possible int overflow */
        r->_mp_size = (xsize > 0) ? rn : -rn;
        _fmpz_demote_val(res);
    }
}

/* todo: 1x1 and nx1 specializations */
/* todo: maybe handle aliasing in radix_mul to simplify? */
void
radix_integer_mul(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, const radix_t radix)
{
    slong xsize = x->size;
    slong ysize = y->size;

    if (xsize == 0 || ysize == 0)
    {
        radix_integer_zero(res, radix);
        return;
    }

    slong xn = FLINT_ABS(xsize);
    slong yn = FLINT_ABS(ysize);
    slong rn = xn + yn;
    slong sgn = xsize ^ ysize;
    nn_ptr rd = radix_integer_fit_limbs(res, rn, radix);
    nn_srcptr xd = x->d;
    nn_srcptr yd = y->d;

    /* todo: squaring check first? */

    if (res == x)
    {
        TMP_INIT;
        TMP_START;
        nn_ptr tmp = TMP_ALLOC(sizeof(ulong) * xn);
        flint_mpn_copyi(tmp, xd, xn);
        if (x == y)
            radix_sqr(rd, tmp, xn, radix);
        else
            radix_mul(rd, tmp, xn, yd, yn, radix);
        TMP_END;
    }
    else if (res == y)
    {
        TMP_INIT;
        TMP_START;
        nn_ptr tmp = TMP_ALLOC(sizeof(ulong) * yn);
        flint_mpn_copyi(tmp, yd, yn);
        radix_mul(rd, xd, xn, tmp, yn, radix);
        TMP_END;
    }
    else
    {
        if (x == y)
            radix_sqr(rd, xd, xn, radix);
        else
            radix_mul(rd, xd, xn, yd, yn, radix);
    }

    rn -= (rd[rn - 1] == 0);
    res->size = (sgn >= 0) ? rn : -rn;
}

void
radix_integer_add(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, const radix_t radix)
{
    slong xsize = x->size;
    slong ysize = y->size;

    if (xsize == 0)
    {
        radix_integer_set(res, y, radix);
        return;
    }

    if (ysize == 0)
    {
        radix_integer_set(res, x, radix);
        return;
    }

    slong xn = FLINT_ABS(xsize);
    slong yn = FLINT_ABS(ysize);
    int xsgnbit = (xsize < 0);
    int ysgnbit = (ysize < 0);
    slong rn = FLINT_MAX(xn, yn) + (xsgnbit == ysgnbit);
    nn_ptr rd = radix_integer_fit_limbs(res, rn, radix);
    nn_srcptr xd = x->d;
    nn_srcptr yd = y->d;
    int rsgnbit = 0;
    ulong cy;

    if (xsgnbit == ysgnbit)
    {
        rsgnbit = xsgnbit;
        if (xn >= yn)
            rd[rn - 1] = cy = radix_add(rd, xd, xn, yd, yn, radix);
        else
            rd[rn - 1] = cy = radix_add(rd, yd, yn, xd, xn, radix);
        rn -= (cy == 0);
    }
    else
    {
        if (xn > yn || (xn == yn && mpn_cmp(xd, yd, xn) >= 0))
        {
            radix_sub(rd, xd, xn, yd, yn, radix);
            rsgnbit = xsgnbit;
        }
        else
        {
            radix_sub(rd, yd, yn, xd, xn, radix);
            rsgnbit = ysgnbit;
        }
        MPN_NORM(rd, rn);
    }

    res->size = rsgnbit ? -rn : rn;
}

void
radix_integer_sub(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, const radix_t radix)
{
    slong xsize = x->size;
    slong ysize = y->size;

    if (xsize == 0)
    {
        radix_integer_neg(res, y, radix);
        return;
    }

    if (ysize == 0)
    {
        radix_integer_set(res, x, radix);
        return;
    }

    slong xn = FLINT_ABS(xsize);
    slong yn = FLINT_ABS(ysize);
    int xsgnbit = (xsize < 0);
    int ysgnbit = (ysize > 0);
    slong rn = FLINT_MAX(xn, yn) + (xsgnbit == ysgnbit);
    nn_ptr rd = radix_integer_fit_limbs(res, rn, radix);
    nn_srcptr xd = x->d;
    nn_srcptr yd = y->d;
    int rsgnbit = 0;
    ulong cy;

    if (xsgnbit == ysgnbit)
    {
        rsgnbit = xsgnbit;
        if (xn >= yn)
            rd[rn - 1] = cy = radix_add(rd, xd, xn, yd, yn, radix);
        else
            rd[rn - 1] = cy = radix_add(rd, yd, yn, xd, xn, radix);
        rn -= (cy == 0);
    }
    else
    {
        if (xn > yn || (xn == yn && mpn_cmp(xd, yd, xn) >= 0))
        {
            radix_sub(rd, xd, xn, yd, yn, radix);
            rsgnbit = xsgnbit;
        }
        else
        {
            radix_sub(rd, yd, yn, xd, xn, radix);
            rsgnbit = ysgnbit;
        }
        MPN_NORM(rd, rn);
    }

    res->size = rsgnbit ? -rn : rn;
}

int
radix_integer_is_normalised(const radix_integer_t x, const radix_t radix)
{
    slong xn = FLINT_ABS(x->size);
    slong i;

    for (i = 0; i < xn; i++)
        if (x->d[i] >= LIMB_RADIX(radix))
            return 0;

    if (xn != 0 && x->d[xn - 1] == 0)
        return 0;

    return 1;
}

void
radix_integer_set_limb(radix_integer_t res, const radix_integer_t x, slong index, ulong c, const radix_t radix)
{
    FLINT_ASSERT(index >= 0);
    FLINT_ASSERT(c < LIMB_RADIX(radix));

    slong xsize = x->size;
    slong xn = FLINT_ABS(xsize);
    slong rn = xn;

    if (index < xn - 1 || (c != 0 && index == xn - 1))
    {
        radix_integer_set(res, x, radix);
        res->d[index] = c;
    }
    else if (c == 0)
    {
        if (index == xn - 1)
        {
            rn = xn - 1;
            MPN_NORM(x->d, rn);
            if (res != x && rn != 0)
                flint_mpn_copyi(radix_integer_fit_limbs(res, rn, radix), x->d, rn);
            res->size = (xsize >= 0) ? rn : -rn;
        }
        else
        {
            radix_integer_set(res, x, radix);
        }
    }
    else
    {
        rn = index + 1;
        nn_ptr rd = radix_integer_fit_limbs(res, rn, radix);
        if (res != x)
            flint_mpn_copyi(rd, x->d, xn);
        flint_mpn_zero(rd + xn, index - xn);
        rd[index] = c;
        res->size = (xsize >= 0) ? rn : -rn;
    }
}

void
radix_integer_lshift_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)
{
    FLINT_ASSERT(n >= 0);

    slong xsize = x->size;

    if (xsize == 0)
    {
        radix_integer_zero(res, radix);
        return;
    }

    slong xn = FLINT_ABS(xsize);
    slong rn = xn + n;

    nn_ptr rd = radix_integer_fit_limbs(res, rn, radix);
    flint_mpn_copyd(rd + n, x->d, xn);
    flint_mpn_zero(rd, n);
    res->size = (xsize > 0) ? rn : -rn;
}

void radix_integer_rshift_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)
{
    FLINT_ASSERT(n >= 0);

    slong xsize = x->size;
    slong xn = FLINT_ABS(xsize);
    slong rn = xn - n;

    if (rn <= 0)
    {
        radix_integer_zero(res, radix);
        return;
    }

    nn_ptr rd;
    if (res == x)
        rd = res->d;
    else
        rd = radix_integer_fit_limbs(res, rn, radix);

    flint_mpn_copyi(rd, x->d + n, rn);
    res->size = (xsize > 0) ? rn : -rn;
}

void
radix_integer_trunc_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)
{
    slong xsize = x->size;
    slong xn = FLINT_ABS(xsize);

    FLINT_ASSERT(n >= 0);

    if (res == x)
    {
        if (n < xn)
        {
            MPN_NORM(res->d, n);
            res->size = (xsize >= 0) ? n : -n;
        }
    }
    else
    {
        if (n < xn)
            MPN_NORM(x->d, n);
        else
            n = xn;

        flint_mpn_copyi(radix_integer_fit_limbs(res, n, radix), x->d, n);
        res->size = (xsize >= 0) ? n : -n;
    }
}

/* todo: optimize */
void
radix_integer_mod_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)
{
    radix_integer_trunc_limbs(res, x, n, radix);

    if (res->size < 0)
    {
        slong rn = FLINT_ABS(res->size);
        radix_neg(res->d, res->d, rn, radix);
        if (rn < n)
        {
            slong i;
            radix_integer_fit_limbs(res, n, radix);
            ulong B1 = LIMB_RADIX(radix) - 1;
            for (i = rn; i < n; i++)
                res->d[i] = B1;
        }
        MPN_NORM(res->d, n);
        res->size = n;
    }
}

void
radix_integer_smod_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)
{
    slong xsize = x->size;
    slong xn = FLINT_ABS(xsize);
    slong rn;

    FLINT_ASSERT(n >= 0);
    rn = FLINT_MIN(n, xn);
    MPN_NORM(x->d, rn);

    if (rn > 0 && rn == n)
    {
        int cmp = radix_cmp_bn_half(x->d, n, radix);

        if (cmp > 0 || (cmp == 0 && xsize > 0 && (LIMB_RADIX(radix) % 2 == 0)))
        {
            if (res == x)
                radix_neg(res->d, res->d, rn, radix);
            else
                radix_neg(radix_integer_fit_limbs(res, n, radix), x->d, rn, radix);

            MPN_NORM(res->d, rn);
            res->size = (xsize > 0) ? -rn : rn;
            return;
        }
    }

    if (res != x)
        flint_mpn_copyi(radix_integer_fit_limbs(res, rn, radix), x->d, rn);
    res->size = (xsize >= 0) ? rn : -rn;
}

void
radix_integer_mullow_limbs(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, slong n, const radix_t radix)
{
    slong xsize = x->size;
    slong ysize = y->size;

    if (xsize == 0 || ysize == 0 || n == 0)
    {
        radix_integer_zero(res, radix);
        return;
    }

    slong xn = FLINT_ABS(xsize);
    slong yn = FLINT_ABS(ysize);

    xn = FLINT_MIN(xn, n);
    yn = FLINT_MIN(yn, n);
    slong rn = FLINT_MIN(xn + yn, n);

    slong sgn = xsize ^ ysize;
    nn_ptr rd = radix_integer_fit_limbs(res, rn, radix);
    nn_srcptr xd = x->d;
    nn_srcptr yd = y->d;

    if (res == x)
    {
        TMP_INIT;
        TMP_START;
        nn_ptr tmp = TMP_ALLOC(sizeof(ulong) * xn);
        flint_mpn_copyi(tmp, xd, xn);
        if (x == y)
            radix_mulmid(rd, tmp, xn, tmp, yn, 0, rn, radix);
        else
            radix_mulmid(rd, tmp, xn, yd, yn, 0, rn, radix);
        TMP_END;
    }
    else if (res == y)
    {
        TMP_INIT;
        TMP_START;
        nn_ptr tmp = TMP_ALLOC(sizeof(ulong) * yn);
        flint_mpn_copyi(tmp, yd, yn);
        radix_mulmid(rd, xd, xn, tmp, yn, 0, rn, radix);
        TMP_END;
    }
    else
    {
        radix_mulmid(rd, xd, xn, yd, yn, 0, rn, radix);
    }

    MPN_NORM(rd, rn);
    res->size = (sgn >= 0) ? rn : -rn;
}

int
radix_integer_invmod_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)
{
    if (n == 0)
    {
        radix_integer_zero(res, radix);
        return 1;
    }

    slong xsize = x->size;
    if (xsize == 0)
        return 0;

    slong xn = FLINT_ABS(xsize);
    xn = FLINT_MIN(xn, n);

    if (xn == 1 && x->d[0] == 1)
    {
        /* preserve sign of x */
        radix_integer_trunc_limbs(res, x, 1, radix);
        return 1;
    }

    slong rn = n;
    int invertible;
    nn_ptr rd = radix_integer_fit_limbs(res, n, radix);
    nn_srcptr xd = x->d;

    if (res == x && rn != 1)
    {
        TMP_INIT;
        TMP_START;
        nn_ptr tmp = TMP_ALLOC(sizeof(ulong) * xn);
        flint_mpn_copyi(tmp, xd, xn);
        invertible = radix_invmod_bn(rd, tmp, xn, n, radix);
        TMP_END;
    }
    else
    {
        invertible = radix_invmod_bn(rd, xd, xn, n, radix);
    }

    MPN_NORM(rd, rn);
    res->size = (xsize > 0) ? rn : -rn;
    return invertible;
}

static ulong
radix_cdivrem(nn_ptr q, nn_ptr r, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)
{
    ulong cy = 0;

    if (b == r)
    {
        nn_ptr tmp;
        TMP_INIT;
        TMP_START;
        tmp = TMP_ALLOC(bn * sizeof(ulong));
        flint_mpn_copyi(tmp, b, bn);
        cy = radix_cdivrem(q, r, a, an, tmp, bn, radix);
        TMP_END;
    }
    else
    {
        radix_divrem(q, r, a, an, b, bn, radix);

        if (!flint_mpn_zero_p(r, bn))
        {
            ulong one = 1;
            radix_sub(r, b, bn, r, bn, radix);
            cy = radix_add(q, q, an - bn + 1, &one, 1, radix);
        }
    }

    return cy;
}

int radix_integer_div(radix_integer_t q,
    const radix_integer_t a, const radix_integer_t b, const radix_t radix)
{
    slong asize, bsize, an, bn, qn;
    nn_srcptr ad, bd;
    nn_ptr qd;
    int exact;

    asize = a->size;
    bsize = b->size;

    an = FLINT_ABS(asize);
    bn = FLINT_ABS(bsize);

    if (bn == 0)
        return 0;

    if (an == 0)
    {
        radix_integer_zero(q, radix);
        return 1;
    }

    qn = an - bn + 1;
    ad = a->d;
    bd = b->d;

    if (an < bn || (an == bn && mpn_cmp(ad, bd, an) < 0))
    {
        radix_integer_zero(q, radix);
        return 0;
    }

    qd = radix_integer_fit_limbs(q, qn, radix);

    exact = radix_div(qd, ad, an, bd, bn, radix);
    if (!exact)
        qn = 0;
    MPN_NORM(qd, qn);
    q->size = ((asize < 0) == (bsize < 0)) ? qn : -qn;

    return exact;
}

void radix_integer_divexact(radix_integer_t q,
    const radix_integer_t a, const radix_integer_t b, const radix_t radix)
{
    slong asize, bsize, an, bn, qn;
    nn_srcptr ad, bd;
    nn_ptr qd;

    asize = a->size;
    bsize = b->size;

    an = FLINT_ABS(asize);
    bn = FLINT_ABS(bsize);

    if (an == 0)
    {
        radix_integer_zero(q, radix);
        return;
    }

    if (bn == 0 || an < bn)
        flint_throw(FLINT_ERROR, "radix_integer_divexact: an < bn");

    qn = an - bn + 1;
    ad = a->d;
    bd = b->d;

    qd = radix_integer_fit_limbs(q, qn, radix);

    radix_divexact(qd, ad, an, bd, bn, radix);
    MPN_NORM(qd, qn);
    q->size = ((asize < 0) == (bsize < 0)) ? qn : -qn;
}

void radix_integer_tdiv_qr(radix_integer_t q, radix_integer_t r,
    const radix_integer_t a, const radix_integer_t b, const radix_t radix)
{
    slong asize, bsize, an, bn, qn, rn;
    nn_srcptr ad, bd;
    nn_ptr qd, rd;

    asize = a->size;
    bsize = b->size;

    an = FLINT_ABS(asize);
    bn = FLINT_ABS(bsize);

    if (bn == 0)
        flint_throw(FLINT_DIVZERO, "radix_integer_tdiv_qr: divide by zero");

    if (an == 0)
    {
        radix_integer_zero(q, radix);
        radix_integer_zero(r, radix);
        return;
    }

    if (an < bn)
    {
        radix_integer_set(r, a, radix);
        radix_integer_zero(q, radix);
        return;
    }

    qn = an - bn + 1;
    rn = bn;

    ad = a->d;
    bd = b->d;

    qd = radix_integer_fit_limbs(q, qn, radix);
    rd = radix_integer_fit_limbs(r, rn, radix);

    radix_divrem(qd, rd, ad, an, bd, bn, radix);

    MPN_NORM(qd, qn);
    MPN_NORM(rd, rn);

    q->size = ((asize < 0) == (bsize < 0)) ? qn : -qn;
    r->size = (asize < 0) ? -rn : rn;
}

void radix_integer_fdiv_qr(radix_integer_t q, radix_integer_t r,
    const radix_integer_t a, const radix_integer_t b, const radix_t radix)
{
    slong asize, bsize, an, bn, qn, rn;
    nn_srcptr ad, bd;
    nn_ptr qd, rd;

    asize = a->size;
    bsize = b->size;

    an = FLINT_ABS(asize);
    bn = FLINT_ABS(bsize);

    if (bn == 0)
        flint_throw(FLINT_DIVZERO, "radix_integer_fdiv_qr: divide by zero");

    if (an == 0)
    {
        radix_integer_zero(q, radix);
        radix_integer_zero(r, radix);
        return;
    }

    if (an < bn)
    {
        if ((asize < 0) == (bsize < 0))
        {
            radix_integer_set(r, a, radix);
            radix_integer_zero(q, radix);
        }
        else
        {
            radix_integer_add(r, a, b, radix);
            radix_integer_neg_one(q, radix);
        }
        return;
    }

    qn = an - bn + 1;
    rn = bn;

    ad = a->d;
    bd = b->d;

    qd = radix_integer_fit_limbs(q, qn, radix);
    rd = radix_integer_fit_limbs(r, rn, radix);

    if ((asize < 0) == (bsize < 0))
    {
        radix_divrem(qd, rd, ad, an, bd, bn, radix);
    }
    else
    {
        ulong cy = radix_cdivrem(qd, rd, ad, an, bd, bn, radix);

        if (cy != 0)
        {
            qd = radix_integer_fit_limbs(q, qn + 1, radix);
            qd[qn] = 1;
            qn++;
        }
    }

    MPN_NORM(qd, qn);
    MPN_NORM(rd, rn);

    q->size = ((asize < 0) == (bsize < 0)) ? qn : -qn;
    r->size = (bsize < 0) ? -rn : rn;
}

void radix_integer_cdiv_qr(radix_integer_t q, radix_integer_t r,
    const radix_integer_t a, const radix_integer_t b, const radix_t radix)
{
    slong asize, bsize, an, bn, qn, rn;
    nn_srcptr ad, bd;
    nn_ptr qd, rd;

    asize = a->size;
    bsize = b->size;

    an = FLINT_ABS(asize);
    bn = FLINT_ABS(bsize);

    if (bn == 0)
        flint_throw(FLINT_DIVZERO, "radix_integer_cdiv_qr: divide by zero");

    if (an == 0)
    {
        radix_integer_zero(q, radix);
        radix_integer_zero(r, radix);
        return;
    }

    if (an < bn)
    {
        if ((asize < 0) != (bsize < 0))
        {
            radix_integer_set(r, a, radix);
            radix_integer_zero(q, radix);
        }
        else
        {
            radix_integer_sub(r, a, b, radix);
            radix_integer_one(q, radix);
        }
        return;
    }

    qn = an - bn + 1;
    rn = bn;

    ad = a->d;
    bd = b->d;

    qd = radix_integer_fit_limbs(q, qn, radix);
    rd = radix_integer_fit_limbs(r, rn, radix);

    if ((asize < 0) != (bsize < 0))
    {
        radix_divrem(qd, rd, ad, an, bd, bn, radix);
    }
    else
    {
        ulong cy = radix_cdivrem(qd, rd, ad, an, bd, bn, radix);

        if (cy != 0)
        {
            qd = radix_integer_fit_limbs(q, qn + 1, radix);
            qd[qn] = 1;
            qn++;
        }
    }

    MPN_NORM(qd, qn);
    MPN_NORM(rd, rn);

    q->size = ((asize < 0) == (bsize < 0)) ? qn : -qn;
    r->size = (bsize < 0) ? rn : -rn;
}

/* ------------------------------------------------------------------------- */
/*    GR wrapper                                                             */
/* ------------------------------------------------------------------------- */

#define GR_RADIX_CTX(ctx) ((radix_struct *) (GR_CTX_DATA_AS_PTR(ctx)))

static void
_gr_radix_integer_ctx_clear(gr_ctx_t ctx)
{
    radix_clear(GR_RADIX_CTX(ctx));
    flint_free(GR_RADIX_CTX(ctx));
}

static int
_gr_radix_integer_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_stream_write(out, "Integers in radix ");
    status |= gr_stream_write_ui(out, DIGIT_RADIX(GR_RADIX_CTX(ctx)));
    status |= gr_stream_write(out, "^");
    status |= gr_stream_write_ui(out, GR_RADIX_CTX(ctx)->exp);
    status |= gr_stream_write(out, " (radix_integer)");
    return status;
}

static void
_gr_radix_integer_init(radix_integer_t res, gr_ctx_t ctx)
{
    radix_integer_init(res, GR_RADIX_CTX(ctx));
}

static void
_gr_radix_integer_clear(radix_integer_t res, gr_ctx_t ctx)
{
    radix_integer_clear(res, GR_RADIX_CTX(ctx));
}

static void _gr_radix_integer_swap(radix_integer_t x, radix_integer_t y, gr_ctx_t ctx)
{
    FLINT_SWAP(radix_integer_struct, *x, *y);
}

static void _gr_radix_integer_set_shallow(radix_integer_t res, const radix_integer_t x, gr_ctx_t ctx)
{
    *res = *x;
}

static int _gr_radix_integer_set(radix_integer_t res, const radix_integer_t x, gr_ctx_t ctx)
{
    radix_integer_set(res, x, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_radix_integer_randtest(radix_integer_t res, flint_rand_t state, gr_ctx_t ctx)
{
    radix_integer_randtest_limbs(res, state, (n_randlimb(state) & 7) ? 1 : 3, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_radix_integer_write(gr_stream_t out, const radix_integer_t x, gr_ctx_t ctx)
{
    slong size = x->size;
    int status = GR_SUCCESS;

    if (size == 0)
    {
        status |= gr_stream_write(out, "0");
    }
    else
    {
        if (DIGIT_RADIX(GR_RADIX_CTX(ctx)) == 10)
        {
            /* todo: chunked output for huge strings */
            status |= gr_stream_write_free(out, radix_get_str_decimal(NULL,
                x->d, FLINT_ABS(size), size < 0, GR_RADIX_CTX(ctx)));
        }
        else
        {
            /* todo: chunked output for huge strings */
            status |= gr_stream_write_free(out, radix_get_str_sum(NULL,
                x->d, FLINT_ABS(size), size < 0, 0, GR_RADIX_CTX(ctx)));
        }
    }

    return status;
}

static int
_gr_radix_integer_zero(radix_integer_t res, gr_ctx_t ctx)
{
    radix_integer_zero(res, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_radix_integer_one(radix_integer_t res, gr_ctx_t ctx)
{
    radix_integer_one(res, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_radix_integer_neg_one(radix_integer_t res, gr_ctx_t ctx)
{
    radix_integer_neg_one(res, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int _gr_radix_integer_set_ui(radix_integer_t res, ulong x, gr_ctx_t ctx)
{
    radix_integer_set_ui(res, x, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int _gr_radix_integer_set_si(radix_integer_t res, slong x, gr_ctx_t ctx)
{
    radix_integer_set_si(res, x, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int _gr_radix_integer_set_fmpz(radix_integer_t res, const fmpz_t x, gr_ctx_t ctx)
{
    radix_integer_set_fmpz(res, x, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int _gr_radix_integer_get_fmpz(fmpz_t res, const radix_integer_t x, gr_ctx_t ctx)
{
    radix_integer_get_fmpz(res, x, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int _gr_radix_integer_get_ui(ulong * res, const radix_integer_t x, gr_ctx_t ctx)
{
    if (x->size < 0)
        return GR_DOMAIN;

    if (x->size <= 1)
    {
        *res = x->size ? x->d[0] : 0;
        return GR_SUCCESS;
    }
    else
    {
        /* todo */
        fmpz_t t;
        fmpz_init(t);
        radix_integer_get_fmpz(t, x, GR_RADIX_CTX(ctx));
        if (fmpz_cmp_ui(t, UWORD_MAX) > 0)
        {
            fmpz_clear(t);
            *res = 0;
            return GR_DOMAIN;
        }
        else
        {
            *res = fmpz_get_ui(t);
            fmpz_clear(t);
            return GR_SUCCESS;
        }
    }
}

static int _gr_radix_integer_get_si(ulong * res, const radix_integer_t x, gr_ctx_t ctx)
{
    if (x->size == 0)
    {
        *res = 0;
        return GR_SUCCESS;
    }
    else
    {
        /* todo */
        fmpz_t t;
        fmpz_init(t);
        radix_integer_get_fmpz(t, x, GR_RADIX_CTX(ctx));
        if (!fmpz_fits_si(t))
        {
            fmpz_clear(t);
            *res = 0;
            return GR_DOMAIN;
        }
        else
        {
            *res = fmpz_get_si(t);
            fmpz_clear(t);
            return GR_SUCCESS;
        }
    }
}

static truth_t _gr_radix_integer_equal(const radix_integer_t x, const radix_integer_t y, gr_ctx_t ctx)
{
    return radix_integer_equal(x, y, GR_RADIX_CTX(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t _gr_radix_integer_is_zero(const radix_integer_t x, gr_ctx_t ctx)
{
    return radix_integer_is_zero(x, GR_RADIX_CTX(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t _gr_radix_integer_is_one(const radix_integer_t x, gr_ctx_t ctx)
{
    return radix_integer_is_one(x, GR_RADIX_CTX(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t _gr_radix_integer_is_neg_one(const radix_integer_t x, gr_ctx_t ctx)
{
    return radix_integer_is_neg_one(x, GR_RADIX_CTX(ctx)) ? T_TRUE : T_FALSE;
}

static int _gr_radix_integer_neg(radix_integer_t res, const radix_integer_t x, gr_ctx_t ctx)
{
    radix_integer_neg(res, x, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int _gr_radix_integer_abs(radix_integer_t res, const radix_integer_t x, gr_ctx_t ctx)
{
    radix_integer_abs(res, x, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int _gr_radix_integer_add(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, gr_ctx_t ctx)
{
    radix_integer_add(res, x, y, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int _gr_radix_integer_sub(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, gr_ctx_t ctx)
{
    radix_integer_sub(res, x, y, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int _gr_radix_integer_mul(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, gr_ctx_t ctx)
{
    radix_integer_mul(res, x, y, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int _gr_radix_integer_div(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, gr_ctx_t ctx)
{
    return radix_integer_div(res, x, y, GR_RADIX_CTX(ctx)) ? GR_SUCCESS : GR_DOMAIN;
}

static int _gr_radix_integer_divexact(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, gr_ctx_t ctx)
{
    if (y->size == 0)
        return GR_DOMAIN;

    radix_integer_divexact(res, x, y, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int _gr_radix_integer_cmp(int * res, const radix_integer_t x, const radix_integer_t y, gr_ctx_t ctx)
{
    *res = radix_integer_cmp(x, y, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int _gr_radix_integer_cmpabs(int * res, const radix_integer_t x, const radix_integer_t y, gr_ctx_t ctx)
{
    *res = radix_integer_cmpabs(x, y, GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}

static int _gr_radix_integer_sgn(radix_integer_t res, const radix_integer_t x, gr_ctx_t ctx)
{
    radix_integer_set_si(res, radix_integer_sgn(x, GR_RADIX_CTX(ctx)), GR_RADIX_CTX(ctx));
    return GR_SUCCESS;
}




int _gr_radix_integer_methods_initialized = 0;


gr_static_method_table _gr_radix_integer_methods;

gr_method_tab_input _gr_radix_integer_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_radix_integer_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_radix_integer_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_THREADSAFE,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_radix_integer_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_radix_integer_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_radix_integer_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_radix_integer_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_radix_integer_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_radix_integer_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_radix_integer_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_radix_integer_one},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) _gr_radix_integer_neg_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_radix_integer_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_radix_integer_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_radix_integer_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_radix_integer_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_radix_integer_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_radix_integer_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_radix_integer_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_radix_integer_set_fmpz},
    {GR_METHOD_GET_SI,          (gr_funcptr) _gr_radix_integer_get_si},
    {GR_METHOD_GET_UI,          (gr_funcptr) _gr_radix_integer_get_ui},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_radix_integer_get_fmpz},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_radix_integer_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_radix_integer_add},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_radix_integer_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_radix_integer_mul},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_radix_integer_div},
    {GR_METHOD_DIVEXACT,        (gr_funcptr) _gr_radix_integer_divexact},
    {GR_METHOD_CMP,             (gr_funcptr) _gr_radix_integer_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) _gr_radix_integer_cmpabs},
    {GR_METHOD_SGN,             (gr_funcptr) _gr_radix_integer_sgn},
    {GR_METHOD_ABS,             (gr_funcptr) _gr_radix_integer_abs},
    {GR_METHOD_FLOOR,           (gr_funcptr) _gr_radix_integer_set},
    {GR_METHOD_CEIL,            (gr_funcptr) _gr_radix_integer_set},
    {GR_METHOD_TRUNC,           (gr_funcptr) _gr_radix_integer_set},
    {GR_METHOD_NINT,            (gr_funcptr) _gr_radix_integer_set},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_radix_integer(gr_ctx_t ctx, ulong n, unsigned int exp)
{
    ctx->which_ring = GR_CTX_RADIX_INTEGER;
    ctx->sizeof_elem = sizeof(radix_integer_struct);
    ctx->size_limit = WORD_MAX;

    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(radix_struct));
    /* todo: error handling */
    radix_init(GR_RADIX_CTX(ctx), n, exp);

    ctx->methods = _gr_radix_integer_methods;

    if (!_gr_radix_integer_methods_initialized)
    {
        gr_method_tab_init(_gr_radix_integer_methods, _gr_radix_integer_methods_input);
        _gr_radix_integer_methods_initialized = 1;
    }
}
