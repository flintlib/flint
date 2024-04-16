/*
    Copyright (C) 2020 William Hart
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "mpn_mod.h"

/* todo: separate method tables, more fixed-width optimizations */
/* todo: fast polynomial multiplication */
/* todo: powering */
/* todo: special functions */

void flint_mpn_mulmod_preinvn_2(mp_ptr r,
        mp_srcptr a, mp_srcptr b,
        mp_srcptr d, mp_srcptr dinv, ulong norm)
{
    mp_limb_t cy, p1, p2, b0, b1, r0, r1;
    mp_limb_t t[10];

    if (norm)
    {
        /* mpn_lshift(b, b, n, norm) */
        b0 = (b[0] << norm);
        b1 = (b[1] << norm) | (b[0] >> (FLINT_BITS - norm));
    }
    else
    {
        b0 = b[0];
        b1 = b[1];
    }

    /* mpn_mul_n(t, a, b, n) */
    umul_ppmm(t[1], t[0], a[0], b0);
    umul_ppmm(t[3], t[2], a[1], b1);
    umul_ppmm(p2, p1, a[0], b1);
    add_sssaaaaaa(t[3], t[2], t[1], t[3], t[2], t[1], 0, p2, p1);
    umul_ppmm(p2, p1, a[1], b0);
    add_sssaaaaaa(t[3], t[2], t[1], t[3], t[2], t[1], 0, p2, p1);

    /* mpn_mul_n(t + 3*n, t + n, dinv, n) */
    umul_ppmm(t[7], t[6], t[2], dinv[0]);
    umul_ppmm(t[9], t[8], t[3], dinv[1]);
    umul_ppmm(p2, p1, t[2], dinv[1]);
    add_sssaaaaaa(t[9], t[8], t[7], t[9], t[8], t[7], 0, p2, p1);
    umul_ppmm(p2, p1, t[3], dinv[0]);
    add_sssaaaaaa(t[9], t[8], t[7], t[9], t[8], t[7], 0, p2, p1);

    /* mpn_add_n(t + 4*n, t + 4*n, t + n, n) */
    add_ssaaaa(t[9], t[8], t[9], t[8], t[3], t[2]);

    /* mpn_mul_n(t + 2*n, t + 4*n, d, n) */
    umul_ppmm(t[5], t[4], t[8], d[0]);
    t[6] = t[9]*d[1];
    umul_ppmm(p2, p1, t[8], d[1]);
    add_ssaaaa(t[6], t[5], t[6], t[5], p2, p1);
    umul_ppmm(p2, p1, t[9], d[0]);
    add_ssaaaa(t[6], t[5], t[6], t[5], p2, p1);

    /* cy = t[n] - t[3*n] - mpn_sub_n(r, t, t + 2*n, n) */
    sub_dddmmmsss(cy, r1, r0, t[2], t[1], t[0], t[6], t[5], t[4]);

    while (cy > 0)
    {
        /* cy -= mpn_sub_n(r, r, d, n) */
        sub_dddmmmsss(cy, r1, r0, cy, r1, r0, 0, d[1], d[0]);
    }

    if ((r1 > d[1]) || (r1 == d[1] && r0 >= d[0]))
    {
        /* mpn_sub_n(r, r, d, n) */
        sub_ddmmss(r1, r0, r1, r0, d[1], d[0]);
    }

    if (norm)
    {
        r[0] = (r0 >> norm) | (r1 << (FLINT_BITS - norm));
        r[1] = (r1 >> norm);
    }
    else
    {
        r[0] = r0;
        r[1] = r1;
    }
}

/* todo: efficient implementation */
char *
_flint_mpn_get_str(mp_srcptr x, mp_size_t n)
{
    fmpz_t t;
    char * s;
    fmpz_init(t);
    fmpz_set_ui_array(t, x, n);
    s = fmpz_get_str(NULL, 10, t);
    fmpz_clear(t);
    return s;
}

int
mpn_mod_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Integers mod ");
    gr_stream_write_free(out, _flint_mpn_get_str(MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx)));
    gr_stream_write(out, " (mpn)");
    return GR_SUCCESS;
}

void
mpn_mod_ctx_clear(gr_ctx_t ctx)
{
    flint_free(MPN_MOD_CTX(ctx));
}

int
mpn_mod_set_ui(mp_ptr res, ulong x, gr_ctx_t ctx)
{
    FLINT_ASSERT(MPN_MOD_CTX_NLIMBS(ctx) >= 2);

    res[0] = x;
    flint_mpn_zero(res + 1, MPN_MOD_CTX_NLIMBS(ctx) - 1);
    return GR_SUCCESS;
}

int
mpn_mod_set_si(mp_ptr res, slong x, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    FLINT_ASSERT(n >= 2);

    if (x >= 0)
    {
        res[0] = x;
        flint_mpn_zero(res + 1, n - 1);
    }
    else
    {
        mpn_sub_1(res, MPN_MOD_CTX_MODULUS(ctx), n, -(ulong) x);
    }

    return GR_SUCCESS;
}

int
mpn_mod_neg_one(mp_ptr res, gr_ctx_t ctx)
{
    return mpn_mod_set_si(res, -1, ctx);
}

/* fixme: flint_mpn_mod_preinvn is misdocumented */

int
mpn_mod_set_mpn(mp_ptr res, mp_srcptr x, mp_size_t xn, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
    mp_size_t rn;

    if (xn < n || (xn == n && mpn_cmp(x, MPN_MOD_CTX_MODULUS(ctx), n) < 0))
    {
        flint_mpn_copyi(res, x, xn);
        flint_mpn_zero(res + xn, n - xn);
    }
    else
    {
        mp_ptr r;
        TMP_INIT;

        TMP_START;
        r = TMP_ALLOC((xn + 1) * sizeof(mp_limb_t));

        if (norm == 0)
        {
            flint_mpn_copyi(r, x, xn);
            rn = xn;
        }
        else
        {
            r[xn] = mpn_lshift(r, x, xn, norm);
            rn = xn + (r[xn] != 0);
        }

        /* todo: good preinv code for nx2 reduction */
        if (n == 2)
        {
            mpn_divrem_2(r + 2, 0, r, rn, MPN_MOD_CTX_MODULUS_NORMED(ctx));

            if (norm == 0)
            {
                res[0] = r[0];
                res[1] = r[1];
            }
            else
            {
                res[0] = (r[0] >> norm) | (r[1] << (FLINT_BITS - norm));
                res[1] = (r[1] >> norm);
            }
        }
        else
        {
            if (norm == 0)
            {
                flint_mpn_mod_preinvn(r, x, rn, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx), MPN_MOD_CTX_MODULUS_PREINV(ctx));
                flint_mpn_copyi(res, r, n);
            }
            else
            {
                flint_mpn_mod_preinvn(r, r, rn, MPN_MOD_CTX_MODULUS_NORMED(ctx), MPN_MOD_CTX_NLIMBS(ctx), MPN_MOD_CTX_MODULUS_PREINV(ctx));
                mpn_rshift(res, r, n, norm);
            }
        }

        TMP_END;
    }

    return GR_SUCCESS;
}

int
mpn_mod_set_fmpz(mp_ptr res, const fmpz_t x, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*x))
        mpn_mod_set_si(res, *x, ctx);
    else
    {
        mp_size_t nd = FLINT_ABS(COEFF_TO_PTR(*x)->_mp_size);
        int neg = COEFF_TO_PTR(*x)->_mp_size < 0;

        mpn_mod_set_mpn(res, COEFF_TO_PTR(*x)->_mp_d, nd, ctx);
        if (neg)
            flint_mpn_negmod_n(res, res, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx));
    }
    return GR_SUCCESS;
}

/* todo: public interface */
typedef struct
{
    fmpz_mod_ctx_struct * ctx;
    truth_t is_prime;
    fmpz a;    /* when used as finite field with defining polynomial x - a */
}
_gr_fmpz_mod_ctx_struct;

#define FMPZ_MOD_CTX(ring_ctx) ((((_gr_fmpz_mod_ctx_struct *)(ring_ctx))->ctx))

int
mpn_mod_set_other(mp_ptr res, gr_ptr v, gr_ctx_t v_ctx, gr_ctx_t ctx)
{
    if (v_ctx == ctx)
        return mpn_mod_set(res, v, ctx);

    if (v_ctx->which_ring == GR_CTX_MPN_MOD)
    {
        mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);

        if (MPN_MOD_CTX_NLIMBS(v_ctx) == n &&
            flint_mpn_equal_p(MPN_MOD_CTX_MODULUS(v_ctx), MPN_MOD_CTX_MODULUS(ctx), n))
        {
            flint_mpn_copyi(res, v, n);
            return GR_SUCCESS;
        }
    }

    if (v_ctx->which_ring == GR_CTX_FMPZ_MOD)
    {
        mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);

        if (fmpz_size(FMPZ_MOD_CTX(v_ctx)->n) == n)
        {
            mp_srcptr vd = COEFF_TO_PTR(*(FMPZ_MOD_CTX(v_ctx)->n))->_mp_d;

            if (flint_mpn_equal_p(vd, MPN_MOD_CTX_MODULUS(ctx), n))
                return mpn_mod_set_fmpz(res, v, ctx);
        }
    }

    return gr_generic_set_other(res, v, v_ctx, ctx);
}

int
mpn_mod_randtest(mp_ptr res, flint_rand_t state, gr_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_ui_array(t, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx));
    fmpz_randtest_mod(t, state, t);
    GR_IGNORE(mpn_mod_set_fmpz(res, t, ctx));
    fmpz_clear(t);
    return GR_SUCCESS;
}

int
mpn_mod_write(gr_stream_t out, mp_srcptr x, gr_ctx_t ctx)
{
    gr_stream_write_free(out, _flint_mpn_get_str(x, MPN_MOD_CTX_NLIMBS(ctx)));
    return GR_SUCCESS;
}

int
mpn_mod_get_fmpz(fmpz_t res, mp_srcptr x, gr_ctx_t ctx)
{
    fmpz_set_ui_array(res, x, MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}


truth_t
mpn_mod_is_neg_one(gr_srcptr x, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    mpn_mod_neg_one(t, ctx);
    return flint_mpn_equal_p(x, t, MPN_MOD_CTX_NLIMBS(ctx)) ? T_TRUE : T_FALSE;
}

int
mpn_mod_neg(mp_ptr res, mp_srcptr x, gr_ctx_t ctx)
{
    if (MPN_MOD_CTX_NLIMBS(ctx) == 2)
        flint_mpn_negmod_2(res, x, MPN_MOD_CTX_MODULUS(ctx));
    else
        flint_mpn_negmod_n(res, x, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

int
mpn_mod_add(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx)
{
    if (MPN_MOD_CTX_NLIMBS(ctx) == 2)
        flint_mpn_addmod_2(res, x, y, MPN_MOD_CTX_MODULUS(ctx));
    else
        flint_mpn_addmod_n(res, x, y, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

int
mpn_mod_sub(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx)
{
    if (MPN_MOD_CTX_NLIMBS(ctx) == 2)
        flint_mpn_submod_2(res, x, y, MPN_MOD_CTX_MODULUS(ctx));
    else
        flint_mpn_submod_n(res, x, y, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

int
mpn_mod_add_ui(mp_ptr res, mp_srcptr x, ulong y, gr_ctx_t ctx)
{
    if (MPN_MOD_CTX_NLIMBS(ctx) == 2)
    {
        mp_limb_t t[2];
        t[0] = y;
        t[1] = 0;
        flint_mpn_addmod_2(res, x, t, MPN_MOD_CTX_MODULUS(ctx));
    }
    else
    {
        flint_mpn_addmod_n_m(res, x, &y, 1, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx));
    }
    return GR_SUCCESS;
}

int
mpn_mod_sub_ui(mp_ptr res, mp_srcptr x, ulong y, gr_ctx_t ctx)
{
    if (MPN_MOD_CTX_NLIMBS(ctx) == 2)
    {
        mp_limb_t t[2];
        t[0] = y;
        t[1] = 0;
        flint_mpn_submod_2(res, x, t, MPN_MOD_CTX_MODULUS(ctx));
    }
    else
    {
        flint_mpn_submod_n_m(res, x, &y, 1, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx));
    }
    return GR_SUCCESS;
}

int
mpn_mod_add_si(mp_ptr res, mp_srcptr x, slong y, gr_ctx_t ctx)
{
    if (y >= 0)
        mpn_mod_add_ui(res, x, (ulong) y, ctx);
    else
        mpn_mod_sub_ui(res, x, -(ulong) y, ctx);
    return GR_SUCCESS;
}

int
mpn_mod_sub_si(mp_ptr res, mp_srcptr x, slong y, gr_ctx_t ctx)
{
    if (y >= 0)
        mpn_mod_sub_ui(res, x, (ulong) y, ctx);
    else
        mpn_mod_add_ui(res, x, -(ulong) y, ctx);
    return GR_SUCCESS;
}

int
mpn_mod_add_fmpz(mp_ptr res, mp_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*y))
    {
        mpn_mod_add_si(res, x, *y, ctx);
    }
    else
    {
        mp_limb_t t[MPN_MOD_MAX_LIMBS];
        mp_srcptr m = MPN_MOD_CTX_MODULUS(ctx);
        mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);

        __mpz_struct * z = COEFF_TO_PTR(*y);
        mp_size_t ssize = z->_mp_size;
        mp_size_t zn = FLINT_ABS(ssize);
        mp_srcptr zd = z->_mp_d;

        if (zn < n || (zn == n && mpn_cmp(zd, m, n) < 0))
        {
            if (ssize > 0)
                flint_mpn_addmod_n_m(res, x, zd, zn, m, n);
            else
                flint_mpn_submod_n_m(res, x, zd, zn, m, n);
        }
        else
        {
            mpn_mod_set_mpn(t, zd, zn, ctx);
            if (ssize > 0)
                mpn_mod_add(res, x, t, ctx);
            else
                mpn_mod_sub(res, x, t, ctx);
        }
    }

    return GR_SUCCESS;
}

int
mpn_mod_sub_fmpz(mp_ptr res, mp_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*y))
    {
        mpn_mod_sub_si(res, x, *y, ctx);
    }
    else
    {
        mp_limb_t t[MPN_MOD_MAX_LIMBS];
        mp_srcptr m = MPN_MOD_CTX_MODULUS(ctx);
        mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);

        __mpz_struct * z = COEFF_TO_PTR(*y);
        mp_size_t ssize = z->_mp_size;
        mp_size_t zn = FLINT_ABS(ssize);
        mp_srcptr zd = z->_mp_d;

        if (zn < n || (zn == n && mpn_cmp(zd, m, n) < 0))
        {
            if (ssize > 0)
                flint_mpn_submod_n_m(res, x, zd, zn, m, n);
            else
                flint_mpn_addmod_n_m(res, x, zd, zn, m, n);
        }
        else
        {
            mpn_mod_set_mpn(t, zd, zn, ctx);
            if (ssize > 0)
                mpn_mod_sub(res, x, t, ctx);
            else
                mpn_mod_add(res, x, t, ctx);
        }
    }

    return GR_SUCCESS;
}


/* should mirror flint_mpn_mulmod_preinvn (backport any improvements) */
int
mpn_mod_mul(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);

    if (n == 2)
    {
        flint_mpn_mulmod_preinvn_2(res, x, y, MPN_MOD_CTX_MODULUS_NORMED(ctx), MPN_MOD_CTX_MODULUS_PREINV(ctx), MPN_MOD_CTX_NORM(ctx));
    }
    else
    {
        mp_limb_t t[5 * MPN_MOD_MAX_LIMBS];
        mp_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
        mp_ptr d = MPN_MOD_CTX_MODULUS_NORMED(ctx);
        mp_ptr dinv = MPN_MOD_CTX_MODULUS_PREINV(ctx);
        mp_limb_t cy;

        if (norm)
        {
            if (x == y)
            {
                flint_mpn_sqr(t, x, n);
                mpn_lshift(t, t, 2 * n, norm);
            }
            else
            {
                mpn_lshift(t + 2 * n, x, n, norm);
                flint_mpn_mul_n(t, t + 2 * n, y, n);
            }
        }
        else
        {
            if (x == y)
                flint_mpn_sqr(t, x, n);
            else
                flint_mpn_mul_n(t, x, y, n);
        }

        flint_mpn_mul_or_mulhigh_n(t + 3 * n, t + n, dinv, n);
        mpn_add_n(t + 4 * n, t + 4 * n, t + n, n);

        flint_mpn_mul_n(t + 2 * n, t + 4 * n, d, n);
        cy = t[n] - t[3 * n] - mpn_sub_n(res, t, t + 2 * n, n);

        while (cy > 0)
            cy -= mpn_sub_n(res, res, d, n);

        if (mpn_cmp(res, d, n) >= 0)
            mpn_sub_n(res, res, d, n);

        if (norm)
            mpn_rshift(res, res, n, norm);
    }

    return GR_SUCCESS;
}

/* todo: check for 0? */
int
mpn_mod_mul_ui(mp_ptr res, mp_srcptr x, ulong y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS + 1];
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_size_t tn;

    t[n] = mpn_mul_1(t, x, n, y);
    tn = n + (t[n] != 0);
    mpn_mod_set_mpn(res, t, tn, ctx);
    return GR_SUCCESS;
}

/* todo: check for 0? */
#define UI_ABS_SI(x) (((slong)(x) < 0) ? (-(ulong)(x)) : ((ulong)(x)))

int
mpn_mod_mul_si(mp_ptr res, mp_srcptr x, slong y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS + 1];
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_size_t tn;

    t[n] = mpn_mul_1(t, x, n, UI_ABS_SI(y));
    tn = n + (t[n] != 0);
    mpn_mod_set_mpn(res, t, tn, ctx);
    if (y < 0)
        mpn_mod_neg(res, res, ctx);

    return GR_SUCCESS;
}

/* todo: check for 0? */
int
mpn_mod_mul_fmpz(mp_ptr res, mp_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*y))
    {
        mpn_mod_mul_si(res, x, *y, ctx);
    }
    else
    {
        mp_limb_t t[2 * MPN_MOD_MAX_LIMBS], cy;
        mp_size_t tn, n = MPN_MOD_CTX_NLIMBS(ctx);

        __mpz_struct * z = COEFF_TO_PTR(*y);
        mp_size_t ssize = z->_mp_size;
        mp_size_t zn = FLINT_ABS(ssize);
        mp_srcptr zd = z->_mp_d;

        if (zn <= n)
        {
            cy = flint_mpn_mul(t, x, n, zd, zn);
            tn = n + zn - (cy == 0);
            mpn_mod_set_mpn(res, t, tn, ctx);
        }
        else
        {
            mpn_mod_set_mpn(t, zd, zn, ctx);
            mpn_mod_mul(res, x, t, ctx);
        }

        if (ssize < 0)
            mpn_mod_neg(res, res, ctx);

    }

    return GR_SUCCESS;
}

int
mpn_mod_addmul(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    mpn_mod_mul(t, x, y, ctx);
    mpn_mod_add(res, res, t, ctx);
    return GR_SUCCESS;
}

int
mpn_mod_addmul_ui(mp_ptr res, mp_srcptr x, ulong y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    mpn_mod_mul_ui(t, x, y, ctx);
    mpn_mod_add(res, res, t, ctx);
    return GR_SUCCESS;
}

int
mpn_mod_addmul_si(mp_ptr res, mp_srcptr x, slong y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    mpn_mod_mul_si(t, x, y, ctx);
    mpn_mod_add(res, res, t, ctx);
    return GR_SUCCESS;
}

int
mpn_mod_addmul_fmpz(mp_ptr res, mp_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    mpn_mod_mul_fmpz(t, x, y, ctx);
    mpn_mod_add(res, res, t, ctx);
    return GR_SUCCESS;
}

int
mpn_mod_submul(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    mpn_mod_mul(t, x, y, ctx);
    mpn_mod_sub(res, res, t, ctx);
    return GR_SUCCESS;
}

int
mpn_mod_submul_ui(mp_ptr res, mp_srcptr x, ulong y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    mpn_mod_mul_ui(t, x, y, ctx);
    mpn_mod_sub(res, res, t, ctx);
    return GR_SUCCESS;
}

int
mpn_mod_submul_si(mp_ptr res, mp_srcptr x, slong y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    mpn_mod_mul_si(t, x, y, ctx);
    mpn_mod_sub(res, res, t, ctx);
    return GR_SUCCESS;
}

int
mpn_mod_submul_fmpz(mp_ptr res, mp_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    mpn_mod_mul_fmpz(t, x, y, ctx);
    mpn_mod_sub(res, res, t, ctx);
    return GR_SUCCESS;
}

int
mpn_mod_inv(mp_ptr res, mp_srcptr x, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_srcptr d = MPN_MOD_CTX_MODULUS(ctx);
    mp_limb_t g[MPN_MOD_MAX_LIMBS];
    mp_limb_t s[MPN_MOD_MAX_LIMBS];
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    mp_limb_t u[MPN_MOD_MAX_LIMBS];
    mp_size_t gsize, ssize;

    if (mpn_mod_is_one(x, ctx) == T_TRUE || mpn_mod_is_neg_one(x, ctx) == T_TRUE)
        return mpn_mod_set(res, x, ctx);

    flint_mpn_copyi(t, x, n);
    flint_mpn_copyi(u, d, n);
    /* todo: does mpn_gcdext allow aliasing? */
    gsize = mpn_gcdext(g, s, &ssize, t, n, u, n);

    if (gsize != 1 || g[0] != 1)
        return GR_DOMAIN;

    flint_mpn_zero(s + FLINT_ABS(ssize), n - FLINT_ABS(ssize));

    if (ssize >= 0)
        flint_mpn_copyi(res, s, n);
    else
        flint_mpn_negmod_n(res, s, d, n);

    return GR_SUCCESS;
}

int
mpn_mod_div(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx)
{
    int status;
    mp_limb_t t[MPN_MOD_MAX_LIMBS];

    status = mpn_mod_inv(t, y, ctx);
    if (status == GR_SUCCESS)
        mpn_mod_mul(res, x, t, ctx);
    else
        mpn_mod_zero(res, ctx);

    return status;
}
