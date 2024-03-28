/*
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

/* for wide add_ssss.... macros. todo; these ought to be provided
   everywhere */
#if FLINT_BITS == 64 && defined(__AVX2__)
#include "crt_helpers.h"
#endif

/* todo: separate method tables, more fixed-width optimizations */
/* todo: fast polynomial & matrix multiplication */
/* todo: powering */
/* todo: special functions */



#define UI_ABS_SI(x) (((slong)(x) < 0) ? (-(ulong)(x)) : ((ulong)(x)))

#define MPN_MOD_MIN_LIMBS 2
#define MPN_MOD_MAX_LIMBS 16

typedef struct
{
    mp_size_t nlimbs;
    mp_limb_t d[MPN_MOD_MAX_LIMBS];
    mp_limb_t dinv[MPN_MOD_MAX_LIMBS];
    mp_limb_t dnormed[MPN_MOD_MAX_LIMBS];
    flint_bitcnt_t norm;
    truth_t is_prime;
}
_mpn_mod_ctx_struct;

#define MPN_MOD_CTX(ctx) ((_mpn_mod_ctx_struct *)(GR_CTX_DATA_AS_PTR(ctx)))
#define MPN_MOD_CTX_NLIMBS(ctx) (MPN_MOD_CTX(ctx)->nlimbs)
#define MPN_MOD_CTX_MODULUS(ctx) (MPN_MOD_CTX(ctx)->d)
#define MPN_MOD_CTX_MODULUS_NORMED(ctx) (MPN_MOD_CTX(ctx)->dnormed)
#define MPN_MOD_CTX_MODULUS_PREINV(ctx) (MPN_MOD_CTX(ctx)->dinv)
#define MPN_MOD_CTX_NORM(ctx) (MPN_MOD_CTX(ctx)->norm)
#define MPN_MOD_CTX_IS_PRIME(ctx) (MPN_MOD_CTX(ctx)->is_prime)

FLINT_FORCE_INLINE
int flint_mpn_equal_p(mp_srcptr x, mp_srcptr y, mp_size_t xsize)
{
    slong i;
    for (i = 0; i < xsize; i++)
    {
        if (x[i] != y[i])
            return 0;
    }
    return 1;
}

FLINT_FORCE_INLINE void
flint_mpn_negmod_n(mp_ptr res, mp_srcptr x, mp_srcptr m, mp_size_t n)
{
    if (flint_mpn_zero_p(x, n))
        flint_mpn_zero(res, n);
    else
        mpn_sub_n(res, m, x, n);
}

FLINT_FORCE_INLINE void
flint_mpn_addmod_n(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m, mp_size_t n)
{
    mp_limb_t cy;
    cy = mpn_add_n(res, x, y, n);
    if (cy || mpn_cmp(res, m, n) >= 0)
        mpn_sub_n(res, res, m, n);
}

FLINT_FORCE_INLINE void
flint_mpn_submod_n(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m, mp_size_t n)
{
    int cmp = (mpn_cmp(x, y, n) < 0);
    mpn_sub_n(res, x, y, n);
    if (cmp)
        mpn_add_n(res, res, m, n);
}

/* assumes m <= n and y < m */
FLINT_FORCE_INLINE void
flint_mpn_addmod_n_m(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_size_t yn, mp_srcptr m, mp_size_t n)
{
    mp_limb_t cy;
    cy = mpn_add(res, x, n, y, yn);
    if (cy || mpn_cmp(res, m, n) >= 0)
        mpn_sub_n(res, res, m, n);
}

/* assumes m <= n and y < m */
FLINT_FORCE_INLINE void
flint_mpn_submod_n_m(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_size_t yn, mp_srcptr m, mp_size_t n)
{
    int cmp = (flint_mpn_zero_p(x + yn, n - yn) && mpn_cmp(x, y, yn) < 0);
    mpn_sub(res, x, n, y, yn);
    if (cmp)
        mpn_add_n(res, res, m, n);
}

FLINT_FORCE_INLINE void
flint_mpn_negmod_2(mp_ptr res, mp_srcptr x, mp_srcptr m)
{
    if (x[0] == 0 && x[1] == 0)
        res[1] = res[0] = 0;
    else
        sub_ddmmss(res[1], res[0], m[1], m[0], x[1], x[0]);
}

FLINT_FORCE_INLINE void
flint_mpn_addmod_2(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m)
{
    mp_limb_t cy;
    mp_limb_t m1 = m[1], m0 = m[0];
    add_sssaaaaaa(cy, res[1], res[0], 0, x[1], x[0], 0, y[1], y[0]);
    if (cy || (res[1] > m1 || (res[1] == m1 && res[0] >= m0)))
        sub_ddmmss(res[1], res[0], res[1], res[0], m1, m0);
}

FLINT_FORCE_INLINE void
flint_mpn_submod_2(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m)
{
    int cmp;
    mp_limb_t m1 = m[1], m0 = m[0];
    mp_limb_t x1 = x[1], x0 = x[0];
    mp_limb_t y1 = y[1], y0 = y[0];
    cmp = (x1 < y1) || (x1 == y1 && x0 < y0);
    sub_ddmmss(res[1], res[0], x1, x0, y1, y0);
    if (cmp)
        add_ssaaaa(res[1], res[0], res[1], res[0], m1, m0);
}

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
_gr_mpn_mod_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Integers mod ");
    gr_stream_write_free(out, _flint_mpn_get_str(MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx)));
    gr_stream_write(out, " (mpn)");
    return GR_SUCCESS;
}

void
_gr_mpn_mod_ctx_clear(gr_ctx_t ctx)
{
    flint_free(MPN_MOD_CTX(ctx));
}

truth_t
_gr_mpn_mod_ctx_is_field(gr_ctx_t ctx)
{
    return MPN_MOD_CTX_IS_PRIME(ctx);
}

void
_gr_mpn_mod_init(mp_ptr x, gr_ctx_t ctx)
{
    flint_mpn_zero(x, MPN_MOD_CTX_NLIMBS(ctx));
}

void
_gr_mpn_mod_clear(mp_ptr x, gr_ctx_t ctx)
{
}

void
_gr_mpn_mod_swap(mp_ptr x, mp_ptr y, gr_ctx_t ctx)
{
    slong i = 0, n = MPN_MOD_CTX_NLIMBS(ctx);
    for (i = 0; i < n; i++)
        FLINT_SWAP(mp_limb_t, x[i], y[i]);
}

int
_gr_mpn_mod_set(mp_ptr res, mp_srcptr x, gr_ctx_t ctx)
{
    flint_mpn_copyi(res, x, MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

int
_gr_mpn_mod_zero(mp_ptr res, gr_ctx_t ctx)
{
    flint_mpn_zero(res, MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

int
_gr_mpn_mod_one(mp_ptr res, gr_ctx_t ctx)
{
    res[0] = 1;
    flint_mpn_zero(res + 1, MPN_MOD_CTX_NLIMBS(ctx) - 1);
    return GR_SUCCESS;
}

int
_gr_mpn_mod_set_ui(mp_ptr res, ulong x, gr_ctx_t ctx)
{
    FLINT_ASSERT(MPN_MOD_CTX_NLIMBS(ctx) >= 2);

    res[0] = x;
    flint_mpn_zero(res + 1, MPN_MOD_CTX_NLIMBS(ctx) - 1);
    return GR_SUCCESS;
}

int
_gr_mpn_mod_set_si(mp_ptr res, slong x, gr_ctx_t ctx)
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
_gr_mpn_mod_neg_one(mp_ptr res, gr_ctx_t ctx)
{
    return _gr_mpn_mod_set_si(res, -1, ctx);
}

/* fixme: flint_mpn_mod_preinvn is misdocumented */

int
_gr_mpn_mod_set_mpn(mp_ptr res, mp_srcptr x, mp_size_t xn, gr_ctx_t ctx)
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
_gr_mpn_mod_set_fmpz(mp_ptr res, const fmpz_t x, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*x))
        _gr_mpn_mod_set_si(res, *x, ctx);
    else
    {
        mp_size_t nd = FLINT_ABS(COEFF_TO_PTR(*x)->_mp_size);
        int neg = COEFF_TO_PTR(*x)->_mp_size < 0;

        _gr_mpn_mod_set_mpn(res, COEFF_TO_PTR(*x)->_mp_d, nd, ctx);
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
_gr_mpn_mod_set_other(mp_ptr res, gr_ptr v, gr_ctx_t v_ctx, gr_ctx_t ctx)
{
    if (v_ctx == ctx)
        return _gr_mpn_mod_set(res, v, ctx);

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
                return _gr_mpn_mod_set_fmpz(res, v, ctx);
        }
    }

    return gr_generic_set_other(res, v, v_ctx, ctx);
}

int
_gr_mpn_mod_randtest(mp_ptr res, flint_rand_t state, gr_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_ui_array(t, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx));
    fmpz_randtest_mod(t, state, t);
    GR_IGNORE(_gr_mpn_mod_set_fmpz(res, t, ctx));
    fmpz_clear(t);
    return GR_SUCCESS;
}

int
_gr_mpn_mod_write(gr_stream_t out, mp_srcptr x, gr_ctx_t ctx)
{
    gr_stream_write_free(out, _flint_mpn_get_str(x, MPN_MOD_CTX_NLIMBS(ctx)));
    return GR_SUCCESS;
}

int
_gr_mpn_mod_get_fmpz(fmpz_t res, mp_srcptr x, gr_ctx_t ctx)
{
    fmpz_set_ui_array(res, x, MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}


truth_t
_gr_mpn_mod_is_zero(mp_srcptr x, gr_ctx_t ctx)
{
    return flint_mpn_zero_p(x, MPN_MOD_CTX_NLIMBS(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_mpn_mod_is_one(mp_srcptr x, gr_ctx_t ctx)
{
    return (x[0] == 1 && flint_mpn_zero_p(x + 1, MPN_MOD_CTX_NLIMBS(ctx) - 1)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_mpn_mod_is_neg_one(gr_srcptr x, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    _gr_mpn_mod_neg_one(t, ctx);
    return flint_mpn_equal_p(x, t, MPN_MOD_CTX_NLIMBS(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_mpn_mod_equal(mp_srcptr x, mp_srcptr y, gr_ctx_t ctx)
{
    return flint_mpn_equal_p(x, y, MPN_MOD_CTX_NLIMBS(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_mpn_mod_neg(mp_ptr res, mp_srcptr x, gr_ctx_t ctx)
{
    if (MPN_MOD_CTX_NLIMBS(ctx) == 2)
        flint_mpn_negmod_2(res, x, MPN_MOD_CTX_MODULUS(ctx));
    else
        flint_mpn_negmod_n(res, x, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

int
_gr_mpn_mod_add(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx)
{
    if (MPN_MOD_CTX_NLIMBS(ctx) == 2)
        flint_mpn_addmod_2(res, x, y, MPN_MOD_CTX_MODULUS(ctx));
    else
        flint_mpn_addmod_n(res, x, y, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

int
_gr_mpn_mod_sub(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx)
{
    if (MPN_MOD_CTX_NLIMBS(ctx) == 2)
        flint_mpn_submod_2(res, x, y, MPN_MOD_CTX_MODULUS(ctx));
    else
        flint_mpn_submod_n(res, x, y, MPN_MOD_CTX_MODULUS(ctx), MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

int
_gr_mpn_mod_add_ui(mp_ptr res, mp_srcptr x, ulong y, gr_ctx_t ctx)
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
_gr_mpn_mod_sub_ui(mp_ptr res, mp_srcptr x, ulong y, gr_ctx_t ctx)
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
_gr_mpn_mod_add_si(mp_ptr res, mp_srcptr x, slong y, gr_ctx_t ctx)
{
    if (y >= 0)
        _gr_mpn_mod_add_ui(res, x, (ulong) y, ctx);
    else
        _gr_mpn_mod_sub_ui(res, x, -(ulong) y, ctx);
    return GR_SUCCESS;
}

int
_gr_mpn_mod_sub_si(mp_ptr res, mp_srcptr x, slong y, gr_ctx_t ctx)
{
    if (y >= 0)
        _gr_mpn_mod_sub_ui(res, x, (ulong) y, ctx);
    else
        _gr_mpn_mod_add_ui(res, x, -(ulong) y, ctx);
    return GR_SUCCESS;
}

int
_gr_mpn_mod_add_fmpz(mp_ptr res, mp_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*y))
    {
        _gr_mpn_mod_add_si(res, x, *y, ctx);
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
            _gr_mpn_mod_set_mpn(t, zd, zn, ctx);
            if (ssize > 0)
                _gr_mpn_mod_add(res, x, t, ctx);
            else
                _gr_mpn_mod_sub(res, x, t, ctx);
        }
    }

    return GR_SUCCESS;
}

int
_gr_mpn_mod_sub_fmpz(mp_ptr res, mp_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*y))
    {
        _gr_mpn_mod_sub_si(res, x, *y, ctx);
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
            _gr_mpn_mod_set_mpn(t, zd, zn, ctx);
            if (ssize > 0)
                _gr_mpn_mod_sub(res, x, t, ctx);
            else
                _gr_mpn_mod_add(res, x, t, ctx);
        }
    }

    return GR_SUCCESS;
}


/* should mirror flint_mpn_mulmod_preinvn (backport any improvements) */
int
_gr_mpn_mod_mul(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx)
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

int
_gr_mpn_mod_sqr(mp_ptr res, mp_srcptr x, gr_ctx_t ctx)
{
    return _gr_mpn_mod_mul(res, x, x, ctx);
}


/* todo: check for 0? */
int
_gr_mpn_mod_mul_ui(mp_ptr res, mp_srcptr x, ulong y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS + 1];
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_size_t tn;

    t[n] = mpn_mul_1(t, x, n, y);
    tn = n + (t[n] != 0);
    _gr_mpn_mod_set_mpn(res, t, tn, ctx);
    return GR_SUCCESS;
}

/* todo: check for 0? */
int
_gr_mpn_mod_mul_si(mp_ptr res, mp_srcptr x, slong y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS + 1];
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_size_t tn;

    t[n] = mpn_mul_1(t, x, n, UI_ABS_SI(y));
    tn = n + (t[n] != 0);
    _gr_mpn_mod_set_mpn(res, t, tn, ctx);
    if (y < 0)
        _gr_mpn_mod_neg(res, res, ctx);

    return GR_SUCCESS;
}

/* todo: check for 0? */
int
_gr_mpn_mod_mul_fmpz(mp_ptr res, mp_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*y))
    {
        _gr_mpn_mod_mul_si(res, x, *y, ctx);
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
            _gr_mpn_mod_set_mpn(res, t, tn, ctx);
        }
        else
        {
            _gr_mpn_mod_set_mpn(t, zd, zn, ctx);
            _gr_mpn_mod_mul(res, x, t, ctx);
        }

        if (ssize < 0)
            _gr_mpn_mod_neg(res, res, ctx);

    }

    return GR_SUCCESS;
}

int
_gr_mpn_mod_addmul(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    _gr_mpn_mod_mul(t, x, y, ctx);
    _gr_mpn_mod_add(res, res, t, ctx);
    return GR_SUCCESS;
}

int
_gr_mpn_mod_addmul_ui(mp_ptr res, mp_srcptr x, ulong y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    _gr_mpn_mod_mul_ui(t, x, y, ctx);
    _gr_mpn_mod_add(res, res, t, ctx);
    return GR_SUCCESS;
}

int
_gr_mpn_mod_addmul_si(mp_ptr res, mp_srcptr x, slong y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    _gr_mpn_mod_mul_si(t, x, y, ctx);
    _gr_mpn_mod_add(res, res, t, ctx);
    return GR_SUCCESS;
}

int
_gr_mpn_mod_addmul_fmpz(mp_ptr res, mp_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    _gr_mpn_mod_mul_fmpz(t, x, y, ctx);
    _gr_mpn_mod_add(res, res, t, ctx);
    return GR_SUCCESS;
}

int
_gr_mpn_mod_submul(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    _gr_mpn_mod_mul(t, x, y, ctx);
    _gr_mpn_mod_sub(res, res, t, ctx);
    return GR_SUCCESS;
}

int
_gr_mpn_mod_submul_ui(mp_ptr res, mp_srcptr x, ulong y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    _gr_mpn_mod_mul_ui(t, x, y, ctx);
    _gr_mpn_mod_sub(res, res, t, ctx);
    return GR_SUCCESS;
}

int
_gr_mpn_mod_submul_si(mp_ptr res, mp_srcptr x, slong y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    _gr_mpn_mod_mul_si(t, x, y, ctx);
    _gr_mpn_mod_sub(res, res, t, ctx);
    return GR_SUCCESS;
}

int
_gr_mpn_mod_submul_fmpz(mp_ptr res, mp_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    _gr_mpn_mod_mul_fmpz(t, x, y, ctx);
    _gr_mpn_mod_sub(res, res, t, ctx);
    return GR_SUCCESS;
}

int
_gr_mpn_mod_inv(mp_ptr res, mp_srcptr x, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_srcptr d = MPN_MOD_CTX_MODULUS(ctx);
    mp_limb_t g[MPN_MOD_MAX_LIMBS];
    mp_limb_t s[MPN_MOD_MAX_LIMBS];
    mp_limb_t t[MPN_MOD_MAX_LIMBS];
    mp_limb_t u[MPN_MOD_MAX_LIMBS];
    mp_size_t gsize, ssize;

    if (_gr_mpn_mod_is_one(x, ctx) == T_TRUE || _gr_mpn_mod_is_neg_one(x, ctx) == T_TRUE)
        return _gr_mpn_mod_set(res, x, ctx);

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
_gr_mpn_mod_div(mp_ptr res, mp_srcptr x, mp_srcptr y, gr_ctx_t ctx)
{
    int status;
    mp_limb_t t[MPN_MOD_MAX_LIMBS];

    status = _gr_mpn_mod_inv(t, y, ctx);
    if (status == GR_SUCCESS)
        _gr_mpn_mod_mul(res, x, t, ctx);
    else
        _gr_mpn_mod_zero(res, ctx);

    return status;
}

int
_gr_mpn_mod_vec_clear(mp_ptr res, slong len, gr_ctx_t ctx)
{
    return GR_SUCCESS;
}

int
_gr_mpn_mod_vec_zero(mp_ptr res, slong len, gr_ctx_t ctx)
{
    flint_mpn_zero(res, len * MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

int
_gr_mpn_mod_vec_set(mp_ptr res, mp_srcptr x, slong len, gr_ctx_t ctx)
{
    flint_mpn_copyi(res, x, len * MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

void
_gr_mpn_mod_vec_swap(mp_ptr vec1, mp_ptr vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    for (i = 0; i < len * n; i++)
        FLINT_SWAP(mp_limb_t, vec1[i], vec2[i]);
}

int
_gr_mpn_mod_vec_neg(mp_ptr res, mp_srcptr x, slong len, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_srcptr d = MPN_MOD_CTX_MODULUS(ctx);
    slong i;

    if (n == 2)
    {
        /* Only read to registers once */
        mp_limb_t dd[2];
        dd[0] = d[0];
        dd[1] = d[1];

        for (i = 0; i < len; i++)
            flint_mpn_negmod_2(res + i * n, x + i * n, dd);
    }
    else
        for (i = 0; i < len; i++)
            flint_mpn_negmod_n(res + i * n, x + i * n, d, n);

    return GR_SUCCESS;
}

int
_gr_mpn_mod_vec_add(mp_ptr res, mp_srcptr x, mp_srcptr y, slong len, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_srcptr d = MPN_MOD_CTX_MODULUS(ctx);
    slong i;

    if (n == 2)
    {
        /* Only read to registers once */
        mp_limb_t dd[2];
        dd[0] = d[0];
        dd[1] = d[1];

        for (i = 0; i < len; i++)
            flint_mpn_addmod_2(res + i * n, x + i * n, y + i * n, dd);
    }
    else
        for (i = 0; i < len; i++)
            flint_mpn_addmod_n(res + i * n, x + i * n, y + i * n, d, n);

    return GR_SUCCESS;
}

int
_gr_mpn_mod_vec_sub(mp_ptr res, mp_srcptr x, mp_srcptr y, slong len, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_srcptr d = MPN_MOD_CTX_MODULUS(ctx);
    slong i;

    if (n == 2)
    {
        /* Only read to registers once */
        mp_limb_t dd[2];
        dd[0] = d[0];
        dd[1] = d[1];

        for (i = 0; i < len; i++)
            flint_mpn_submod_2(res + i * n, x + i * n, y + i * n, dd);
    }
    else
        for (i = 0; i < len; i++)
            flint_mpn_submod_n(res + i * n, x + i * n, y + i * n, d, n);

    return GR_SUCCESS;
}

int
_gr_mpn_mod_vec_mul(mp_ptr res, mp_srcptr x, mp_srcptr y, slong len, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    slong i;

    if (n == 2)
    {
        mp_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
        mp_srcptr dnormed = MPN_MOD_CTX_MODULUS_NORMED(ctx);
        mp_srcptr dinv = MPN_MOD_CTX_MODULUS_PREINV(ctx);

        for (i = 0; i < len; i++)
            flint_mpn_mulmod_preinvn_2(res + i * n, x + i * n, y + i * n, dnormed, dinv, norm);
    }
    else
        for (i = 0; i < len; i++)
            _gr_mpn_mod_mul(res + i * n, x + i * n, y + i * n, ctx);

    return GR_SUCCESS;
}

/* todo: worth it to check for special cases (0, 1)? */
/* todo: shoup multiplication */
int
_gr_mpn_mod_vec_mul_scalar(mp_ptr res, mp_srcptr x, slong len, mp_srcptr y, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    slong i;

    if (n == 2)
    {
        mp_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
        mp_srcptr dnormed = MPN_MOD_CTX_MODULUS_NORMED(ctx);
        mp_srcptr dinv = MPN_MOD_CTX_MODULUS_PREINV(ctx);

        for (i = 0; i < len; i++)
            flint_mpn_mulmod_preinvn_2(res + i * n, x + i * n, y, dnormed, dinv, norm);
    }
    else
        for (i = 0; i < len; i++)
            _gr_mpn_mod_mul(res + i * n, x + i * n, y, ctx);

    return GR_SUCCESS;
}

int
_gr_mpn_mod_scalar_mul_vec(mp_ptr res, mp_srcptr y, mp_srcptr x, slong len, gr_ctx_t ctx)
{
    return _gr_mpn_mod_vec_mul_scalar(res, x, len, y, ctx);
}

/* todo: worth it to check for special cases (0, 1)? */
/* todo: shoup multiplication */
int
_gr_mpn_mod_vec_addmul_scalar(mp_ptr res, mp_srcptr x, slong len, mp_srcptr y, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    slong i;

    if (n == 2)
    {
        mp_limb_t t[2];
        mp_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
        mp_srcptr dnormed = MPN_MOD_CTX_MODULUS_NORMED(ctx);
        mp_srcptr dinv = MPN_MOD_CTX_MODULUS_PREINV(ctx);
        mp_srcptr d = MPN_MOD_CTX_MODULUS(ctx);

        for (i = 0; i < len; i++)
        {
            flint_mpn_mulmod_preinvn_2(t, x + i * n, y, dnormed, dinv, norm);
            flint_mpn_addmod_2(res + i * n, res + i * n, t, d);
        }
    }
    else
    {
        mp_limb_t t[MPN_MOD_MAX_LIMBS];

        for (i = 0; i < len; i++)
        {
            _gr_mpn_mod_mul(t, x + i * n, y, ctx);
            _gr_mpn_mod_add(res + i * n, res + i * n, t, ctx);
        }
    }

    return GR_SUCCESS;
}

/* todo: length 1 */
int
_gr_mpn_mod_vec_dot(mp_ptr res, mp_srcptr initial, int subtract, mp_srcptr vec1, mp_srcptr vec2, slong len, gr_ctx_t ctx)
{
    mp_limb_t s[2 * MPN_MOD_MAX_LIMBS + 1];
    mp_limb_t t[2 * MPN_MOD_MAX_LIMBS];
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_size_t sn;
    slong i;

    if (len <= 0)
    {
        if (initial == NULL)
            flint_mpn_zero(res, n);
        else
            flint_mpn_copyi(res, initial, n);
        return GR_SUCCESS;
    }

    switch (n)
    {
        case 2:
            FLINT_MPN_MUL_2X2(s[3], s[2], s[1], s[0], vec1[1], vec1[0], vec2[1], vec2[0]);
            s[4] = 0;
            for (i = 1; i < len; i++)
            {
                FLINT_MPN_MUL_2X2(t[3], t[2], t[1], t[0], vec1[2 * i + 1], vec1[2 * i], vec2[2 * i + 1], vec2[2 * i]);
#if defined(add_sssssaaaaaaaaaa)
                add_sssssaaaaaaaaaa(s[4], s[3], s[2], s[1], s[0], s[4], s[3], s[2], s[1], s[0], 0, t[3], t[2], t[1], t[0]);
#else
                s[2 * n] += mpn_add_n(s, s, t, 2 * n);
#endif
            }

            break;

        default:
            flint_mpn_mul_n(s, vec1, vec2, n);
            s[2 * n] = 0;
            for (i = 1; i < len; i++)
            {
                flint_mpn_mul_n(t, vec1 + i * n, vec2 + i * n, n);
                s[2 * n] += mpn_add_n(s, s, t, 2 * n);
            }
    }

    sn = 2 * n + (s[2 * n] != 0);

    if (initial == NULL)
    {
        if (!subtract)
        {
            _gr_mpn_mod_set_mpn(res, s, sn, ctx);
        }
        else
        {
            _gr_mpn_mod_set_mpn(t, s, sn, ctx);
            _gr_mpn_mod_neg(res, t, ctx);
        }
    }
    else
    {
        /* todo: consider setting initial as a starting value for the sum */
        _gr_mpn_mod_set_mpn(t, s, sn, ctx);

        if (!subtract)
            _gr_mpn_mod_add(res, initial, t, ctx);
        else
            _gr_mpn_mod_sub(res, initial, t, ctx);
    }

    return GR_SUCCESS;
}

int
_gr_mpn_mod_vec_dot_rev(mp_ptr res, mp_srcptr initial, int subtract, mp_srcptr vec1, mp_srcptr vec2, slong len, gr_ctx_t ctx)
{
    mp_limb_t s[2 * MPN_MOD_MAX_LIMBS + 1];
    mp_limb_t t[2 * MPN_MOD_MAX_LIMBS];
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_size_t sn;
    slong i;

    if (len <= 0)
    {
        if (initial == NULL)
            flint_mpn_zero(res, n);
        else
            flint_mpn_copyi(res, initial, n);
        return GR_SUCCESS;
    }

    switch (n)
    {
        case 2:
            FLINT_MPN_MUL_2X2(s[3], s[2], s[1], s[0], vec1[1], vec1[0], vec2[2 * len - 1], vec2[2 * len - 2]);
            s[4] = 0;
            for (i = 1; i < len; i++)
            {
                FLINT_MPN_MUL_2X2(t[3], t[2], t[1], t[0], vec1[2 * i + 1], vec1[2 * i], vec2[2 * (len - i - 1) + 1], vec2[2 * (len - i - 1)]);
#if defined(add_sssssaaaaaaaaaa)
                add_sssssaaaaaaaaaa(s[4], s[3], s[2], s[1], s[0], s[4], s[3], s[2], s[1], s[0], 0, t[3], t[2], t[1], t[0]);
#else
                s[2 * n] += mpn_add_n(s, s, t, 2 * n);
#endif
            }

            break;

        default:
            flint_mpn_mul_n(s, vec1, vec2 + (len - 1) * n, n);
            s[2 * n] = 0;
            for (i = 1; i < len; i++)
            {
                flint_mpn_mul_n(t, vec1 + i * n, vec2 + (len - i - 1) * n, n);
                s[2 * n] += mpn_add_n(s, s, t, 2 * n);
            }
    }

    sn = 2 * n + (s[2 * n] != 0);

    if (initial == NULL)
    {
        if (!subtract)
        {
            _gr_mpn_mod_set_mpn(res, s, sn, ctx);
        }
        else
        {
            _gr_mpn_mod_set_mpn(t, s, sn, ctx);
            _gr_mpn_mod_neg(res, t, ctx);
        }
    }
    else
    {
        /* todo: consider setting initial as a starting value for the sum */
        _gr_mpn_mod_set_mpn(t, s, sn, ctx);

        if (!subtract)
            _gr_mpn_mod_add(res, initial, t, ctx);
        else
            _gr_mpn_mod_sub(res, initial, t, ctx);
    }

    return GR_SUCCESS;
}

/* todo: retune this when arithmetic is more optimized */
int
_gr_mpn_mod_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
{
    slong cutoff;

    if (MPN_MOD_CTX_NLIMBS(ctx) == 2)
        cutoff = 14;
    else
        cutoff = 4;

    if (gr_mat_nrows(A, ctx) <= 2 * cutoff || gr_mat_ncols(A, ctx) <= 2 * cutoff)
        return gr_mat_lu_classical(rank, P, LU, A, rank_check, ctx);
    else
        return gr_mat_lu_recursive(rank, P, LU, A, rank_check, cutoff, ctx);
}


int _mpn_mod_methods_initialized = 0;

gr_static_method_table _mpn_mod_methods;

gr_method_tab_input _mpn_mod_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_mpn_mod_ctx_write},
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_mpn_mod_ctx_clear},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) _gr_mpn_mod_ctx_is_field},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) _gr_mpn_mod_ctx_is_field},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) _gr_mpn_mod_ctx_is_field},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_mpn_mod_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_mpn_mod_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_mpn_mod_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_mpn_mod_set},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_mpn_mod_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_mpn_mod_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_mpn_mod_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_mpn_mod_one},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) _gr_mpn_mod_neg_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_mpn_mod_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_mpn_mod_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_mpn_mod_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_mpn_mod_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_mpn_mod_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_mpn_mod_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_mpn_mod_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_mpn_mod_set_fmpz},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_mpn_mod_set_other},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_mpn_mod_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_mpn_mod_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_mpn_mod_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_mpn_mod_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_mpn_mod_add_fmpz},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_mpn_mod_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_mpn_mod_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_mpn_mod_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_mpn_mod_sub_fmpz},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_mpn_mod_mul},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_mpn_mod_mul_si},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_mpn_mod_mul_ui},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_mpn_mod_mul_fmpz},
    {GR_METHOD_ADDMUL,          (gr_funcptr) _gr_mpn_mod_addmul},
    {GR_METHOD_ADDMUL_SI,       (gr_funcptr) _gr_mpn_mod_addmul_si},
    {GR_METHOD_ADDMUL_UI,       (gr_funcptr) _gr_mpn_mod_addmul_ui},
    {GR_METHOD_ADDMUL_FMPZ,     (gr_funcptr) _gr_mpn_mod_addmul_fmpz},
    {GR_METHOD_SUBMUL,          (gr_funcptr) _gr_mpn_mod_submul},
    {GR_METHOD_SUBMUL_SI,       (gr_funcptr) _gr_mpn_mod_submul_si},
    {GR_METHOD_SUBMUL_UI,       (gr_funcptr) _gr_mpn_mod_submul_ui},
    {GR_METHOD_SUBMUL_FMPZ,     (gr_funcptr) _gr_mpn_mod_submul_fmpz},

/*
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_mpn_mod_mul_two},
*/
    {GR_METHOD_SQR,             (gr_funcptr) _gr_mpn_mod_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_mpn_mod_div},
/*
    {GR_METHOD_DIV_NONUNIQUE,   (gr_funcptr) _gr_mpn_mod_div_nonunique},
    {GR_METHOD_DIVIDES,         (gr_funcptr) _gr_mpn_mod_divides},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_mpn_mod_is_invertible},
*/


    {GR_METHOD_INV,             (gr_funcptr) _gr_mpn_mod_inv},
/*
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_mpn_mod_pow_si},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_mpn_mod_pow_ui},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_mpn_mod_pow_fmpz},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_mpn_mod_sqrt},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_mpn_mod_is_square},
*/

    {GR_METHOD_VEC_INIT,        (gr_funcptr) _gr_mpn_mod_vec_zero},
    {GR_METHOD_VEC_CLEAR,       (gr_funcptr) _gr_mpn_mod_vec_clear},
    {GR_METHOD_VEC_SET,         (gr_funcptr) _gr_mpn_mod_vec_set},
    {GR_METHOD_VEC_SWAP,        (gr_funcptr) _gr_mpn_mod_vec_swap},
    {GR_METHOD_VEC_ZERO,        (gr_funcptr) _gr_mpn_mod_vec_zero},
    {GR_METHOD_VEC_NEG,         (gr_funcptr) _gr_mpn_mod_vec_neg},
    {GR_METHOD_VEC_ADD,         (gr_funcptr) _gr_mpn_mod_vec_add},
    {GR_METHOD_VEC_SUB,         (gr_funcptr) _gr_mpn_mod_vec_sub},
    {GR_METHOD_VEC_MUL,         (gr_funcptr) _gr_mpn_mod_vec_mul},
    {GR_METHOD_VEC_MUL_SCALAR,  (gr_funcptr) _gr_mpn_mod_vec_mul_scalar},
    {GR_METHOD_VEC_ADDMUL_SCALAR,    (gr_funcptr) _gr_mpn_mod_vec_addmul_scalar},
    {GR_METHOD_SCALAR_MUL_VEC,  (gr_funcptr) _gr_mpn_mod_scalar_mul_vec},

    {GR_METHOD_VEC_DOT,         (gr_funcptr) _gr_mpn_mod_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _gr_mpn_mod_vec_dot_rev},

/*
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_mpn_mod_poly_mullow},
    {GR_METHOD_POLY_INV_SERIES, (gr_funcptr) _gr_mpn_mod_poly_inv_series},
    {GR_METHOD_POLY_DIV_SERIES, (gr_funcptr) _gr_mpn_mod_poly_div_series},
    {GR_METHOD_POLY_DIVREM,     (gr_funcptr) _gr_mpn_mod_poly_divrem},
    {GR_METHOD_POLY_ROOTS,      (gr_funcptr) _gr_mpn_mod_roots_gr_poly},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_mpn_mod_mat_mul},
*/

    {GR_METHOD_MAT_LU,          (gr_funcptr) _gr_mpn_mod_mat_lu},

/*
    {GR_METHOD_MAT_DET,         (gr_funcptr) _gr_mpn_mod_mat_det},
*/
    {0,                         (gr_funcptr) NULL},
};

/* todo: allow passing modulus as an mpn */
int
gr_ctx_init_mpn_mod(gr_ctx_t ctx, const fmpz_t n)
{
    mp_size_t s;
    mp_srcptr nptr;
    mp_bitcnt_t norm;

    s = fmpz_size(n);

    if (s < MPN_MOD_MIN_LIMBS || s > MPN_MOD_MAX_LIMBS || fmpz_sgn(n) < 0)
        return GR_DOMAIN;

    ctx->which_ring = GR_CTX_MPN_MOD;
    ctx->sizeof_elem = s * sizeof(mp_limb_t);

    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(_mpn_mod_ctx_struct));

    MPN_MOD_CTX_NLIMBS(ctx) = s;

    nptr = COEFF_TO_PTR(*n)->_mp_d;
    flint_mpn_copyi(MPN_MOD_CTX_MODULUS(ctx), nptr, s);
    MPN_MOD_CTX_NORM(ctx) = norm = flint_clz(nptr[s - 1]);

    if (norm == 0)
        flint_mpn_copyi(MPN_MOD_CTX_MODULUS_NORMED(ctx), nptr, s);
    else
        mpn_lshift(MPN_MOD_CTX_MODULUS_NORMED(ctx), nptr, s, norm);

    flint_mpn_preinvn(MPN_MOD_CTX_MODULUS_PREINV(ctx), MPN_MOD_CTX_MODULUS_NORMED(ctx), s);

    MPN_MOD_CTX_IS_PRIME(ctx) = T_UNKNOWN;

    ctx->size_limit = WORD_MAX;

    ctx->methods = _mpn_mod_methods;

    if (!_mpn_mod_methods_initialized)
    {
        gr_method_tab_init(_mpn_mod_methods, _mpn_mod_methods_input);
        _mpn_mod_methods_initialized = 1;
    }

    return GR_SUCCESS;
}

/* todo: have a generic interface for this */
void
gr_ctx_mpn_mod_set_primality(gr_ctx_t ctx, truth_t is_prime)
{
    MPN_MOD_CTX_IS_PRIME(ctx) = is_prime;
}
