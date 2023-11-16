/*
    Copyright (C) 2015 Vladimir Glazachev
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_generic.h"

/*
    Returns smallest integer k satisfying:
        log(n) < (k * (k + 1) * 2^(2 * k)) / (2^(k + 1) - k - 2) + 1
*/
static ulong
sliding_select_k(ulong bits)
{
    if (bits <= 8)  return 1;
    if (bits <= 24) return 2;
    if (bits <= 69) return 3;
    if (bits <= 196) return 4;
    if (bits <= 538) return 5;
    if (bits <= 1433) return 6;
    if (bits <= 3714) return 7;
    if (bits <= 9399) return 8;
    if (bits <= 23290) return 9;
    if (bits <= 56651) return 10;
    return 11;
}

#define MPN_TSTBIT(x, xn, b) (((b) / FLINT_BITS >= ((xn))) ? 0 : ((x[(b) / FLINT_BITS] >> ((b) % FLINT_BITS)) & UWORD(1)))

/* todo: avoid swaps (or perform pointer swaps) */
/* note: supports aliasing */
static int
_gr_pow_mpn_sliding(gr_ptr f, gr_srcptr g, mp_srcptr exp, mp_size_t en, gr_ctx_t ctx)
{
    slong h, k, value;
    slong i, j, alloc;
    gr_ptr temp;
    fmpz * g_powers;
    slong sz = ctx->sizeof_elem;
    mp_bitcnt_t ebits;
    int status = GR_SUCCESS;

    ebits = (en - 1) * FLINT_BITS + FLINT_BIT_COUNT(exp[en - 1]);

    /* selects optimal k value for n */
    k = sliding_select_k(ebits);

    /*
        g_powers store odd powers of g up to 2^k - 1;
        g_powers[(i + 1) / 2] = g^i
    */
    alloc = (WORD(1) << (k - 1)) + 1;

    GR_TMP_INIT_VEC(g_powers, alloc + 1, ctx);
    temp = GR_ENTRY(g_powers, alloc, sz);

    /* temp = g * g */
    status |= gr_sqr(temp, g, ctx);

    status |= gr_one(GR_ENTRY(g_powers, 0, sz), ctx);
    status |= gr_set(GR_ENTRY(g_powers, 1, sz), g, ctx);

    /* sets g_powers[i] = g^2 * g_powers[i - 1] */
    for (i = 2; i <= WORD(1) << (k - 1); i++)
        status |= gr_mul(GR_ENTRY(g_powers, i, sz), GR_ENTRY(g_powers, i - 1, sz), temp, ctx);

    status |= gr_one(f, ctx);

    i = ebits - 1;

    /* working with pow = (e_l, e_{l-1}, ... , e_0) in 2 base */
    while (i >= 0)
    {
        if (!MPN_TSTBIT(exp, en, i))
        {
            status |= gr_sqr(temp, f, ctx);
            gr_swap(f, temp, ctx);
            i--;
        }
        else
        {
            /*
                finds length of chain; chain is length of
                longest bitstring less then k ending on 1
            */
            j = FLINT_MAX(i - k + 1, 0);
            while (MPN_TSTBIT(exp, en, j) == 0 && j <= i)
                j++;

            /* f = f^(2^(i - j + 1)) */
            for (h = 0; h < i - j + 1; h++)
            {
                status |= gr_sqr(temp, f, ctx);
                gr_swap(f, temp, ctx);
            }

            /* value = binary number (e_i, ... , e_j) */
            value = 0;
            for (h = 0; h < i - j + 1; h++)
                value += MPN_TSTBIT(exp, en, j + h) << h;

            /* f = f * g^value */
            status |= gr_mul(temp, f, GR_ENTRY(g_powers, (value + 1) / 2, sz), ctx);
            gr_swap(f, temp, ctx);

            /* increase i */
            i = j - 1;
        }
    }

    GR_TMP_CLEAR_VEC(g_powers, alloc + 1, ctx);

    return status;
}

int
gr_generic_pow_fmpz_sliding(gr_ptr f, gr_srcptr g, const fmpz_t pow, gr_ctx_t ctx)
{
    if (fmpz_sgn(pow) < 0)
        return GR_UNABLE;

    if (fmpz_is_zero(pow))
        return gr_one(f, ctx);

    if (!COEFF_IS_MPZ(*pow))
    {
        ulong t = *pow;
        return _gr_pow_mpn_sliding(f, g, &t, 1, ctx);
    }
    else
    {
        return _gr_pow_mpn_sliding(f, g, COEFF_TO_PTR(*pow)->_mp_d, COEFF_TO_PTR(*pow)->_mp_size, ctx);
    }
}

int
gr_generic_pow_ui_sliding(gr_ptr f, gr_srcptr g, ulong pow, gr_ctx_t ctx)
{
    if (pow == 0)
        return gr_one(f, ctx);
    else
        return _gr_pow_mpn_sliding(f, g, &pow, 1, ctx);
}

/* Assumes exp >= 2; res and tmp not not aliased with x. */
int
_gr_generic_pow_ui_binexp(gr_ptr res, gr_ptr tmp, gr_srcptr x, ulong exp, gr_ctx_t ctx)
{
    gr_ptr R, S, T;
    gr_method_unary_op sqr = GR_UNARY_OP(ctx, SQR);
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    int status;
    int zeros;
    ulong bit;

    status = GR_SUCCESS;

    /* Determine parity due to swaps */
    zeros = 0;
    bit = exp;
    while (bit > 1)
    {
        zeros += !(bit & 1);
        bit >>= 1;
    }

    if (zeros % 2)
    {
        R = res;
        S = tmp;
    }
    else
    {
        R = tmp;
        S = res;
    }

    bit = UWORD(1) << (FLINT_BIT_COUNT(exp) - 2);

    status |= sqr(R, x, ctx);

    if (bit & exp)
    {
        status |= mul(S, R, x, ctx);
        T = R;
        R = S;
        S = T;
    }

    while (bit >>= 1)
    {
        status |= sqr(S, R, ctx);

        if (bit & exp)
        {
            status |= mul(R, S, x, ctx);
        }
        else
        {
            T = R;
            R = S;
            S = T;
        }
    }

    return status;
}

static int
gr_generic_pow3(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status;

    if (res == x)
    {
        gr_ptr t;
        GR_TMP_INIT(t, ctx);
        status = gr_sqr(t, x, ctx);
        status |= gr_mul(res, t, x, ctx);
        GR_TMP_CLEAR(t, ctx);
    }
    else
    {
        status = gr_sqr(res, x, ctx);
        status |= gr_mul(res, res, x, ctx);
    }

    return status;
}

static int
gr_generic_pow4(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status;
    status = gr_sqr(res, x, ctx);
    status |= gr_sqr(res, res, ctx);
    return status;
}

/* todo: optimize swaps */
int
gr_generic_pow_fmpz_binexp(gr_ptr res, gr_srcptr x, const fmpz_t exp, gr_ctx_t ctx)
{
    gr_ptr t, u;
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_unary_op sqr = GR_UNARY_OP(ctx, SQR);
    gr_method_swap_op swap = GR_SWAP_OP(ctx, SWAP);
    int status;
    slong i;

    if (*exp == 0)
        return gr_one(res, ctx);
    else if (*exp == 1)
        return gr_set(res, x, ctx);
    else if (*exp == 2)
        return gr_sqr(res, x, ctx);
    else if (*exp == 3)
        return gr_generic_pow3(res, x, ctx);
    else if (*exp == 4)
        return gr_generic_pow4(res, x, ctx);

    if (fmpz_sgn(exp) < 0)
        return GR_UNABLE;

    status = GR_SUCCESS;

    GR_TMP_INIT2(t, u, ctx);

    status |= gr_set(t, x, ctx);

    for (i = fmpz_bits(exp) - 2; i >= 0; i--)
    {
        status |= sqr(u, t, ctx);

        if (fmpz_tstbit(exp, i))
            status |= mul(t, u, x, ctx);
        else
            swap(t, u, ctx);
    }

    swap(res, t, ctx);

    GR_TMP_CLEAR2(t, u, ctx);

    return status;
}

int
gr_generic_pow_ui_binexp(gr_ptr res, gr_srcptr x, ulong e, gr_ctx_t ctx)
{
    int status;
    gr_ptr t, u;


    if (e == 0)
        return gr_one(res, ctx);
    else if (e == 1)
        return gr_set(res, x, ctx);
    else if (e == 2)
        return gr_sqr(res, x, ctx);
    else if (e == 3)
        return gr_generic_pow3(res, x, ctx);
    else if (e == 4)
        return gr_generic_pow4(res, x, ctx);

    if (res == x)
    {
        GR_TMP_INIT2(t, u, ctx);
        status = gr_set(u, x, ctx);
        status |= _gr_generic_pow_ui_binexp(res, t, u, e, ctx);
        GR_TMP_CLEAR2(t, u, ctx);
    }
    else
    {
        GR_TMP_INIT(t, ctx);
        status = _gr_generic_pow_ui_binexp(res, t, x, e, ctx);
        GR_TMP_CLEAR(t, ctx);
    }

    return status;
}

int
gr_generic_pow_ui(gr_ptr res, gr_srcptr x, ulong e, gr_ctx_t ctx)
{
    return gr_generic_pow_ui_binexp(res, x, e, ctx);
}

int
gr_generic_pow_si(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx)
{
    if (e >= 0)
    {
        return gr_pow_ui(res, x, e, ctx);
    }
    else
    {
        int status;

        status = gr_inv(res, x, ctx);

        if (status == GR_SUCCESS && e != -1)
            status = gr_pow_ui(res, res, -e, ctx);

        return status;
    }
}

int
gr_generic_pow_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t e, gr_ctx_t ctx)
{
    int status;

    if (fmpz_sgn(e) < 0)
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_neg(f, e);

        /* todo: some heuristic for when we want to invert before/after powering */
        status = gr_inv(res, x, ctx);
        if (status == GR_SUCCESS)
            status = gr_generic_pow_fmpz(res, res, f, ctx);

        fmpz_clear(f);
        return status;
    }

    return gr_generic_pow_fmpz_binexp(res, x, e, ctx);
}
