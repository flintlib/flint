/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpq.h"
#include "arith.h"
#include "arb.h"
#include "gr_special.h"

/* todo: maybe add a second table for a second word */

#if FLINT_BITS == 64
#define FAC_TAB_SIZE 21
#else
#define FAC_TAB_SIZE 13
#endif

static const mp_limb_t fac_tab[] =
{
  UWORD(1), UWORD(1), UWORD(2), UWORD(6), UWORD(24), UWORD(120), UWORD(720), UWORD(5040), UWORD(40320), UWORD(362880),
  UWORD(3628800), UWORD(39916800), UWORD(479001600),
#if FLINT_BITS == 64
  UWORD(6227020800), UWORD(87178291200), UWORD(1307674368000), UWORD(20922789888000),
  UWORD(355687428096000), UWORD(6402373705728000), UWORD(121645100408832000),
  UWORD(2432902008176640000),
#endif
};

/* xxx: don't define here */
#define NMOD_CTX_REF(ring_ctx) (((nmod_t *)((ring_ctx))))
#define NMOD_CTX(ring_ctx) (*NMOD_CTX_REF(ring_ctx))

int
gr_generic_fac_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_FMPZ)
    {
        fmpz_fac_ui(res, n);
        return GR_SUCCESS;
    }
    else if (n < FAC_TAB_SIZE)
    {
        return gr_set_ui(res, fac_tab[n], ctx);
    }
    else if (gr_ctx_is_finite_characteristic(ctx) == T_TRUE)
    {
        gr_method_binary_op_ui mul_ui = GR_BINARY_OP_UI(ctx, MUL_UI);
        int status = GR_SUCCESS;
        ulong i;

        if (ctx->which_ring == GR_CTX_NMOD)
        {
            ulong m = NMOD_CTX(ctx).n;
            /* todo: multipoint eval also for other rings */
            if (n >= m)
                return gr_zero(res, ctx);

            if (n >= 2000000)
                return gr_set_ui(res, n_factorial_fast_mod2_preinv(n, NMOD_CTX(ctx).n, NMOD_CTX(ctx).ninv), ctx);
        }

        status |= gr_set_ui(res, fac_tab[FAC_TAB_SIZE - 1], ctx);

        i = FAC_TAB_SIZE;

        for ( ; i + 8 <= FLINT_MIN(n, UWORD(1) << (FLINT_BITS / 8)); i += 8)
            status |= mul_ui(res, res, i*(i+1)*(i+2)*(i+3)*(i+4)*(i+5)*(i+6)*(i+7), ctx);
        for ( ; i + 6 <= FLINT_MIN(n, UWORD(1) << (FLINT_BITS / 6)); i += 6)
            status |= mul_ui(res, res, i*(i+1)*(i+2)*(i+3)*(i+4)*(i+5), ctx);
        for ( ; i + 4 <= FLINT_MIN(n, UWORD(1) << (FLINT_BITS / 4)); i += 4)
            status |= mul_ui(res, res, i*(i+1)*(i+2)*(i+3), ctx);
        for ( ; i + 3 <= FLINT_MIN(n, UWORD(1) << (FLINT_BITS / 3)); i += 3)
            status |= mul_ui(res, res, i*(i+1)*(i+2), ctx);
        for ( ; i + 2 <= FLINT_MIN(n, UWORD(1) << (FLINT_BITS / 2)); i += 2)
            status |= mul_ui(res, res, i*(i+1), ctx);
        for ( ; i <= n; i += 1)
            status |= mul_ui(res, res, i, ctx);

        return status;
    }
    else if (gr_ctx_has_real_prec(ctx) == T_TRUE)  /* compute via arb */
    {
        gr_ctx_t RR;
        arb_t t;
        slong prec;
        int status = GR_SUCCESS;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));
        gr_ctx_init_real_arb(RR, prec);
        arb_init(t);
        status |= gr_fac_ui(t, n, RR);
        status |= gr_set_other(res, t, RR, ctx);
        arb_clear(t);
        gr_ctx_clear(RR);
        return status;
    }
    else        /* compute via fmpz */
    {
        int status;
        fmpz_t t;
        fmpz_init(t);
        fmpz_fac_ui(t, n);
        status = gr_set_fmpz(res, t, ctx);
        fmpz_clear(t);
        return status;
    }
}

int
gr_generic_fac_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
{
    if (*n >= 0 && *n <= COEFF_MAX)
    {
        return gr_fac_ui(res, *n, ctx);
    }
    else if (fmpz_sgn(n) < 0)
    {
        return GR_DOMAIN;
    }
    else if (gr_ctx_has_real_prec(ctx) == T_TRUE)
    {
        gr_ctx_t RR;
        arb_t t;
        slong prec;
        int status = GR_SUCCESS;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));
        gr_ctx_init_real_arb(RR, prec);
        arb_init(t);
        status |= gr_fac_fmpz(t, n, RR);
        status |= gr_set_other(res, t, RR, ctx);
        arb_clear(t);
        gr_ctx_clear(RR);
        return status;
    }
    else
    {
        return GR_UNABLE;
    }
}

int
gr_generic_fac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    ulong n;

    if (gr_get_ui(&n, x, ctx) == GR_SUCCESS)
    {
        return gr_fac_ui(res, n, ctx);
    }
    else
    {
        return gr_add_ui(res, x, 1, ctx) | gr_gamma(res, res, ctx);
    }
}

int
gr_generic_fac_vec(gr_ptr res, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op_ui mul_ui = GR_BINARY_OP_UI(ctx, MUL_UI);
    gr_method_unary_op_ui set_ui = GR_UNARY_OP_UI(ctx, SET_UI);
    int status = GR_SUCCESS;
    slong i;
    slong sz = ctx->sizeof_elem;

    for (i = 0; i < FLINT_MIN(FAC_TAB_SIZE, len); i++)
        status |= set_ui(GR_ENTRY(res, i, sz), fac_tab[i], ctx);

    /* todo: consider increasing numerical precision */
    for (; i < len; i++)
        status |= mul_ui(GR_ENTRY(res, i, sz), GR_ENTRY(res, i - 1, sz), i, ctx);

    return status;
}

int
gr_generic_rfac_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
{
    return gr_fac_ui(res, n, ctx) | gr_inv(res, res, ctx);
}

int
gr_generic_rfac_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
{
    if (0 <= *n && *n <= COEFF_MAX)
    {
        return gr_rfac_ui(res, *n, ctx);
    }
    else if (fmpz_sgn(n) < 0)
    {
        return gr_zero(res, ctx);
    }
    else
    {
        fmpz_t t;
        int status;
        fmpz_init(t);
        fmpz_add_ui(t, n, 1);
        status = gr_set_fmpz(res, t, ctx) | gr_rgamma(res, res, ctx);
        fmpz_clear(t);
        return status;
    }
}

int
gr_generic_rfac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    ulong n;

    if (gr_get_ui(&n, x, ctx) == GR_SUCCESS)
    {
        return gr_rfac_ui(res, n, ctx);
    }
    else
    {
        return gr_add_ui(res, x, 1, ctx) | gr_rgamma(res, res, ctx);
    }
}

int
gr_generic_rfac_vec(gr_ptr res, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op_ui mul_ui = GR_BINARY_OP_UI(ctx, MUL_UI);
    int status = GR_SUCCESS;
    slong i;
    slong sz = ctx->sizeof_elem;

    /* todo: increase numerical precision */
    if (len >= 3)
    {
        status = gr_rfac_ui(GR_ENTRY(res, len - 1, sz), len - 1, ctx);

        if (status == GR_SUCCESS)
        {
            for (i = len - 2; i >= 2; i--)
                status |= mul_ui(GR_ENTRY(res, i, sz), GR_ENTRY(res, i + 1, sz), i + 1, ctx);
        }
    }

    if (len >= 2) status |= gr_one(GR_ENTRY(res, 1, sz), ctx);
    if (len >= 1) status |= gr_one(res, ctx);

    /* make 1/2 exact */
    if (len >= 3 && gr_ctx_has_real_prec(ctx) == T_TRUE)
        status |= gr_mul_2exp_si(GR_ENTRY(res, 2, sz), res, -1, ctx);

    return status;
}

/* todo: specializations (fmpz) */
/* todo: finite characteristic */
int
gr_generic_doublefac_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (n % 2 == 0)
    {
        status |= gr_fac_ui(res, n / 2, ctx);
        status |= gr_mul_2exp_si(res, res, n / 2, ctx);
    }
    else
    {
        gr_ptr t;
        GR_TMP_INIT(t, ctx);
        status |= gr_doublefac_ui(t, n - 1, ctx);
        status |= gr_fac_ui(res, n, ctx);
        status |= gr_div(res, res, t, ctx);
        GR_TMP_CLEAR(t, ctx);

        if (status != GR_SUCCESS)
            status = GR_UNABLE;
    }

    return status;
}

int
gr_generic_doublefac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ptr t, u, v;

    ulong n;

    if (gr_get_ui(&n, x, ctx) == GR_SUCCESS)
        return gr_doublefac_ui(res, n, ctx);

    /* 2^(x/2) (pi/2)^((cospi(x)-1)/4) gamma(1+x/2) */

    GR_TMP_INIT3(t, u, v, ctx);

    status |= gr_cos_pi(t, x, ctx);
    status |= gr_sub_ui(t, t, 1, ctx);
    status |= gr_mul_2exp_si(t, t, -2, ctx);
    status |= gr_pi(u, ctx);
    status |= gr_mul_2exp_si(u, u, -1, ctx);
    status |= gr_pow(u, u, t, ctx);

    status |= gr_mul_2exp_si(t, x, -1, ctx);
    status |= gr_set_ui(v, 2, ctx);
    status |= gr_pow(t, v, t, ctx);
    status |= gr_mul(t, t, u, ctx);

    status |= gr_mul_2exp_si(u, x, -1, ctx);
    status |= gr_add_ui(u, u, 1, ctx);
    status |= gr_gamma(u, u, ctx);

    status |= gr_mul(res, t, u, ctx);

    GR_TMP_CLEAR3(t, u, v, ctx);

    if (status != GR_SUCCESS)
        status = GR_UNABLE;

    return status;
}

/* todo: any non-stupid algorithm; at least use rising2 if nothing else */
/* todo: specializations */
int
gr_generic_harmonic_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
{
    gr_ptr s, t;
    ulong i;
    int status = GR_SUCCESS;

    if (n <= 1)
        return gr_set_ui(res, n, ctx);

    if (n > 100 && gr_ctx_has_real_prec(ctx) == T_TRUE)
    {
        gr_ptr t;
        GR_TMP_INIT(t, ctx);
        status |= gr_set_ui(t, n, ctx);
        status |= gr_add_ui(t, t, 1, ctx);
        status |= gr_digamma(t, t, ctx);
        status |= gr_euler(res, ctx);
        status |= gr_add(res, res, t, ctx);
        GR_TMP_CLEAR(t, ctx);
        return status;
    }

    if (n <= 100 || gr_ctx_is_finite_characteristic(ctx) == T_FALSE)
    {
        fmpq_t t;
        fmpq_init(t);
        fmpq_harmonic_ui(t, n);
        status = gr_set_fmpq(res, t, ctx);
        fmpq_clear(t);
        return status;
    }

    GR_TMP_INIT2(s, t, ctx);

    for (i = n; i >= 1 && status == GR_SUCCESS; i--)
    {
        status |= gr_set_ui(t, i, ctx);
        status |= gr_inv(t, t, ctx);
        status |= gr_add(s, s, t, ctx);
    }

    gr_swap(res, s, ctx);

    GR_TMP_CLEAR2(s, t, ctx);

    return status;
}

int
gr_generic_harmonic(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    ulong n;

    if (gr_get_ui(&n, x, ctx) == GR_SUCCESS)
    {
        return gr_harmonic_ui(res, n, ctx);
    }
    else
    {
        gr_ptr t;
        int status = GR_SUCCESS;
        GR_TMP_INIT(t, ctx);
        status |= gr_add_ui(t, x, 1, ctx);
        status |= gr_digamma(t, t, ctx);
        status |= gr_euler(res, ctx);
        status |= gr_add(res, res, t, ctx);
        GR_TMP_CLEAR(t, ctx);

        if (status != GR_SUCCESS)
            status = GR_UNABLE;

        return status;
    }
}

/* todo */
int
gr_generic_beta(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    gr_ptr t, u, v;
    int status = GR_SUCCESS;

    GR_TMP_INIT3(t, u, v, ctx);
    status |= gr_gamma(t, x, ctx);
    status |= gr_gamma(u, y, ctx);
    status |= gr_add(v, x, y, ctx);
    status |= gr_rgamma(v, v, ctx);
    status |= gr_mul(res, t, u, ctx);
    status |= gr_mul(res, res, v, ctx);
    GR_TMP_CLEAR3(t, u, v, ctx);

    if (status != GR_SUCCESS)
        status = GR_UNABLE;

    return status;
}
