/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include "arf.h"
#include "gr.h"

typedef struct
{
    slong prec;
}
_mpf_ctx_struct;

#define GR_MPF_CTX(ctx) ((_mpf_ctx_struct *)(ctx))
#define GR_MPF_CTX_PREC(ctx) (GR_MPF_CTX(ctx)->prec)

static int
_gr_mpf_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Floating-point numbers with prec = ");
    gr_stream_write_si(out, GR_MPF_CTX_PREC(ctx));
    gr_stream_write(out, " (mpf)");
    return GR_SUCCESS;
}

static void
_gr_mpf_init(mpf_t res, gr_ctx_t ctx)
{
    mpf_init2(res, GR_MPF_CTX_PREC(ctx));
}

static void
_gr_mpf_clear(mpf_t res, gr_ctx_t ctx)
{
    mpf_clear(res);
}

static void _gr_mpf_swap(mpf_t x, mpf_t y, gr_ctx_t ctx)
{
    mpf_swap(x, y);
}

static void _gr_mpf_set_shallow(mpf_t res, const mpf_t x, gr_ctx_t ctx)
{
    *res = *x;
}

static int _gr_mpf_set(mpf_t res, const mpf_t x, gr_ctx_t ctx)
{
    mpf_set(res, x);
    return GR_SUCCESS;
}

static int _gr_mpf_ctx_set_real_prec(gr_ctx_t ctx, slong prec)
{
    GR_MPF_CTX_PREC(ctx) = prec;
    return GR_SUCCESS;
}

static int _gr_mpf_ctx_get_real_prec(slong * res, gr_ctx_t ctx)
{
    *res = GR_MPF_CTX_PREC(ctx);
    return GR_SUCCESS;
}

static void
arf_get_mpf(mpf_t res, const arf_t x)
{
    mpfr_t t;
    mpfr_init2(t, FLINT_MAX(arf_bits(x), 2));
    arf_get_mpfr(t, x, MPFR_RNDD);
    mpfr_get_f(res, t, MPFR_RNDD);
    mpfr_clear(t);
}

static void
arf_set_mpf(arf_t res, const mpf_t x)
{
    mpfr_t t;
    mpfr_init2(t, mpf_get_prec(x) + FLINT_BITS);
    mpfr_set_f(t, x, MPFR_RNDD);
    arf_set_mpfr(res, t);
    mpfr_clear(t);
}

static int
_gr_mpf_randtest(mpf_t res, flint_rand_t state, gr_ctx_t ctx)
{
    gr_ctx_t arf_ctx;
    arf_t t;
    int status;
    gr_ctx_init_real_float_arf(arf_ctx, GR_MPF_CTX_PREC(ctx));
    arf_init(t);
    status = gr_randtest(t, state, arf_ctx);
    arf_get_mpf(res, t);
    arf_clear(t);
    gr_ctx_clear(arf_ctx);
    return status;
}

static int
_gr_mpf_write(gr_stream_t out, const mpf_t x, gr_ctx_t ctx)
{
    gr_ctx_t arf_ctx;
    arf_t t;
    int status;
    gr_ctx_init_real_float_arf(arf_ctx, GR_MPF_CTX_PREC(ctx));
    arf_init(t);
    arf_set_mpf(t, x);
    status = gr_write(out, t, arf_ctx);
    arf_clear(t);
    gr_ctx_clear(arf_ctx);
    return status;
}

static int
_gr_mpf_zero(mpf_t res, gr_ctx_t ctx)
{
    mpf_set_ui(res, 0);
    return GR_SUCCESS;
}

static int
_gr_mpf_one(mpf_t res, gr_ctx_t ctx)
{
    mpf_set_ui(res, 1);
    return GR_SUCCESS;
}

static int
_gr_mpf_neg_one(mpf_t res, gr_ctx_t ctx)
{
    mpf_set_si(res, -1);
    return GR_SUCCESS;
}

static int _gr_mpf_set_ui(mpf_t res, ulong x, gr_ctx_t ctx)
{
    mpz_t t;
    t->_mp_size = t->_mp_alloc = (x != 0);
    t->_mp_d = &x;
    mpf_set_z(res, t);
    return GR_SUCCESS;
}

static int _gr_mpf_set_si(mpf_t res, slong x, gr_ctx_t ctx)
{
    mpz_t t;
    ulong d;
    if (x >= 0)
    {
        d = x;
        t->_mp_size = t->_mp_alloc = (x != 0);
    }
    else
    {
        d = -x;
        t->_mp_size = -1;
        t->_mp_alloc = 1;
        t->_mp_d = &d;
    }
    t->_mp_d = &d;
    mpf_set_z(res, t);
    return GR_SUCCESS;
}

static int _gr_mpf_set_fmpz(mpf_t res, const fmpz_t x, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*x))
        return _gr_mpf_set_si(res, *x, ctx);
    else
    {
        mpf_set_z(res, COEFF_TO_PTR(*x));
        return GR_SUCCESS;
    }
}

static int _gr_mpf_set_d(mpf_t res, double x, gr_ctx_t ctx)
{
    if (isfinite(x))
    {
        mpf_set_d(res, x);
        return GR_SUCCESS;
    }
    else
    {
        return GR_UNABLE;
    }
}

static int _gr_mpf_get_fmpz(fmpz_t res, const mpf_t x, gr_ctx_t ctx)
{
    slong exp, xn;

    if (x->_mp_size == 0)
    {
        fmpz_zero(res);
        return GR_SUCCESS;
    }

    exp = x->_mp_exp;
    xn = FLINT_ABS(x->_mp_size);

    if (exp <= 0)
        return GR_DOMAIN;

    if (exp >= xn)
    {
        fmpz_set_ui_array(res, x->_mp_d, xn);
        fmpz_mul_2exp(res, res, (exp - xn) * FLINT_BITS);
    }
    else
    {
        if (!mpn_zero_p(x->_mp_d, xn - exp))
            return GR_DOMAIN;

        fmpz_set_ui_array(res, x->_mp_d + (xn - exp), exp);
    }

    if (x->_mp_size < 0)
        fmpz_neg(res, res);

    return GR_SUCCESS;
}

static truth_t _gr_mpf_equal(const mpf_t x, const mpf_t y, gr_ctx_t ctx)
{
    return (mpf_cmp(x, y) == 0) ? T_TRUE : T_FALSE;
}

static truth_t _gr_mpf_is_zero(const mpf_t x, gr_ctx_t ctx)
{
    return (mpf_cmp_ui(x, 0) == 0) ? T_TRUE : T_FALSE;
}

static truth_t _gr_mpf_is_one(const mpf_t x, gr_ctx_t ctx)
{
    return (mpf_cmp_si(x, 1) == 0) ? T_TRUE : T_FALSE;
}

static truth_t _gr_mpf_is_neg_one(const mpf_t x, gr_ctx_t ctx)
{
    return (mpf_cmp_si(x, -1) == 0) ? T_TRUE : T_FALSE;
}

static int _gr_mpf_neg(mpf_t res, const mpf_t x, gr_ctx_t ctx)
{
    mpf_neg(res, x);
    return GR_SUCCESS;
}

static int _gr_mpf_abs(mpf_t res, const mpf_t x, gr_ctx_t ctx)
{
    mpf_abs(res, x);
    return GR_SUCCESS;
}

static int _gr_mpf_floor(mpf_t res, const mpf_t x, gr_ctx_t ctx)
{
    mpf_floor(res, x);
    return GR_SUCCESS;
}

static int _gr_mpf_ceil(mpf_t res, const mpf_t x, gr_ctx_t ctx)
{
    mpf_ceil(res, x);
    return GR_SUCCESS;
}

static int _gr_mpf_trunc(mpf_t res, const mpf_t x, gr_ctx_t ctx)
{
    mpf_trunc(res, x);
    return GR_SUCCESS;
}

static int _gr_mpf_add(mpf_t res, const mpf_t x, const mpf_t y, gr_ctx_t ctx)
{
    mpf_add(res, x, y);
    return GR_SUCCESS;
}

static int _gr_mpf_sub(mpf_t res, const mpf_t x, const mpf_t y, gr_ctx_t ctx)
{
    mpf_sub(res, x, y);
    return GR_SUCCESS;
}

static int _gr_mpf_mul(mpf_t res, const mpf_t x, const mpf_t y, gr_ctx_t ctx)
{
    mpf_mul(res, x, y);
    return GR_SUCCESS;
}

static int _gr_mpf_mul_2exp_si(mpf_t res, const mpf_t x, slong y, gr_ctx_t ctx)
{
    if (y >= 0)
        mpf_mul_2exp(res, x, y);
    else
        mpf_div_2exp(res, x, -y);
    return GR_SUCCESS;
}

static int _gr_mpf_inv(mpf_t res, const mpf_t x, gr_ctx_t ctx)
{
    mpf_ui_div(res, 1, x);
    return GR_SUCCESS;
}

static int _gr_mpf_cmp(int * res, const mpf_t x, const mpf_t y, gr_ctx_t ctx)
{
    int cmp;

    cmp = mpf_cmp(x, y);

    if (cmp > 0) cmp = 1;
    if (cmp < 0) cmp = -1;
    *res = cmp;
    return GR_SUCCESS;
}

static int _gr_mpf_sgn(mpf_t res, const mpf_t x, gr_ctx_t ctx)
{
    mpf_set_si(res, mpf_sgn(x));
    return GR_SUCCESS;
}

static int _gr_mpf_pi(mpf_t res, gr_ctx_t ctx)
{
    mpfr_t t;
    mpfr_init2(t, GR_MPF_CTX_PREC(ctx));
    mpfr_const_pi(t, MPFR_RNDN);
    mpfr_get_f(res, t, MPFR_RNDN);
    mpfr_clear(t);
    return GR_SUCCESS;
}

static int _gr_mpf_div(mpf_t res, const mpf_t x, const mpf_t y, gr_ctx_t ctx)
{
    mpf_div(res, x, y);
    return GR_SUCCESS;
}

static int _gr_mpf_get_d_2exp_si(double * res, slong * exp, const mpf_t x, gr_ctx_t ctx)
{
    long e;
    *res = mpf_get_d_2exp(&e, x);
    *exp = e;
    return GR_SUCCESS;
}



int _gr_mpf_methods_initialized = 0;


gr_static_method_table _gr_mpf_methods;

gr_method_tab_input _gr_mpf_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_mpf_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_HAS_REAL_PREC, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_REAL_PREC, (gr_funcptr) _gr_mpf_ctx_set_real_prec},
    {GR_METHOD_CTX_GET_REAL_PREC, (gr_funcptr) _gr_mpf_ctx_get_real_prec},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_mpf_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_mpf_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_mpf_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_mpf_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_mpf_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_mpf_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_mpf_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_mpf_one},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) _gr_mpf_neg_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_mpf_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_mpf_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_mpf_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_mpf_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_mpf_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_mpf_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_mpf_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_mpf_set_fmpz},
    {GR_METHOD_SET_D,           (gr_funcptr) _gr_mpf_set_d},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_mpf_get_fmpz},
    {GR_METHOD_GET_D_2EXP_SI,   (gr_funcptr) _gr_mpf_get_d_2exp_si},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_mpf_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_mpf_add},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_mpf_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_mpf_mul},
    {GR_METHOD_INV,             (gr_funcptr) _gr_mpf_inv},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_mpf_div},
    {GR_METHOD_MUL_2EXP_SI,     (gr_funcptr) _gr_mpf_mul_2exp_si},
    {GR_METHOD_CMP,             (gr_funcptr) _gr_mpf_cmp},
    {GR_METHOD_SGN,             (gr_funcptr) _gr_mpf_sgn},
    {GR_METHOD_ABS,             (gr_funcptr) _gr_mpf_abs},
    {GR_METHOD_FLOOR,           (gr_funcptr) _gr_mpf_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) _gr_mpf_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) _gr_mpf_trunc},
    {GR_METHOD_PI,              (gr_funcptr) _gr_mpf_pi},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_mpf(gr_ctx_t ctx, slong prec)
{
    ctx->which_ring = GR_CTX_MPF;
    ctx->sizeof_elem = sizeof(__mpf_struct);
    ctx->size_limit = WORD_MAX;

    GR_MPF_CTX_PREC(ctx) = FLINT_MAX(1, prec);

    ctx->methods = _gr_mpf_methods;

    if (!_gr_mpf_methods_initialized)
    {
        gr_method_tab_init(_gr_mpf_methods, _gr_mpf_methods_input);
        _gr_mpf_methods_initialized = 1;
    }
}

