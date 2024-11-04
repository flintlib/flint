/*
    Copyright (C) 2023, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "fmpq.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_generic.h"

#ifdef __GNUC__
# define strcmp __builtin_strcmp
#else
# include <string.h>
#endif

static const char * default_var = "x";

static void _gr_gr_series_mod_ctx_clear(gr_ctx_t ctx)
{
    if (SERIES_MOD_CTX(ctx)->var != default_var)
        flint_free(SERIES_MOD_CTX(ctx)->var);
}

static int _gr_gr_series_mod_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Power series over ");
    gr_ctx_write(out, SERIES_MOD_ELEM_CTX(ctx));
    gr_stream_write(out, " mod ");
    gr_stream_write(out, SERIES_MOD_CTX(ctx)->var);
    gr_stream_write(out, "^");
    gr_stream_write_si(out, SERIES_MOD_N(ctx));
    return GR_SUCCESS;
}

static truth_t
_gr_gr_series_mod_ctx_is_ring(gr_ctx_t ctx)
{
    if (SERIES_MOD_N(ctx) == 0)
        return T_TRUE;

    return gr_ctx_is_ring(SERIES_MOD_ELEM_CTX(ctx));
}

static truth_t
_gr_gr_series_mod_ctx_is_commutative_ring(gr_ctx_t ctx)
{
    if (SERIES_MOD_N(ctx) == 0)
        return T_TRUE;

    return gr_ctx_is_commutative_ring(SERIES_MOD_ELEM_CTX(ctx));
}

static truth_t
_gr_gr_series_mod_ctx_is_integral_domain(gr_ctx_t ctx)
{
    if (SERIES_MOD_N(ctx) != 1)
        return T_FALSE;

    return gr_ctx_is_integral_domain(SERIES_MOD_ELEM_CTX(ctx));
}

static truth_t
_gr_gr_series_mod_ctx_is_field(gr_ctx_t ctx)
{
    if (SERIES_MOD_N(ctx) != 1)
        return T_FALSE;

    return gr_ctx_is_field(SERIES_MOD_ELEM_CTX(ctx));
}

static int _gr_gr_series_mod_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
{
    slong len;
    len = strlen(s);

    if (SERIES_MOD_CTX(ctx)->var == default_var)
        SERIES_MOD_CTX(ctx)->var = NULL;

    SERIES_MOD_CTX(ctx)->var = flint_realloc(SERIES_MOD_CTX(ctx)->var, len + 1);
    memcpy(SERIES_MOD_CTX(ctx)->var, s, len + 1);
    return GR_SUCCESS;
}

static int _gr_gr_series_mod_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
{
    return _gr_gr_series_mod_ctx_set_gen_name(ctx, s[0]);
}

static int
_gr_gr_series_mod_gens_recursive(gr_vec_t vec, gr_ctx_t ctx)
{
    int status;
    gr_vec_t vec1;
    slong i, n;

    /* Get generators of the element ring */
    gr_vec_init(vec1, 0, SERIES_MOD_ELEM_CTX(ctx));
    status = gr_gens_recursive(vec1, SERIES_MOD_ELEM_CTX(ctx));
    n = vec1->length;

    gr_vec_set_length(vec, n + 1, ctx);

    /* Promote to polynomials */
    for (i = 0; i < n; i++)
        status |= gr_poly_set_scalar(gr_vec_entry_ptr(vec, i, ctx),
                gr_vec_entry_srcptr(vec1, i, SERIES_MOD_ELEM_CTX(ctx)),
                SERIES_MOD_ELEM_CTX(ctx));

    status |= gr_poly_gen(gr_vec_entry_ptr(vec, n, ctx), SERIES_MOD_ELEM_CTX(ctx));

    gr_vec_clear(vec1, SERIES_MOD_ELEM_CTX(ctx));

    return status;
}

static void _gr_gr_series_mod_init(gr_poly_t res, gr_ctx_t ctx)
{
    gr_poly_init(res, SERIES_MOD_ELEM_CTX(ctx));
}

static void _gr_gr_series_mod_clear(gr_poly_t res, gr_ctx_t ctx)
{
    gr_poly_clear(res, SERIES_MOD_ELEM_CTX(ctx));
}

static void _gr_gr_series_mod_swap(gr_poly_t x, gr_poly_t y, gr_ctx_t ctx)
{
    gr_poly_swap(x, y, SERIES_MOD_ELEM_CTX(ctx));
}

static int _gr_gr_series_mod_randtest(gr_poly_t res, flint_rand_t state, gr_ctx_t ctx)
{
    return gr_poly_randtest(res, state, SERIES_MOD_N(ctx), SERIES_MOD_ELEM_CTX(ctx));
}

static int _gr_gr_series_mod_write(gr_stream_t out, const gr_poly_t x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_poly_write(out, x, SERIES_MOD_CTX(ctx)->var, SERIES_MOD_ELEM_CTX(ctx));
    gr_stream_write(out, " (mod ");
    gr_stream_write(out, SERIES_MOD_CTX(ctx)->var);
    gr_stream_write(out, "^");
    gr_stream_write_si(out, SERIES_MOD_N(ctx));
    gr_stream_write(out, ")");
    return status;
}

static int _gr_gr_series_mod_zero(gr_poly_t res, gr_ctx_t ctx)
{
    return gr_poly_zero(res, SERIES_MOD_ELEM_CTX(ctx));
}

static int _gr_gr_series_mod_one(gr_poly_t res, gr_ctx_t ctx)
{
    return (SERIES_MOD_N(ctx) == 0) ? GR_SUCCESS : gr_poly_one(res, SERIES_MOD_ELEM_CTX(ctx));
}

static int _gr_gr_series_mod_gen(gr_poly_t res, gr_ctx_t ctx)
{
    return (SERIES_MOD_N(ctx) <= 1) ? gr_poly_zero(res, ctx) : gr_poly_gen(res, SERIES_MOD_ELEM_CTX(ctx));
}

static int _gr_gr_series_mod_set(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
{
    return gr_poly_set(res, x, SERIES_MOD_ELEM_CTX(ctx));
}

static void _gr_gr_series_mod_set_shallow(gr_poly_t res, const gr_poly_t x, gr_ctx_t FLINT_UNUSED(ctx))
{
    *res = *x;
}

static int _gr_gr_series_mod_set_si(gr_poly_t res, slong c, gr_ctx_t ctx)
{
    return (SERIES_MOD_N(ctx) == 0) ? GR_SUCCESS : gr_poly_set_si(res, c, SERIES_MOD_ELEM_CTX(ctx));
}

static int _gr_gr_series_mod_set_ui(gr_poly_t res, ulong c, gr_ctx_t ctx)
{
    return (SERIES_MOD_N(ctx) == 0) ? GR_SUCCESS : gr_poly_set_ui(res, c, SERIES_MOD_ELEM_CTX(ctx));
}

static int _gr_gr_series_mod_set_fmpz(gr_poly_t res, const fmpz_t c, gr_ctx_t ctx)
{
    return (SERIES_MOD_N(ctx) == 0) ? GR_SUCCESS : gr_poly_set_fmpz(res, c, SERIES_MOD_ELEM_CTX(ctx));
}

static int _gr_gr_series_mod_set_fmpq(gr_poly_t res, const fmpq_t c, gr_ctx_t ctx)
{
    return (SERIES_MOD_N(ctx) == 0) ? GR_SUCCESS : gr_poly_set_fmpq(res, c, SERIES_MOD_ELEM_CTX(ctx));
}

static truth_t _gr_gr_series_mod_is_zero(const gr_poly_t x, gr_ctx_t ctx)
{
    return (SERIES_MOD_N(ctx) == 0) ? T_TRUE : gr_poly_is_zero(x, SERIES_MOD_ELEM_CTX(ctx));
}

static truth_t _gr_gr_series_mod_is_one(const gr_poly_t x, gr_ctx_t ctx)
{
    return (SERIES_MOD_N(ctx) == 0) ? T_TRUE : gr_poly_is_one(x, SERIES_MOD_ELEM_CTX(ctx));
}

static truth_t _gr_gr_series_mod_equal(const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
{
    return gr_poly_equal(x, y, SERIES_MOD_ELEM_CTX(ctx));
}

static int _gr_gr_series_mod_neg(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
{
    return gr_poly_neg(res, x, SERIES_MOD_ELEM_CTX(ctx));
}

static int _gr_gr_series_mod_add(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
{
    return gr_poly_add(res, x, y, SERIES_MOD_ELEM_CTX(ctx));
}

static int _gr_gr_series_mod_sub(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
{
    return gr_poly_sub(res, x, y, SERIES_MOD_ELEM_CTX(ctx));
}

static int _gr_gr_series_mod_mul(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
{
    return gr_poly_mullow(res, x, y, SERIES_MOD_N(ctx), SERIES_MOD_ELEM_CTX(ctx));
}

static int _gr_gr_series_mod_inv(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
{
    return gr_poly_inv_series(res, x, SERIES_MOD_N(ctx), SERIES_MOD_ELEM_CTX(ctx));
}

static int _gr_gr_series_mod_div(gr_poly_t res, const gr_poly_t x, const gr_poly_t y, gr_ctx_t ctx)
{
    return gr_poly_div_series(res, x, y, SERIES_MOD_N(ctx), SERIES_MOD_ELEM_CTX(ctx));
}

#define UNARY_POLY_WRAPPER(func) \
static int \
_gr_gr_series_mod_ ## func(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx) \
{ \
    return gr_poly_ ## func ## _series(res, x, SERIES_MOD_N(ctx), SERIES_MOD_ELEM_CTX(ctx)); \
} \

UNARY_POLY_WRAPPER(exp)
UNARY_POLY_WRAPPER(log)
UNARY_POLY_WRAPPER(rsqrt)
UNARY_POLY_WRAPPER(tan)
UNARY_POLY_WRAPPER(asin)
UNARY_POLY_WRAPPER(acos)
UNARY_POLY_WRAPPER(atan)
UNARY_POLY_WRAPPER(asinh)
UNARY_POLY_WRAPPER(acosh)
UNARY_POLY_WRAPPER(atanh)


/* fixme: gr_poly_sqrt_series does not deal with leading zeros */
static int _gr_gr_series_mod_sqrt(gr_poly_t res, const gr_poly_t x, gr_ctx_t ctx)
{
    slong n = SERIES_MOD_N(ctx);

    if (n == 0)
        return GR_SUCCESS;

    if (x->length == 0)
        return gr_poly_zero(res, SERIES_MOD_ELEM_CTX(ctx));

    if (gr_is_zero(x->coeffs, SERIES_MOD_ELEM_CTX(ctx)) != T_FALSE)
        return GR_UNABLE;

    return gr_poly_sqrt_series(res, x, n, SERIES_MOD_ELEM_CTX(ctx));
}


typedef struct
{
    gr_poly_struct poly;
    slong error;
}
gr_series_struct;

static int
_set_truncate_poly(gr_poly_t res, const gr_poly_t x, gr_ctx_t x_elem_ctx, slong n, gr_ctx_t elem_ctx)
{
    if (x_elem_ctx == elem_ctx)
    {
        return gr_poly_truncate(res, x, n, elem_ctx);
    }
    else
    {
        if (gr_poly_length(x, x_elem_ctx) <= n)
        {
            return gr_poly_set_gr_poly_other(res, x, x_elem_ctx, elem_ctx);
        }
        else
        {
            int status = GR_SUCCESS;
            gr_poly_t t;
            gr_poly_init(t, x_elem_ctx);
            status |= gr_poly_truncate(t, x, n, x_elem_ctx);
            status |= gr_poly_set_gr_poly_other(res, x, x_elem_ctx, elem_ctx);
            gr_poly_clear(t, x_elem_ctx);
            return status;
        }
    }
}

static int
_gr_gr_series_mod_set_other(gr_poly_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (x_ctx == ctx)
    {
        return gr_poly_set(res, x, ctx);
    }
    else if (x_ctx == SERIES_MOD_ELEM_CTX(ctx))
    {
        return (SERIES_MOD_N(ctx) == 0) ? GR_SUCCESS : gr_poly_set_scalar(res, x, SERIES_MOD_ELEM_CTX(ctx));
    }
    else if (x_ctx->which_ring == GR_CTX_SERIES_MOD_GR_POLY && !strcmp(SERIES_MOD_CTX(x_ctx)->var, SERIES_MOD_CTX(ctx)->var))
    {
        if (SERIES_MOD_N(ctx) <= SERIES_MOD_N(x_ctx))   /* set() does modular reduction */
            return _set_truncate_poly(res, x, SERIES_MOD_ELEM_CTX(x_ctx), SERIES_MOD_N(ctx), SERIES_MOD_ELEM_CTX(ctx));
        else
            return GR_DOMAIN;   /* set() doesn't do lift */
    }
    else if (x_ctx->which_ring == GR_CTX_GR_POLY && !strcmp(POLYNOMIAL_CTX(x_ctx)->var, SERIES_MOD_CTX(ctx)->var))
    {
        return _set_truncate_poly(res, x, POLYNOMIAL_ELEM_CTX(x_ctx), SERIES_MOD_N(ctx), SERIES_MOD_ELEM_CTX(ctx));
    }
    else if (x_ctx->which_ring == GR_CTX_GR_SERIES && !strcmp(SERIES_CTX(x_ctx)->var, SERIES_MOD_CTX(ctx)->var))
    {
        /* all coefficients below x^n must be known */
        if (((const gr_series_struct *) x)->error < SERIES_MOD_N(ctx))
            return GR_UNABLE;
        else
            return _set_truncate_poly(res, &((const gr_series_struct *) x)->poly, SERIES_ELEM_CTX(x_ctx), SERIES_MOD_N(ctx), SERIES_MOD_ELEM_CTX(ctx));
    }
    else
    {
        int status = GR_SUCCESS;

        gr_poly_fit_length(res, 1, SERIES_MOD_ELEM_CTX(ctx));
        status = gr_set_other(res->coeffs, x, x_ctx, SERIES_MOD_ELEM_CTX(ctx));

        if (status == GR_SUCCESS)
        {
            _gr_poly_set_length(res, 1, SERIES_MOD_ELEM_CTX(ctx));
            _gr_poly_normalise(res, SERIES_MOD_ELEM_CTX(ctx));
        }
        else
            _gr_poly_set_length(res, 0, SERIES_MOD_ELEM_CTX(ctx));

        status |= gr_poly_truncate(res, res, SERIES_MOD_N(ctx), SERIES_MOD_ELEM_CTX(ctx));
        return status;
    }

    return GR_UNABLE;
}



int _gr_series_mod_methods_initialized = 0;

gr_static_method_table _gr_series_mod_methods;

gr_method_tab_input _gr_series_mod_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) _gr_gr_series_mod_ctx_clear},
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) _gr_gr_series_mod_ctx_write},
    {GR_METHOD_CTX_SET_GEN_NAME, (gr_funcptr) _gr_gr_series_mod_ctx_set_gen_name},
    {GR_METHOD_CTX_SET_GEN_NAMES, (gr_funcptr) _gr_gr_series_mod_ctx_set_gen_names},
    {GR_METHOD_CTX_IS_RING, (gr_funcptr) _gr_gr_series_mod_ctx_is_ring},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) _gr_gr_series_mod_ctx_is_commutative_ring},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN, (gr_funcptr) _gr_gr_series_mod_ctx_is_integral_domain},
    {GR_METHOD_CTX_IS_FIELD, (gr_funcptr) _gr_gr_series_mod_ctx_is_field},
    {GR_METHOD_INIT,        (gr_funcptr) _gr_gr_series_mod_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) _gr_gr_series_mod_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) _gr_gr_series_mod_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) _gr_gr_series_mod_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) _gr_gr_series_mod_randtest},
    {GR_METHOD_WRITE,       (gr_funcptr) _gr_gr_series_mod_write},
    {GR_METHOD_ZERO,        (gr_funcptr) _gr_gr_series_mod_zero},
    {GR_METHOD_ONE,         (gr_funcptr) _gr_gr_series_mod_one},
    {GR_METHOD_IS_ZERO,     (gr_funcptr) _gr_gr_series_mod_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) _gr_gr_series_mod_is_one},
    {GR_METHOD_EQUAL,       (gr_funcptr) _gr_gr_series_mod_equal},
    {GR_METHOD_GEN,         (gr_funcptr) _gr_gr_series_mod_gen},
    {GR_METHOD_GENS,        (gr_funcptr) gr_generic_gens_single},
    {GR_METHOD_GENS_RECURSIVE,  (gr_funcptr) _gr_gr_series_mod_gens_recursive},
    {GR_METHOD_SET,         (gr_funcptr) _gr_gr_series_mod_set},
    {GR_METHOD_SET_UI,      (gr_funcptr) _gr_gr_series_mod_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) _gr_gr_series_mod_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) _gr_gr_series_mod_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) _gr_gr_series_mod_set_fmpq},
    {GR_METHOD_SET_OTHER,   (gr_funcptr) _gr_gr_series_mod_set_other},
    {GR_METHOD_SET_STR,     (gr_funcptr) gr_generic_set_str_balance_additions},
    {GR_METHOD_NEG,         (gr_funcptr) _gr_gr_series_mod_neg},
    {GR_METHOD_ADD,         (gr_funcptr) _gr_gr_series_mod_add},
    {GR_METHOD_SUB,         (gr_funcptr) _gr_gr_series_mod_sub},
    {GR_METHOD_MUL,         (gr_funcptr) _gr_gr_series_mod_mul},
    {GR_METHOD_INV,         (gr_funcptr) _gr_gr_series_mod_inv},
    {GR_METHOD_DIV,         (gr_funcptr) _gr_gr_series_mod_div},
    {GR_METHOD_SQRT,        (gr_funcptr) _gr_gr_series_mod_sqrt},
    {GR_METHOD_RSQRT,       (gr_funcptr) _gr_gr_series_mod_rsqrt},
    {GR_METHOD_EXP,         (gr_funcptr) _gr_gr_series_mod_exp},
    {GR_METHOD_LOG,         (gr_funcptr) _gr_gr_series_mod_log},
    {GR_METHOD_TAN,         (gr_funcptr) _gr_gr_series_mod_tan},
    {GR_METHOD_ASIN,        (gr_funcptr) _gr_gr_series_mod_asin},
    {GR_METHOD_ACOS,        (gr_funcptr) _gr_gr_series_mod_acos},
    {GR_METHOD_ATAN,        (gr_funcptr) _gr_gr_series_mod_atan},
    {GR_METHOD_ASINH,       (gr_funcptr) _gr_gr_series_mod_asinh},
    {GR_METHOD_ACOSH,       (gr_funcptr) _gr_gr_series_mod_acosh},
    {GR_METHOD_ATANH,       (gr_funcptr) _gr_gr_series_mod_atanh},
/*
    {GR_METHOD_GAMMA,       (gr_funcptr) _gr_gr_series_mod_gamma},
    {GR_METHOD_RGAMMA,      (gr_funcptr) _gr_gr_series_mod_rgamma},
    {GR_METHOD_LGAMMA,      (gr_funcptr) _gr_gr_series_mod_lgamma},
    {GR_METHOD_DIGAMMA,     (gr_funcptr) _gr_gr_series_mod_digamma},
    {GR_METHOD_ERF,         (gr_funcptr) _gr_gr_series_mod_erf},
    {GR_METHOD_ERFC,        (gr_funcptr) _gr_gr_series_mod_erfc},
    {GR_METHOD_ERFI,        (gr_funcptr) _gr_gr_series_mod_erfi},
    {GR_METHOD_FRESNEL,     (gr_funcptr) _gr_gr_series_mod_fresnel},
    {GR_METHOD_FRESNEL_S,   (gr_funcptr) _gr_gr_series_mod_fresnel_s},
    {GR_METHOD_FRESNEL_C,   (gr_funcptr) _gr_gr_series_mod_fresnel_c},
    {GR_METHOD_AIRY,        (gr_funcptr) _gr_gr_series_mod_airy},
    {GR_METHOD_AIRY_AI,        (gr_funcptr) _gr_gr_series_mod_airy_ai},
    {GR_METHOD_AIRY_AI_PRIME,  (gr_funcptr) _gr_gr_series_mod_airy_ai_prime},
    {GR_METHOD_AIRY_BI,        (gr_funcptr) _gr_gr_series_mod_airy_bi},
    {GR_METHOD_AIRY_BI_PRIME,  (gr_funcptr) _gr_gr_series_mod_airy_bi_prime},
    {GR_METHOD_EXP_INTEGRAL_EI,       (gr_funcptr) _gr_gr_series_mod_exp_integral_ei},
    {GR_METHOD_COS_INTEGRAL,          (gr_funcptr) _gr_gr_series_mod_cos_integral},
    {GR_METHOD_COSH_INTEGRAL,         (gr_funcptr) _gr_gr_series_mod_cosh_integral},
    {GR_METHOD_SIN_INTEGRAL,          (gr_funcptr) _gr_gr_series_mod_sin_integral},
    {GR_METHOD_SINH_INTEGRAL,         (gr_funcptr) _gr_gr_series_mod_sinh_integral},
    {GR_METHOD_LOG_INTEGRAL,          (gr_funcptr) _gr_gr_series_mod_log_integral},
    {GR_METHOD_GAMMA_UPPER,           (gr_funcptr) _gr_gr_series_mod_gamma_upper},
    {GR_METHOD_GAMMA_LOWER,           (gr_funcptr) _gr_gr_series_mod_gamma_lower},
    {GR_METHOD_BETA_LOWER,            (gr_funcptr) _gr_gr_series_mod_beta_lower},
    {GR_METHOD_HYPGEOM_PFQ,           (gr_funcptr) _gr_gr_series_mod_hypgeom_pfq},
    {GR_METHOD_HURWITZ_ZETA,          (gr_funcptr) _gr_gr_series_mod_hurwitz_zeta},
    {GR_METHOD_POLYLOG,               (gr_funcptr) _gr_gr_series_mod_polylog},
    {GR_METHOD_DIRICHLET_L,           (gr_funcptr) _gr_gr_series_mod_dirichlet_l},
    {GR_METHOD_DIRICHLET_HARDY_Z,     (gr_funcptr) _gr_gr_series_mod_dirichlet_hardy_z},
    {GR_METHOD_DIRICHLET_HARDY_THETA, (gr_funcptr) _gr_gr_series_mod_dirichlet_hardy_theta},
    {GR_METHOD_JACOBI_THETA,          (gr_funcptr) _gr_gr_series_mod_jacobi_theta},
    {GR_METHOD_JACOBI_THETA_1,          (gr_funcptr) _gr_gr_series_mod_jacobi_theta_1},
    {GR_METHOD_JACOBI_THETA_2,          (gr_funcptr) _gr_gr_series_mod_jacobi_theta_2},
    {GR_METHOD_JACOBI_THETA_3,          (gr_funcptr) _gr_gr_series_mod_jacobi_theta_3},
    {GR_METHOD_JACOBI_THETA_4,          (gr_funcptr) _gr_gr_series_mod_jacobi_theta_4},
    {GR_METHOD_AGM1,                   (gr_funcptr) _gr_gr_series_mod_agm1},
    {GR_METHOD_ELLIPTIC_K,             (gr_funcptr) _gr_gr_series_mod_elliptic_k},
    {GR_METHOD_WEIERSTRASS_P,          (gr_funcptr) _gr_gr_series_mod_weierstrass_p},
*/
    {0,                     (gr_funcptr) NULL},
};

void
gr_ctx_init_series_mod_gr_poly(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)
{
    ctx->which_ring = GR_CTX_SERIES_MOD_GR_POLY;
    ctx->sizeof_elem = sizeof(gr_poly_struct);
    ctx->size_limit = WORD_MAX;

    SERIES_MOD_CTX(ctx)->base_ring = (gr_ctx_struct *) base_ring;
    SERIES_MOD_CTX(ctx)->var = (char *) default_var;
    SERIES_MOD_CTX(ctx)->n = FLINT_MAX(0, n);

    ctx->methods = _gr_series_mod_methods;

    if (!_gr_series_mod_methods_initialized)
    {
        gr_method_tab_init(_gr_series_mod_methods, _gr_series_mod_methods_input);
        _gr_series_mod_methods_initialized = 1;
    }
}
