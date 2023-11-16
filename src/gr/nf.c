/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "fexpr.h"
#include "fmpz.h"
#include "nf.h"
#include "nf_elem.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_vec.h"
#include "gr_poly.h"

typedef struct
{
    nf_struct * nf;
    char * var;
}
_gr_nf_ctx_t;

#define NF_CTX(ctx) (((_gr_nf_ctx_t *)(ctx))->nf)
#define NF_VAR(ctx) (((_gr_nf_ctx_t *)(ctx))->var)

static const char * default_var = "a";

int
_gr_nf_ctx_write(gr_stream_t out, const gr_ctx_t ctx)
{
    gr_stream_write(out, "Number field ");
    gr_stream_write_free(out, fmpq_poly_get_str_pretty(NF_CTX(ctx)->pol, NF_VAR(ctx)));
    return GR_SUCCESS;
}

int _gr_nf_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
{
    slong len;
    len = strlen(s);

    if (NF_VAR(ctx) == default_var)
        NF_VAR(ctx) = NULL;

    NF_VAR(ctx) = flint_realloc(NF_VAR(ctx), len + 1);
    memcpy(NF_VAR(ctx), s, len + 1);
    return GR_SUCCESS;
}

void
_gr_nf_ctx_clear(gr_ctx_t ctx)
{
    nf_clear(NF_CTX(ctx));
    flint_free(NF_CTX(ctx));
    if (NF_VAR(ctx) != default_var)
        flint_free(NF_VAR(ctx));
}

void
_gr_nf_init(nf_elem_t x, const gr_ctx_t ctx)
{
    nf_elem_init(x, NF_CTX(ctx));
}

void
_gr_nf_clear(nf_elem_t x, const gr_ctx_t ctx)
{
    nf_elem_clear(x, NF_CTX(ctx));
}

void
_gr_nf_swap(nf_elem_t x, nf_elem_t y, const gr_ctx_t ctx)
{
    nf_elem_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

void
_gr_nf_set_shallow(nf_elem_t res, const nf_elem_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

int
_gr_nf_randtest(nf_elem_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    switch (n_randint(state, 10))
    {
        case 0:
            nf_elem_randtest(res, state, 100, NF_CTX(ctx));
            break;
        default:
            nf_elem_randtest(res, state, 10, NF_CTX(ctx));
    }

    return GR_SUCCESS;
}

int
_gr_nf_write(gr_stream_t out, const nf_elem_t x, const gr_ctx_t ctx)
{
    gr_stream_write_free(out, nf_elem_get_str_pretty(x, NF_VAR(ctx), NF_CTX(ctx)));
    return GR_SUCCESS;
}

int
_gr_nf_gen(nf_elem_t x, const gr_ctx_t ctx)
{
    nf_elem_gen(x, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_zero(nf_elem_t x, const gr_ctx_t ctx)
{
    nf_elem_zero(x, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_one(nf_elem_t x, const gr_ctx_t ctx)
{
    nf_elem_one(x, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_set_si(nf_elem_t res, slong v, const gr_ctx_t ctx)
{
    nf_elem_set_si(res, v, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_set_ui(nf_elem_t res, ulong v, const gr_ctx_t ctx)
{
    nf_elem_set_ui(res, v, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_set_fmpz(nf_elem_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    nf_elem_set_fmpz(res, v, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_set_fmpq(nf_elem_t res, const fmpq_t v, const gr_ctx_t ctx)
{
    nf_elem_set_fmpq(res, v, NF_CTX(ctx));
    return GR_SUCCESS;
}

int gr_generic_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t xctx, gr_ctx_t ctx);

int
_gr_nf_set_other(nf_elem_t res, gr_ptr v, gr_ctx_t v_ctx, gr_ctx_t ctx)
{
    if (v_ctx->which_ring == GR_CTX_NF && ctx->which_ring == GR_CTX_NF)
    {
        if (fmpq_poly_equal(NF_CTX(v_ctx)->pol, NF_CTX(ctx)->pol))
        {
            nf_elem_set(res, v, NF_CTX(ctx));
            return GR_SUCCESS;
        }
/*
        else if (nf_elem_is_rational(v, NF_CTX(v_ctx)))
        {
            fmpq_t t;
            fmpq_init(t);
            nf_elem_get_fmpq(t, v, NF_CTX(v_ctx));
            nf_elem_set_fmpq(res, t, NF_CTX(ctx));
            fmpq_clear(t);
        }
*/
        else
        {
            return GR_UNABLE;
        }
    }

    return gr_generic_set_other(res, v, v_ctx, ctx);
}

void fexpr_set_nf_elem(fexpr_t res, const nf_elem_t a, const nf_t nf, const fexpr_t var);

int
_gr_nf_get_fexpr(fexpr_t res, const nf_elem_t a, const gr_ctx_t ctx)
{
    fexpr_t var;
    fexpr_init(var);
    fexpr_set_symbol_str(var, NF_VAR(ctx));
    fexpr_set_nf_elem(res, a, NF_CTX(ctx), var);
    fexpr_clear(var);
    return GR_SUCCESS;
}

int
_gr_nf_set_fexpr(nf_elem_t res, fexpr_vec_t inp, gr_vec_t out, const fexpr_t expr, gr_ctx_t ctx)
{
    fexpr_t var;
    nf_elem_t gen;
    int status;

    fexpr_init(var);

    fexpr_set_symbol_str(var, NF_VAR(ctx));
    nf_elem_init(gen, NF_CTX(ctx));
    nf_elem_gen(gen, NF_CTX(ctx));

    /* todo: pop after use? */
    fexpr_vec_append(inp, var);
    GR_MUST_SUCCEED(gr_vec_append(out, gen, ctx));

    fexpr_clear(var);
    nf_elem_clear(gen, NF_CTX(ctx));

    status = gr_generic_set_fexpr(res, inp, out, expr, ctx);

    return status;
}

truth_t
_gr_nf_is_zero(const nf_elem_t x, const gr_ctx_t ctx)
{
    return nf_elem_is_zero(x, NF_CTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_nf_is_one(const nf_elem_t x, const gr_ctx_t ctx)
{
    return nf_elem_is_one(x, NF_CTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_nf_is_neg_one(const nf_elem_t x, const gr_ctx_t ctx)
{
    return nf_elem_equal_si(x, -1, NF_CTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_nf_equal(const nf_elem_t x, const nf_elem_t y, const gr_ctx_t ctx)
{
    return nf_elem_equal(x, y, NF_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_nf_set(nf_elem_t res, const nf_elem_t x, const gr_ctx_t ctx)
{
    nf_elem_set(res, x, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_neg(nf_elem_t res, const nf_elem_t x, const gr_ctx_t ctx)
{
    nf_elem_neg(res, x, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_add(nf_elem_t res, const nf_elem_t x, const nf_elem_t y, const gr_ctx_t ctx)
{
    nf_elem_add(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_add_si(nf_elem_t res, const nf_elem_t x, slong y, const gr_ctx_t ctx)
{
    nf_elem_add_si(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

/*
int
_gr_nf_add_ui(nf_elem_t res, const nf_elem_t x, ulong y, const gr_ctx_t ctx)
{
    nf_elem_add_ui(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}
*/

int
_gr_nf_add_fmpz(nf_elem_t res, const nf_elem_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    nf_elem_add_fmpz(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_add_fmpq(nf_elem_t res, const nf_elem_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    nf_elem_add_fmpq(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_sub(nf_elem_t res, const nf_elem_t x, const nf_elem_t y, const gr_ctx_t ctx)
{
    nf_elem_sub(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_sub_si(nf_elem_t res, const nf_elem_t x, slong y, const gr_ctx_t ctx)
{
    nf_elem_sub_si(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

/*
int
_gr_nf_sub_ui(nf_elem_t res, const nf_elem_t x, ulong y, const gr_ctx_t ctx)
{
    nf_elem_sub_ui(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}
*/

int
_gr_nf_sub_fmpz(nf_elem_t res, const nf_elem_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    nf_elem_sub_fmpz(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_sub_fmpq(nf_elem_t res, const nf_elem_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    nf_elem_sub_fmpq(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_mul(nf_elem_t res, const nf_elem_t x, const nf_elem_t y, const gr_ctx_t ctx)
{
    nf_elem_mul(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_mul_si(nf_elem_t res, const nf_elem_t x, slong y, const gr_ctx_t ctx)
{
    nf_elem_scalar_mul_si(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

/*
int
_gr_nf_mul_ui(nf_elem_t res, const nf_elem_t x, ulong y, const gr_ctx_t ctx)
{
    nf_elem_scalar_mul_ui(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}
*/

int
_gr_nf_mul_fmpz(nf_elem_t res, const nf_elem_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    nf_elem_scalar_mul_fmpz(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_mul_fmpq(nf_elem_t res, const nf_elem_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    nf_elem_scalar_mul_fmpq(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_mul_two(nf_elem_t res, const nf_elem_t x, const gr_ctx_t ctx)
{
    nf_elem_add(res, x, x, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_sqr(nf_elem_t res, const nf_elem_t x, const gr_ctx_t ctx)
{
    nf_elem_mul(res, x, x, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_inv(nf_elem_t res, const nf_elem_t x, const gr_ctx_t ctx)
{
    if (nf_elem_is_zero(x, NF_CTX(ctx)))
        return GR_DOMAIN;

    nf_elem_inv(res, x, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_div(nf_elem_t res, const nf_elem_t x, const nf_elem_t y, const gr_ctx_t ctx)
{
    if (nf_elem_is_zero(y, NF_CTX(ctx)))
        return GR_DOMAIN;
    nf_elem_div(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_div_si(nf_elem_t res, const nf_elem_t x, slong y, const gr_ctx_t ctx)
{
    if (y == 0)
        return GR_DOMAIN;
    nf_elem_scalar_div_si(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

/*
int
_gr_nf_div_ui(nf_elem_t res, const nf_elem_t x, ulong y, const gr_ctx_t ctx)
{
    if (y == 0)
        return GR_DOMAIN;
    nf_elem_scalar_div_ui(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}
*/

int
_gr_nf_div_fmpz(nf_elem_t res, const nf_elem_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
        return GR_DOMAIN;
    nf_elem_scalar_div_fmpz(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_div_fmpq(nf_elem_t res, const nf_elem_t x, const fmpq_t y, const gr_ctx_t ctx)
{
    if (fmpq_is_zero(y))
        return GR_DOMAIN;
    nf_elem_scalar_div_fmpq(res, x, y, NF_CTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_nf_is_invertible(const nf_elem_t x, const gr_ctx_t ctx)
{
    return nf_elem_is_zero(x, NF_CTX(ctx)) ? T_FALSE : T_TRUE;
}

int
_gr_nf_pow_ui(nf_elem_t res, const nf_elem_t x, ulong exp, const gr_ctx_t ctx)
{
    nf_elem_pow(res, x, exp, NF_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_nf_numerator(nf_elem_t res, const nf_elem_t x, const gr_ctx_t ctx)
{
    if (NF_CTX(ctx)->flag & NF_LINEAR)
    {
        fmpz_set(LNF_ELEM_NUMREF(res), LNF_ELEM_NUMREF(x));
        fmpz_one(LNF_ELEM_DENREF(res));
    }
    else if (NF_CTX(ctx)->flag & NF_QUADRATIC)
    {
        fmpz_set(QNF_ELEM_NUMREF(res), QNF_ELEM_NUMREF(x));
        fmpz_set(QNF_ELEM_NUMREF(res) + 1, QNF_ELEM_NUMREF(x) + 1);
        fmpz_one(QNF_ELEM_DENREF(res));
    }
    else
    {
        fmpq_poly_set(NF_ELEM(res), NF_ELEM(x));
        fmpz_one(NF_ELEM_DENREF(res));
    }

    return GR_SUCCESS;
}

int
_gr_nf_denominator(nf_elem_t res, const nf_elem_t x, const gr_ctx_t ctx)
{
    if (NF_CTX(ctx)->flag & NF_LINEAR)
    {
        fmpz_set(LNF_ELEM_NUMREF(res), LNF_ELEM_DENREF(x));
        fmpz_one(LNF_ELEM_DENREF(res));
    }
    else if (NF_CTX(ctx)->flag & NF_QUADRATIC)
    {
        fmpz_set(QNF_ELEM_NUMREF(res), QNF_ELEM_DENREF(x));
        fmpz_zero(QNF_ELEM_NUMREF(res) + 1);
        fmpz_one(QNF_ELEM_DENREF(res));
    }
    else
    {
        fmpq_poly_set_fmpz(NF_ELEM(res), NF_ELEM_DENREF(x));
        fmpz_one(NF_ELEM_DENREF(res));
    }

    return GR_SUCCESS;
}

/* todo: dot products, without intermediate reductions? */
/* todo: polynomial multiplication, etc. */


int _nf_methods_initialized = 0;

gr_static_method_table _nf_methods;

gr_method_tab_input _nf_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_nf_ctx_write},
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_nf_ctx_clear},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_GEN_NAME, (gr_funcptr) _gr_nf_ctx_set_gen_name},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_nf_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_nf_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_nf_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_nf_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_nf_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_nf_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_nf_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_nf_one},
    {GR_METHOD_GEN,             (gr_funcptr) _gr_nf_gen},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_nf_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_nf_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_nf_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_nf_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_nf_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_nf_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_nf_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_nf_set_fmpz},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_nf_set_other},

    {GR_METHOD_SET_FEXPR,       (gr_funcptr) _gr_nf_set_fexpr},
    {GR_METHOD_GET_FEXPR,       (gr_funcptr) _gr_nf_get_fexpr},


    {GR_METHOD_NEG,             (gr_funcptr) _gr_nf_neg},

    {GR_METHOD_ADD,             (gr_funcptr) _gr_nf_add},
/*    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_nf_add_ui}, */
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_nf_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_nf_add_fmpz},
    {GR_METHOD_ADD_FMPQ,        (gr_funcptr) _gr_nf_add_fmpq},

    {GR_METHOD_SUB,             (gr_funcptr) _gr_nf_sub},
/*    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_nf_sub_ui}, */
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_nf_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_nf_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,        (gr_funcptr) _gr_nf_sub_fmpq},

    {GR_METHOD_MUL,             (gr_funcptr) _gr_nf_mul},
/*    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_nf_mul_ui}, */
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_nf_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_nf_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,        (gr_funcptr) _gr_nf_mul_fmpq},

    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_nf_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_nf_sqr},

    {GR_METHOD_DIV,             (gr_funcptr) _gr_nf_div},
/*    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_nf_div_ui}, */
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_nf_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_nf_div_fmpz},
    {GR_METHOD_DIV_FMPQ,        (gr_funcptr) _gr_nf_div_fmpq},

    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_nf_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_nf_inv},

    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_nf_pow_ui},

    {GR_METHOD_NUMERATOR,       (gr_funcptr) _gr_nf_numerator},
    {GR_METHOD_DENOMINATOR,     (gr_funcptr) _gr_nf_denominator},

    {0,                         (gr_funcptr) NULL},
};

/* todo: verify irreducibility? */
void
gr_ctx_init_nf(gr_ctx_t ctx, const fmpq_poly_t poly)
{
    ctx->which_ring = GR_CTX_NF;
    ctx->sizeof_elem = sizeof(nf_elem_struct);
    ctx->size_limit = WORD_MAX;

    NF_CTX(ctx) = flint_malloc(sizeof(nf_struct));
    nf_init(NF_CTX(ctx), poly);
    NF_VAR(ctx) = (char *) default_var;

    ctx->methods = _nf_methods;

    if (!_nf_methods_initialized)
    {
        gr_method_tab_init(_nf_methods, _nf_methods_input);
        _nf_methods_initialized = 1;
    }
}

void
gr_ctx_init_nf_fmpz_poly(gr_ctx_t ctx, const fmpz_poly_t poly)
{
    fmpq_poly_t f;
    fmpz one = 1;
    f->coeffs = poly->coeffs;
    f->alloc = poly->alloc;
    f->length = poly->length;
    *f->den = one;

    gr_ctx_init_nf(ctx, f);
}

void
_gr_ctx_init_nf_from_ref(gr_ctx_t ctx, const void * nfctx)
{
    ctx->which_ring = GR_CTX_NF;
    ctx->sizeof_elem = sizeof(nf_elem_struct);
    ctx->size_limit = WORD_MAX;

    NF_CTX(ctx) = (nf_struct *) nfctx;
    NF_VAR(ctx) = (char *) default_var;

    ctx->methods = _nf_methods;

    if (!_nf_methods_initialized)
    {
        gr_method_tab_init(_nf_methods, _nf_methods_input);
        _nf_methods_initialized = 1;
    }
}
