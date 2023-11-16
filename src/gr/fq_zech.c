/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fq_nmod.h"
#include "fq_zech.h"
#include "fq_zech_vec.h"
#include "fq_zech_poly.h"
#include "fq_zech_poly_factor.h"
#include "fq_zech_mat.h"
#include "gr.h"
#include "gr_vec.h"

#define FQ_CTX(ring_ctx) ((fq_zech_ctx_struct *)(GR_CTX_DATA_AS_PTR(ring_ctx)))

/* todo: lots of inlining */

void
_gr_fq_zech_ctx_clear(gr_ctx_t ctx)
{
    fq_zech_ctx_clear(FQ_CTX(ctx));
    flint_free(GR_CTX_DATA_AS_PTR(ctx));
}

int
_gr_fq_zech_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Finite field (fq_zech)");
    return GR_SUCCESS;
}

void
_gr_fq_zech_init(fq_zech_t x, const gr_ctx_t ctx)
{
    fq_zech_init(x, FQ_CTX(ctx));
}

void
_gr_fq_zech_clear(fq_zech_t x, const gr_ctx_t ctx)
{
    fq_zech_clear(x, FQ_CTX(ctx));
}

void
_gr_fq_zech_swap(fq_zech_t x, fq_zech_t y, const gr_ctx_t ctx)
{
    fq_zech_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

void
_gr_fq_zech_set_shallow(fq_zech_t res, const fq_zech_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

int
_gr_fq_zech_randtest(fq_zech_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    fq_zech_randtest(res, state, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_write(gr_stream_t out, const fq_zech_t x, const gr_ctx_t ctx)
{
    gr_stream_write_free(out, fq_zech_get_str_pretty(x, FQ_CTX(ctx)));
    return GR_SUCCESS;
}

int
_gr_fq_zech_zero(fq_zech_t x, const gr_ctx_t ctx)
{
    fq_zech_zero(x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_one(fq_zech_t x, const gr_ctx_t ctx)
{
    fq_zech_one(x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_set_si(fq_zech_t res, slong v, const gr_ctx_t ctx)
{
    fq_zech_set_si(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_set_ui(fq_zech_t res, ulong v, const gr_ctx_t ctx)
{
    fq_zech_set_ui(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_set_fmpz(fq_zech_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    fq_zech_set_fmpz(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_fq_zech_is_zero(const fq_zech_t x, const gr_ctx_t ctx)
{
    return fq_zech_is_zero(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fq_zech_is_one(const fq_zech_t x, const gr_ctx_t ctx)
{
    return fq_zech_is_one(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fq_zech_equal(const fq_zech_t x, const fq_zech_t y, const gr_ctx_t ctx)
{
    return fq_zech_equal(x, y, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fq_zech_set(fq_zech_t res, const fq_zech_t x, const gr_ctx_t ctx)
{
    fq_zech_set(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_neg(fq_zech_t res, const fq_zech_t x, const gr_ctx_t ctx)
{
    fq_zech_neg(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_add(fq_zech_t res, const fq_zech_t x, const fq_zech_t y, const gr_ctx_t ctx)
{
    fq_zech_add(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_sub(fq_zech_t res, const fq_zech_t x, const fq_zech_t y, const gr_ctx_t ctx)
{
    fq_zech_sub(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_mul(fq_zech_t res, const fq_zech_t x, const fq_zech_t y, const gr_ctx_t ctx)
{
    fq_zech_mul(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_mul_si(fq_zech_t res, const fq_zech_t x, slong y, const gr_ctx_t ctx)
{
    fq_zech_mul_si(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_mul_ui(fq_zech_t res, const fq_zech_t x, ulong y, const gr_ctx_t ctx)
{
    fq_zech_mul_ui(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_mul_fmpz(fq_zech_t res, const fq_zech_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fq_zech_mul_fmpz(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_inv(fq_zech_t res, const fq_zech_t x, const gr_ctx_t ctx)
{
    if (fq_zech_is_zero(x, FQ_CTX(ctx)))
    {
        return GR_DOMAIN;
    }
    else
    {
        fq_zech_inv(res, x, FQ_CTX(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_fq_zech_div(fq_zech_t res, const fq_zech_t x, const fq_zech_t y, const gr_ctx_t ctx)
{
    if (fq_zech_is_zero(y, FQ_CTX(ctx)))
    {
        return GR_DOMAIN;
    }
    else
    {
        fq_zech_t t;
        fq_zech_init(t, FQ_CTX(ctx));
        fq_zech_inv(t, y, FQ_CTX(ctx));
        fq_zech_mul(res, x, t, FQ_CTX(ctx));
        fq_zech_clear(t, FQ_CTX(ctx));
        return GR_SUCCESS;
    }
}


truth_t
_gr_fq_zech_is_invertible(const fq_zech_t x, const gr_ctx_t ctx)
{
    return (!fq_zech_is_zero(x, FQ_CTX(ctx))) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fq_zech_is_square(const fq_zech_t x, const gr_ctx_t ctx)
{
    return fq_zech_is_square(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fq_zech_sqrt(fq_zech_t res, const fq_zech_t x, const gr_ctx_t ctx)
{
    if (fq_zech_sqrt(res, x, FQ_CTX(ctx)))
    {
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_ctx_fq_zech_prime(fmpz_t p, gr_ctx_t ctx)
{
    fmpz_set(p, fq_zech_ctx_prime(FQ_CTX(ctx)));
    return GR_SUCCESS;
}

int
_gr_ctx_fq_zech_degree(slong * deg, gr_ctx_t ctx)
{
    *deg = fq_zech_ctx_degree(FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ctx_fq_zech_order(fmpz_t q, gr_ctx_t ctx)
{
    fq_zech_ctx_order(q, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_gen(gr_ptr res, gr_ctx_t ctx)
{
    fq_zech_gen(res, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_frobenius(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx)
{
    fq_zech_frobenius(res, x, e, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_multiplicative_order(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    int ret;
    ret = fq_zech_multiplicative_order(res, x, FQ_CTX(ctx));

    if (ret == 1)
        return GR_SUCCESS;

    /* todo: better solution? */
    return GR_DOMAIN;
}

int
_gr_fq_zech_norm(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_zech_norm(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_trace(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_zech_trace(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_fq_zech_is_primitive(gr_srcptr x, gr_ctx_t ctx)
{
    return fq_zech_is_primitive(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fq_zech_pth_root(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_zech_pth_root(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

void
_gr_fq_zech_vec_init(fq_zech_struct * vec, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        fq_zech_init(vec + i, FQ_CTX(ctx));
}

void
_gr_fq_zech_vec_clear(fq_zech_struct * vec, slong len, gr_ctx_t ctx)
{
}

void
_gr_fq_zech_vec_swap(fq_zech_struct * vec1, fq_zech_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        fq_zech_swap(vec1 + i, vec2 + i, FQ_CTX(ctx));
}

int
_gr_fq_zech_vec_set(fq_zech_struct * res, const fq_zech_struct * vec, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        res[i] = vec[i];

    return GR_SUCCESS;
}

int
_gr_fq_zech_vec_normalise(slong * res, const fq_zech_struct * vec, slong len, gr_ctx_t ctx)
{
    while (len > 0 && fq_zech_is_zero(vec + len - 1, FQ_CTX(ctx)))
        len--;

    res[0] = len;
    return GR_SUCCESS;
}

slong
_gr_fq_zech_vec_normalise_weak(const fq_zech_struct * vec, slong len, gr_ctx_t ctx)
{
    while (len > 0 && fq_zech_is_zero(vec + len - 1, FQ_CTX(ctx)))
        len--;

    return len;
}

int
_gr_fq_zech_vec_mul_scalar(fq_zech_struct * res, const fq_zech_struct * vec, slong len, const fq_zech_t x, gr_ctx_t ctx)
{
    _fq_zech_vec_scalar_mul_fq_zech(res, vec, len, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_vec_addmul_scalar(fq_zech_struct * res, const fq_zech_struct * vec, slong len, const fq_zech_t x, gr_ctx_t ctx)
{
    /* the poly method checks for special cases. worth it? */
    _fq_zech_poly_scalar_addmul_fq_zech(res, vec, len, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_zech_vec_submul_scalar(fq_zech_struct * res, const fq_zech_struct * vec, slong len, const fq_zech_t x, gr_ctx_t ctx)
{
    /* the poly method checks for special cases. worth it? */
    _fq_zech_poly_scalar_submul_fq_zech(res, vec, len, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}


/* todo: _fq_zech_poly_mullow should do the right thing */
int
_gr_fq_zech_poly_mullow(fq_zech_struct * res,
    const fq_zech_struct * poly1, slong len1,
    const fq_zech_struct * poly2, slong len2, slong n, gr_ctx_t ctx)
{
    if (len1 + len2 - 1 == n)
    {
        if (poly1 == poly2 && len1 == len2)
            _fq_zech_poly_sqr(res, poly1, len1, FQ_CTX(ctx));
        else if (len1 >= len2)
            _fq_zech_poly_mul(res, poly1, len1, poly2, len2, FQ_CTX(ctx));
        else
            _fq_zech_poly_mul(res, poly2, len2, poly1, len1, FQ_CTX(ctx));
    }
    else
    {
        if (len1 >= len2)
            _fq_zech_poly_mullow(res, poly1, len1, poly2, len2, n, FQ_CTX(ctx));
        else
            _fq_zech_poly_mullow(res, poly2, len2, poly1, len1, n, FQ_CTX(ctx));
    }

    return GR_SUCCESS;
}

/* todo: also need the _other version ... ? */
/* todo: implement generically */

int
_gr_fq_zech_roots_gr_poly(gr_vec_t roots, gr_vec_t mult, const fq_zech_poly_t poly, int flags, gr_ctx_t ctx)
{
    if (poly->length == 0)
        return GR_DOMAIN;

    {
        gr_ctx_t ZZ;
        fq_zech_poly_factor_t fac;
        slong i, num;

        gr_ctx_init_fmpz(ZZ);
        fq_zech_poly_factor_init(fac, FQ_CTX(ctx));
        fq_zech_poly_roots(fac, poly, 1, FQ_CTX(ctx));

        num = fac->num;

        gr_vec_set_length(roots, num, ctx);
        gr_vec_set_length(mult, num, ZZ);

        for (i = 0; i < num; i++)
        {
            fq_zech_neg(gr_vec_entry_ptr(roots, i, ctx), fac->poly[i].coeffs, FQ_CTX(ctx));

            /* work around flint bug: factors can be non-monic */
            if (!fq_zech_is_one(fac->poly[i].coeffs + 1, FQ_CTX(ctx)))
                fq_zech_div(gr_vec_entry_ptr(roots, i, ctx), gr_vec_entry_ptr(roots, i, ctx), fac->poly[i].coeffs + 1, FQ_CTX(ctx));

            fmpz_set_ui(((fmpz *) mult->entries) + i, fac->exp[i]);
        }

        fq_zech_poly_factor_clear(fac, FQ_CTX(ctx));
        gr_ctx_clear(ZZ);
    }

    return GR_SUCCESS;
}

int
_gr_fq_zech_mat_mul(fq_zech_mat_t res, const fq_zech_mat_t x, const fq_zech_mat_t y, gr_ctx_t ctx)
{
    fq_zech_mat_mul(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int _fq_zech_methods_initialized = 0;

gr_static_method_table _fq_zech_methods;

gr_method_tab_input _fq_zech_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_fq_zech_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_fq_zech_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_fq_zech_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_fq_zech_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_fq_zech_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_fq_zech_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_fq_zech_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_fq_zech_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_fq_zech_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_fq_zech_one},
    {GR_METHOD_GEN,             (gr_funcptr) _gr_fq_zech_gen},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_fq_zech_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_fq_zech_is_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_fq_zech_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_fq_zech_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_fq_zech_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_fq_zech_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_fq_zech_set_fmpz},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_fq_zech_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_fq_zech_add},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_fq_zech_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_fq_zech_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_fq_zech_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_fq_zech_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_fq_zech_mul_fmpz},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fq_zech_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fq_zech_inv},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_fq_zech_div},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_fq_zech_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_fq_zech_sqrt},

    {GR_METHOD_CTX_FQ_PRIME,            (gr_funcptr) _gr_ctx_fq_zech_prime},
    {GR_METHOD_CTX_FQ_DEGREE,           (gr_funcptr) _gr_ctx_fq_zech_degree},
    {GR_METHOD_CTX_FQ_ORDER,            (gr_funcptr) _gr_ctx_fq_zech_order},
    {GR_METHOD_FQ_FROBENIUS,            (gr_funcptr) _gr_fq_zech_frobenius},
    {GR_METHOD_FQ_MULTIPLICATIVE_ORDER, (gr_funcptr) _gr_fq_zech_multiplicative_order},
    {GR_METHOD_FQ_NORM,                 (gr_funcptr) _gr_fq_zech_norm},
    {GR_METHOD_FQ_TRACE,                (gr_funcptr) _gr_fq_zech_trace},
    {GR_METHOD_FQ_IS_PRIMITIVE,         (gr_funcptr) _gr_fq_zech_is_primitive},
    {GR_METHOD_FQ_PTH_ROOT,             (gr_funcptr) _gr_fq_zech_pth_root},

    {GR_METHOD_VEC_INIT,            (gr_funcptr) _gr_fq_zech_vec_init},
    {GR_METHOD_VEC_CLEAR,           (gr_funcptr) _gr_fq_zech_vec_clear},
    {GR_METHOD_VEC_SET,             (gr_funcptr) _gr_fq_zech_vec_set},
    {GR_METHOD_VEC_SWAP,            (gr_funcptr) _gr_fq_zech_vec_swap},
    {GR_METHOD_VEC_NORMALISE,       (gr_funcptr) _gr_fq_zech_vec_normalise},
    {GR_METHOD_VEC_NORMALISE_WEAK,  (gr_funcptr) _gr_fq_zech_vec_normalise_weak},
    {GR_METHOD_VEC_MUL_SCALAR,            (gr_funcptr) _gr_fq_zech_vec_mul_scalar},
    {GR_METHOD_VEC_ADDMUL_SCALAR,            (gr_funcptr) _gr_fq_zech_vec_addmul_scalar},
    {GR_METHOD_VEC_SUBMUL_SCALAR,            (gr_funcptr) _gr_fq_zech_vec_submul_scalar},

    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_fq_zech_poly_mullow},
    {GR_METHOD_POLY_ROOTS,      (gr_funcptr) _gr_fq_zech_roots_gr_poly},

    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_fq_zech_mat_mul},
    {0,                         (gr_funcptr) NULL},
};

void
_gr_ctx_init_fq_zech_from_ref(gr_ctx_t ctx, const void * fq_zech_ctx)
{
    ctx->which_ring = GR_CTX_FQ_ZECH;
    ctx->sizeof_elem = sizeof(fq_zech_struct);
    GR_CTX_DATA_AS_PTR(ctx) = (fq_zech_ctx_struct *) fq_zech_ctx;
    ctx->size_limit = WORD_MAX;
    ctx->methods = _fq_zech_methods;

    if (!_fq_zech_methods_initialized)
    {
        gr_method_tab_init(_fq_zech_methods, _fq_zech_methods_input);
        _fq_zech_methods_initialized = 1;
    }
}

void
gr_ctx_init_fq_zech(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var)
{
    fq_zech_ctx_struct * fq_zech_ctx;

    fq_zech_ctx = flint_malloc(sizeof(fq_zech_ctx_struct));
    fq_zech_ctx_init(fq_zech_ctx, p, d, var == NULL ? "a" : var);
    _gr_ctx_init_fq_zech_from_ref(ctx, fq_zech_ctx);
}
