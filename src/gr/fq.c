/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "gr.h"
#include "fq.h"
#include "fq_poly.h"
#include "fq_mat.h"
#include "fq_poly_factor.h"
#include "gr_vec.h"

#define FQ_CTX(ring_ctx) ((fq_ctx_struct *)(GR_CTX_DATA_AS_PTR(ring_ctx)))

void
_gr_fq_ctx_clear(gr_ctx_t ctx)
{
    fq_ctx_clear(FQ_CTX(ctx));
    flint_free(GR_CTX_DATA_AS_PTR(ctx));
}

int
_gr_fq_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Finite field (fq)");
    return GR_SUCCESS;
}

void
_gr_fq_init(fq_t x, const gr_ctx_t ctx)
{
    fq_init(x, FQ_CTX(ctx));
}

void
_gr_fq_clear(fq_t x, const gr_ctx_t ctx)
{
    fq_clear(x, FQ_CTX(ctx));
}

void
_gr_fq_swap(fq_t x, fq_t y, const gr_ctx_t ctx)
{
    fq_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

void
_gr_fq_set_shallow(fq_t res, const fq_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

int
_gr_fq_randtest(fq_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    fq_randtest(res, state, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_write(gr_stream_t out, const fq_t x, const gr_ctx_t ctx)
{
    gr_stream_write_free(out, fq_get_str_pretty(x, FQ_CTX(ctx)));
    return GR_SUCCESS;
}

int
_gr_fq_zero(fq_t x, const gr_ctx_t ctx)
{
    fq_zero(x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_one(fq_t x, const gr_ctx_t ctx)
{
    fq_one(x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_set_si(fq_t res, slong v, const gr_ctx_t ctx)
{
    fq_set_si(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_set_ui(fq_t res, ulong v, const gr_ctx_t ctx)
{
    fq_set_ui(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_set_fmpz(fq_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    fq_set_fmpz(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_set_other(fq_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            fq_set_fmpz(res, x, FQ_CTX(ctx));
            return GR_SUCCESS;
    }

    return GR_UNABLE;
}

truth_t
_gr_fq_is_zero(const fq_t x, const gr_ctx_t ctx)
{
    return fq_is_zero(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fq_is_one(const fq_t x, const gr_ctx_t ctx)
{
    return fq_is_one(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fq_equal(const fq_t x, const fq_t y, const gr_ctx_t ctx)
{
    return fq_equal(x, y, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fq_set(fq_t res, const fq_t x, const gr_ctx_t ctx)
{
    fq_set(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_neg(fq_t res, const fq_t x, const gr_ctx_t ctx)
{
    fq_neg(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_add(fq_t res, const fq_t x, const fq_t y, const gr_ctx_t ctx)
{
    fq_add(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_sub(fq_t res, const fq_t x, const fq_t y, const gr_ctx_t ctx)
{
    fq_sub(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_mul(fq_t res, const fq_t x, const fq_t y, const gr_ctx_t ctx)
{
    fq_mul(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_mul_si(fq_t res, const fq_t x, slong y, const gr_ctx_t ctx)
{
    fq_mul_si(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_mul_ui(fq_t res, const fq_t x, ulong y, const gr_ctx_t ctx)
{
    fq_mul_ui(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_mul_fmpz(fq_t res, const fq_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fq_mul_fmpz(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_si_mul(fq_t res, slong y, const fq_t x, const gr_ctx_t ctx)
{
    fq_mul_si(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_ui_mul(fq_t res, ulong y, const fq_t x, const gr_ctx_t ctx)
{
    fq_mul_ui(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_fmpz_mul(fq_t res, const fmpz_t y, const fq_t x, const gr_ctx_t ctx)
{
    fq_mul_fmpz(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_inv(fq_t res, const fq_t x, const gr_ctx_t ctx)
{
    if (fq_is_zero(x, FQ_CTX(ctx)))
    {
        return GR_DOMAIN;
    }
    else
    {
        fq_inv(res, x, FQ_CTX(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_fq_div(fq_t res, const fq_t x, const fq_t y, const gr_ctx_t ctx)
{
    if (fq_is_zero(y, FQ_CTX(ctx)))
    {
        return GR_DOMAIN;
    }
    else
    {
        fq_t t;
        fq_init(t, FQ_CTX(ctx));
        fq_inv(t, y, FQ_CTX(ctx));
        fq_mul(res, x, t, FQ_CTX(ctx));
        fq_clear(t, FQ_CTX(ctx));
        return GR_SUCCESS;
    }
}


truth_t
_gr_fq_is_invertible(const fq_t x, const gr_ctx_t ctx)
{
    return (!fq_is_zero(x, FQ_CTX(ctx))) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fq_is_square(const fq_t x, const gr_ctx_t ctx)
{
    return fq_is_square(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fq_sqrt(fq_t res, const fq_t x, const gr_ctx_t ctx)
{
    if (fq_sqrt(res, x, FQ_CTX(ctx)))
    {
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_ctx_fq_prime(fmpz_t p, gr_ctx_t ctx)
{
    fmpz_set(p, fq_ctx_prime(FQ_CTX(ctx)));
    return GR_SUCCESS;
}

int
_gr_ctx_fq_degree(slong * deg, gr_ctx_t ctx)
{
    *deg = fq_ctx_degree(FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ctx_fq_order(fmpz_t q, gr_ctx_t ctx)
{
    fq_ctx_order(q, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_gen(gr_ptr res, gr_ctx_t ctx)
{
    fq_gen(res, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_frobenius(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx)
{
    fq_frobenius(res, x, e, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_multiplicative_order(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    int ret;
    ret = fq_multiplicative_order(res, x, FQ_CTX(ctx));

    if (ret == 1)
        return GR_SUCCESS;

    /* todo: better solution? */
    return GR_DOMAIN;
}

int
_gr_fq_norm(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_norm(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_trace(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_trace(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_fq_is_primitive(gr_srcptr x, gr_ctx_t ctx)
{
    return fq_is_primitive(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fq_pth_root(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_pth_root(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_vec_dot(fq_struct * res, const fq_struct * initial, int subtract, const fq_struct * vec1, const fq_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    fmpz * s, * t;
    slong slen, tlen, len1, len2;
    slong plen;

    if (len <= 0)
    {
        if (initial == NULL)
            fq_zero(res, FQ_CTX(ctx));
        else
            fq_set(res, initial, FQ_CTX(ctx));
        return GR_SUCCESS;
    }

    plen = FQ_CTX(ctx)->modulus->length;

    t = GR_TMP_ALLOC((4 * plen) * sizeof(fmpz));
    s = t + 2 * plen;

    for (i = 0; i < 4 * plen; i++)
        fmpz_init(t + i);

    len1 = vec1[0].length;
    len2 = vec2[0].length;

    if (len1 == 0 || len2 == 0)
    {
        slen = 0;
    }
    else
    {
        slen = len1 + len2 - 1;
        if (len1 >= len2)
            _fmpz_poly_mul(s, vec1[0].coeffs, len1, vec2[0].coeffs, len2);
        else
            _fmpz_poly_mul(s, vec2[0].coeffs, len2, vec1[0].coeffs, len1);
    }

    for (i = 1; i < len; i++)
    {
        len1 = vec1[i].length;
        len2 = vec2[i].length;

        if (len1 != 0 && len2 != 0)
        {
            tlen = len1 + len2 - 1;
            if (len1 >= len2)
                _fmpz_poly_mul(t, vec1[i].coeffs, len1, vec2[i].coeffs, len2);
            else
                _fmpz_poly_mul(t, vec2[i].coeffs, len2, vec1[i].coeffs, len1);

            _fmpz_poly_add(s, s, slen, t, tlen);
            slen = FLINT_MAX(slen, tlen);
        }
    }

    if (initial == NULL)
    {
        if (subtract)
            _fmpz_vec_neg(s, s, slen);
    }
    else
    {
        len2 = initial->length;

        if (subtract)
            _fmpz_poly_sub(s, initial->coeffs, len2, s, slen);
        else
            _fmpz_poly_add(s, initial->coeffs, len2, s, slen);

        slen = FLINT_MAX(slen, len2);
    }

    while (slen > 0 && fmpz_is_zero(s + slen - 1))
        slen--;

    _fq_reduce(s, slen, FQ_CTX(ctx));
    slen = FLINT_MIN(slen, plen - 1);

    while (slen > 0 && fmpz_is_zero(s + slen - 1))
        slen--;

    fmpz_poly_fit_length(res, slen);
    _fmpz_vec_set(res->coeffs, s, slen); /* todo: swap */
    _fmpz_poly_set_length(res, slen);

    for (i = 0; i < 4 * plen; i++)
        fmpz_clear(t + i);

    GR_TMP_FREE(t, (4 * plen) * sizeof(fmpz));

    return GR_SUCCESS;
}

int
_gr_fq_vec_dot_rev(fq_struct * res, const fq_struct * initial, int subtract, const fq_struct * vec1, const fq_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    fmpz * s, * t;
    slong slen, tlen, len1, len2;
    slong plen;

    if (len <= 0)
    {
        if (initial == NULL)
            fq_zero(res, FQ_CTX(ctx));
        else
            fq_set(res, initial, FQ_CTX(ctx));
        return GR_SUCCESS;
    }

    plen = FQ_CTX(ctx)->modulus->length;

    t = GR_TMP_ALLOC((4 * plen) * sizeof(fmpz));
    s = t + 2 * plen;

    for (i = 0; i < 4 * plen; i++)
        fmpz_init(t + i);

    len1 = vec1[0].length;
    len2 = vec2[len - 1].length;

    if (len1 == 0 || len2 == 0)
    {
        slen = 0;
    }
    else
    {
        slen = len1 + len2 - 1;
        if (len1 >= len2)
            _fmpz_poly_mul(s, vec1[0].coeffs, len1, vec2[len - 1].coeffs, len2);
        else
            _fmpz_poly_mul(s, vec2[len - 1].coeffs, len2, vec1[0].coeffs, len1);
    }

    for (i = 1; i < len; i++)
    {
        len1 = vec1[i].length;
        len2 = vec2[len - 1 - i].length;

        if (len1 != 0 && len2 != 0)
        {
            tlen = len1 + len2 - 1;
            if (len1 >= len2)
                _fmpz_poly_mul(t, vec1[i].coeffs, len1, vec2[len - 1 - i].coeffs, len2);
            else
                _fmpz_poly_mul(t, vec2[len - 1 - i].coeffs, len2, vec1[i].coeffs, len1);

            _fmpz_poly_add(s, s, slen, t, tlen);
            slen = FLINT_MAX(slen, tlen);
        }
    }

    if (initial == NULL)
    {
        if (subtract)
            _fmpz_vec_neg(s, s, slen);
    }
    else
    {
        len2 = initial->length;

        if (subtract)
            _fmpz_poly_sub(s, initial->coeffs, len2, s, slen);
        else
            _fmpz_poly_add(s, initial->coeffs, len2, s, slen);

        slen = FLINT_MAX(slen, len2);
    }

    while (slen > 0 && fmpz_is_zero(s + slen - 1))
        slen--;

    _fq_reduce(s, slen, FQ_CTX(ctx));
    slen = FLINT_MIN(slen, plen - 1);

    while (slen > 0 && fmpz_is_zero(s + slen - 1))
        slen--;

    fmpz_poly_fit_length(res, slen);
    _fmpz_vec_set(res->coeffs, s, slen); /* todo: swap */
    _fmpz_poly_set_length(res, slen);

    for (i = 0; i < 4 * plen; i++)
        fmpz_clear(t + i);

    GR_TMP_FREE(t, (4 * plen) * sizeof(fmpz));

    return GR_SUCCESS;
}

/* todo: _fq_poly_mullow should do the right thing */
int
_gr_fq_poly_mullow(fq_struct * res,
    const fq_struct * poly1, slong len1,
    const fq_struct * poly2, slong len2, slong n, gr_ctx_t ctx)
{
    if (len1 + len2 - 1 == n)
    {
        if (poly1 == poly2 && len1 == len2)
            _fq_poly_sqr(res, poly1, len1, FQ_CTX(ctx));
        else if (len1 >= len2)
            _fq_poly_mul(res, poly1, len1, poly2, len2, FQ_CTX(ctx));
        else
            _fq_poly_mul(res, poly2, len2, poly1, len1, FQ_CTX(ctx));
    }
    else
    {
        if (len1 >= len2)
            _fq_poly_mullow(res, poly1, len1, poly2, len2, n, FQ_CTX(ctx));
        else
            _fq_poly_mullow(res, poly2, len2, poly1, len1, n, FQ_CTX(ctx));
    }

    return GR_SUCCESS;
}

/* todo: also need the _other version ... ? */
/* todo: implement generically */

int
_gr_fq_roots_gr_poly(gr_vec_t roots, gr_vec_t mult, const fq_poly_t poly, int flags, gr_ctx_t ctx)
{
    if (poly->length == 0)
        return GR_DOMAIN;

    {
        gr_ctx_t ZZ;
        fq_poly_factor_t fac;
        slong i, num;

        gr_ctx_init_fmpz(ZZ);
        fq_poly_factor_init(fac, FQ_CTX(ctx));
        fq_poly_roots(fac, poly, 1, FQ_CTX(ctx));

        num = fac->num;

        gr_vec_set_length(roots, num, ctx);
        gr_vec_set_length(mult, num, ZZ);

        for (i = 0; i < num; i++)
        {
            fq_neg(gr_vec_entry_ptr(roots, i, ctx), fac->poly[i].coeffs, FQ_CTX(ctx));

            /* work around flint bug: factors can be non-monic */
            if (!fq_is_one(fac->poly[i].coeffs + 1, FQ_CTX(ctx)))
                fq_div(gr_vec_entry_ptr(roots, i, ctx), gr_vec_entry_ptr(roots, i, ctx), fac->poly[i].coeffs + 1, FQ_CTX(ctx));

            fmpz_set_ui(((fmpz *) mult->entries) + i, fac->exp[i]);
        }

        fq_poly_factor_clear(fac, FQ_CTX(ctx));
        gr_ctx_clear(ZZ);
    }

    return GR_SUCCESS;
}

int
_gr_fq_mat_mul(fq_mat_t res, const fq_mat_t x, const fq_mat_t y, gr_ctx_t ctx)
{
    fq_mat_mul(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int _fq_methods_initialized = 0;

gr_static_method_table _fq_methods;

gr_method_tab_input _fq_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_fq_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_fq_ctx_write},
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
    {GR_METHOD_INIT,            (gr_funcptr) _gr_fq_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_fq_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_fq_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_fq_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_fq_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_fq_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_fq_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_fq_one},
    {GR_METHOD_GEN,                     (gr_funcptr) _gr_fq_gen},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_fq_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_fq_is_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_fq_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_fq_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_fq_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_fq_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_fq_set_fmpz},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_fq_set_other},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_fq_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_fq_add},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_fq_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_fq_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_fq_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_fq_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_fq_mul_fmpz},
/*
    todo ...
    {GR_METHOD_SI_MUL,          (gr_funcptr) _gr_fq_si_mul},
    {GR_METHOD_UI_MUL,          (gr_funcptr) _gr_fq_ui_mul},
    {GR_METHOD_FMPZ_MUL,        (gr_funcptr) _gr_fq_fmpz_mul},
    {GR_METHOD_MUL_OTHER,        (gr_funcptr) _gr_fq_mul_other},
    {GR_METHOD_OTHER_MUL,        (gr_funcptr) _gr_fq_other_mul},
*/

    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fq_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fq_inv},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_fq_div},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_fq_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_fq_sqrt},

    {GR_METHOD_CTX_FQ_PRIME,            (gr_funcptr) _gr_ctx_fq_prime},
    {GR_METHOD_CTX_FQ_DEGREE,           (gr_funcptr) _gr_ctx_fq_degree},
    {GR_METHOD_CTX_FQ_ORDER,            (gr_funcptr) _gr_ctx_fq_order},
    {GR_METHOD_FQ_FROBENIUS,            (gr_funcptr) _gr_fq_frobenius},
    {GR_METHOD_FQ_MULTIPLICATIVE_ORDER, (gr_funcptr) _gr_fq_multiplicative_order},
    {GR_METHOD_FQ_NORM,                 (gr_funcptr) _gr_fq_norm},
    {GR_METHOD_FQ_TRACE,                (gr_funcptr) _gr_fq_trace},
    {GR_METHOD_FQ_IS_PRIMITIVE,         (gr_funcptr) _gr_fq_is_primitive},
    {GR_METHOD_FQ_PTH_ROOT,             (gr_funcptr) _gr_fq_pth_root},

    {GR_METHOD_VEC_DOT,         (gr_funcptr) _gr_fq_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _gr_fq_vec_dot_rev},
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_fq_poly_mullow},
    {GR_METHOD_POLY_ROOTS,      (gr_funcptr) _gr_fq_roots_gr_poly},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_fq_mat_mul},
    {0,                         (gr_funcptr) NULL},
};

void
_gr_ctx_init_fq_from_ref(gr_ctx_t ctx, const void * fq_ctx)
{
    ctx->which_ring = GR_CTX_FQ;
    ctx->sizeof_elem = sizeof(fq_struct);
    GR_CTX_DATA_AS_PTR(ctx) = (fq_ctx_struct *) fq_ctx;
    ctx->size_limit = WORD_MAX;
    ctx->methods = _fq_methods;

    if (!_fq_methods_initialized)
    {
        gr_method_tab_init(_fq_methods, _fq_methods_input);
        _fq_methods_initialized = 1;
    }
}

void
gr_ctx_init_fq(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var)
{
    fq_ctx_struct * fq_ctx;
    fq_ctx = flint_malloc(sizeof(fq_ctx_struct));
    fq_ctx_init(fq_ctx, p, d, var == NULL ? "a" : var);
    _gr_ctx_init_fq_from_ref(ctx, fq_ctx);
}
