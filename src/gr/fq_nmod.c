/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_poly.h"
#include "fmpz.h"
#include "fq.h"
#include "fq_nmod.h"
#include "fq_nmod_mat.h"
#include "fq_nmod_poly.h"
#include "fq_nmod_poly_factor.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"

#define FQ_CTX(ring_ctx) ((fq_nmod_ctx_struct *)(GR_CTX_DATA_AS_PTR(ring_ctx)))

void
_gr_fq_nmod_ctx_clear(gr_ctx_t ctx)
{
    fq_nmod_ctx_clear(FQ_CTX(ctx));
    flint_free(GR_CTX_DATA_AS_PTR(ctx));
}

int
_gr_fq_nmod_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Finite field (fq_nmod)");
    return GR_SUCCESS;
}

void
_gr_fq_nmod_init(fq_nmod_t x, const gr_ctx_t ctx)
{
    fq_nmod_init(x, FQ_CTX(ctx));
}

void
_gr_fq_nmod_clear(fq_nmod_t x, const gr_ctx_t ctx)
{
    fq_nmod_clear(x, FQ_CTX(ctx));
}

void
_gr_fq_nmod_swap(fq_nmod_t x, fq_nmod_t y, const gr_ctx_t ctx)
{
    fq_nmod_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

void
_gr_fq_nmod_set_shallow(fq_nmod_t res, const fq_nmod_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

int
_gr_fq_nmod_randtest(fq_nmod_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    fq_nmod_randtest(res, state, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_write(gr_stream_t out, const fq_nmod_t x, const gr_ctx_t ctx)
{
    gr_stream_write_free(out, fq_nmod_get_str_pretty(x, FQ_CTX(ctx)));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_zero(fq_nmod_t x, const gr_ctx_t ctx)
{
    fq_nmod_zero(x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_one(fq_nmod_t x, const gr_ctx_t ctx)
{
    fq_nmod_one(x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_set_si(fq_nmod_t res, slong v, const gr_ctx_t ctx)
{
    fq_nmod_set_si(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_set_ui(fq_nmod_t res, ulong v, const gr_ctx_t ctx)
{
    fq_nmod_set_ui(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_set_fmpz(fq_nmod_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    fq_nmod_set_fmpz(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_fq_nmod_is_zero(const fq_nmod_t x, const gr_ctx_t ctx)
{
    return fq_nmod_is_zero(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fq_nmod_is_one(const fq_nmod_t x, const gr_ctx_t ctx)
{
    return fq_nmod_is_one(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fq_nmod_equal(const fq_nmod_t x, const fq_nmod_t y, const gr_ctx_t ctx)
{
    return fq_nmod_equal(x, y, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fq_nmod_set(fq_nmod_t res, const fq_nmod_t x, const gr_ctx_t ctx)
{
    fq_nmod_set(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_neg(fq_nmod_t res, const fq_nmod_t x, const gr_ctx_t ctx)
{
    fq_nmod_neg(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_add(fq_nmod_t res, const fq_nmod_t x, const fq_nmod_t y, const gr_ctx_t ctx)
{
    fq_nmod_add(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_sub(fq_nmod_t res, const fq_nmod_t x, const fq_nmod_t y, const gr_ctx_t ctx)
{
    fq_nmod_sub(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_mul(fq_nmod_t res, const fq_nmod_t x, const fq_nmod_t y, const gr_ctx_t ctx)
{
    fq_nmod_mul(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_mul_si(fq_nmod_t res, const fq_nmod_t x, slong y, const gr_ctx_t ctx)
{
    fq_nmod_mul_si(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_mul_ui(fq_nmod_t res, const fq_nmod_t x, ulong y, const gr_ctx_t ctx)
{
    fq_nmod_mul_ui(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_mul_fmpz(fq_nmod_t res, const fq_nmod_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fq_nmod_mul_fmpz(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_inv(fq_nmod_t res, const fq_nmod_t x, const gr_ctx_t ctx)
{
    if (fq_nmod_is_zero(x, FQ_CTX(ctx)))
    {
        return GR_DOMAIN;
    }
    else
    {
        fq_nmod_inv(res, x, FQ_CTX(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_fq_nmod_div(fq_nmod_t res, const fq_nmod_t x, const fq_nmod_t y, const gr_ctx_t ctx)
{
    if (fq_nmod_is_zero(y, FQ_CTX(ctx)))
    {
        return GR_DOMAIN;
    }
    else
    {
        fq_nmod_t t;
        fq_nmod_init(t, FQ_CTX(ctx));
        fq_nmod_inv(t, y, FQ_CTX(ctx));
        fq_nmod_mul(res, x, t, FQ_CTX(ctx));
        fq_nmod_clear(t, FQ_CTX(ctx));
        return GR_SUCCESS;
    }
}


truth_t
_gr_fq_nmod_is_invertible(const fq_nmod_t x, const gr_ctx_t ctx)
{
    return (!fq_nmod_is_zero(x, FQ_CTX(ctx))) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fq_nmod_is_square(const fq_nmod_t x, const gr_ctx_t ctx)
{
    return fq_nmod_is_square(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fq_nmod_sqrt(fq_nmod_t res, const fq_nmod_t x, const gr_ctx_t ctx)
{
    if (fq_nmod_sqrt(res, x, FQ_CTX(ctx)))
    {
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_ctx_fq_nmod_prime(fmpz_t p, gr_ctx_t ctx)
{
    fmpz_set(p, fq_nmod_ctx_prime(FQ_CTX(ctx)));
    return GR_SUCCESS;
}

int
_gr_ctx_fq_nmod_degree(slong * deg, gr_ctx_t ctx)
{
    *deg = fq_nmod_ctx_degree(FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ctx_fq_nmod_order(fmpz_t q, gr_ctx_t ctx)
{
    fq_nmod_ctx_order(q, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_gen(gr_ptr res, gr_ctx_t ctx)
{
    fq_nmod_gen(res, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_frobenius(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx)
{
    fq_nmod_frobenius(res, x, e, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_multiplicative_order(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    int ret;
    ret = fq_nmod_multiplicative_order(res, x, FQ_CTX(ctx));

    if (ret == 1)
        return GR_SUCCESS;

    /* todo: better solution? */
    return GR_DOMAIN;
}

int
_gr_fq_nmod_norm(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_nmod_norm(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_nmod_trace(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_nmod_trace(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_fq_nmod_is_primitive(gr_srcptr x, gr_ctx_t ctx)
{
    return fq_nmod_is_primitive(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fq_nmod_pth_root(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_nmod_pth_root(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

/* todo: basecase multiplication without reductions */
int
__gr_fq_nmod_vec_dot(fq_nmod_struct * res, const fq_nmod_struct * initial, int subtract, const fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    mp_ptr s, t;
    slong slen, tlen, len1, len2;
    slong plen;
    nmod_t mod;

    if (len <= 0)
    {
        if (initial == NULL)
            fq_nmod_zero(res, FQ_CTX(ctx));
        else
            fq_nmod_set(res, initial, FQ_CTX(ctx));
        return GR_SUCCESS;
    }

    plen = FQ_CTX(ctx)->modulus->length;

    t = GR_TMP_ALLOC((4 * plen) * sizeof(mp_limb_t));
    s = t + 2 * plen;

    mod = FQ_CTX(ctx)->mod;

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
            _nmod_poly_mul(s, vec1[0].coeffs, len1, vec2[0].coeffs, len2, mod);
        else
            _nmod_poly_mul(s, vec2[0].coeffs, len2, vec1[0].coeffs, len1, mod);
    }

    for (i = 1; i < len; i++)
    {
        len1 = vec1[i].length;
        len2 = vec2[i].length;

        if (len1 != 0 && len2 != 0)
        {
            tlen = len1 + len2 - 1;
            if (len1 >= len2)
                _nmod_poly_mul(t, vec1[i].coeffs, len1, vec2[i].coeffs, len2, mod);
            else
                _nmod_poly_mul(t, vec2[i].coeffs, len2, vec1[i].coeffs, len1, mod);

            _nmod_poly_add(s, s, slen, t, tlen, mod);
            slen = FLINT_MAX(slen, tlen);
        }
    }

    if (initial == NULL)
    {
        if (subtract)
            _nmod_vec_neg(s, s, slen, mod);
    }
    else
    {
        len2 = initial->length;

        if (subtract)
            _nmod_poly_sub(s, initial->coeffs, len2, s, slen, mod);
        else
            _nmod_poly_add(s, initial->coeffs, len2, s, slen, mod);

        slen = FLINT_MAX(slen, len2);
    }

    while (slen > 0 && s[slen - 1] == 0)
        slen--;

    _fq_nmod_reduce(s, slen, FQ_CTX(ctx));
    slen = FLINT_MIN(slen, plen - 1);

    while (slen > 0 && s[slen - 1] == 0)
        slen--;

    nmod_poly_fit_length(res, slen);
    _nmod_vec_set(res->coeffs, s, slen);
    _nmod_poly_set_length(res, slen);

    GR_TMP_FREE(t, (4 * plen) * sizeof(mp_limb_t));

    return GR_SUCCESS;
}

/* todo: basecase multiplication without reductions */
int
__gr_fq_nmod_vec_dot_rev(fq_nmod_struct * res, const fq_nmod_struct * initial, int subtract, const fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    mp_ptr s, t;
    slong slen, tlen, len1, len2;
    slong plen;
    nmod_t mod;

    if (len <= 0)
    {
        if (initial == NULL)
            fq_nmod_zero(res, FQ_CTX(ctx));
        else
            fq_nmod_set(res, initial, FQ_CTX(ctx));
        return GR_SUCCESS;
    }

    plen = FQ_CTX(ctx)->modulus->length;

    t = GR_TMP_ALLOC((4 * plen) * sizeof(mp_limb_t));
    s = t + 2 * plen;

    mod = FQ_CTX(ctx)->mod;

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
            _nmod_poly_mul(s, vec1[0].coeffs, len1, vec2[len - 1].coeffs, len2, mod);
        else
            _nmod_poly_mul(s, vec2[len - 1].coeffs, len2, vec1[0].coeffs, len1, mod);
    }

    for (i = 1; i < len; i++)
    {
        len1 = vec1[i].length;
        len2 = vec2[len - 1 - i].length;

        if (len1 != 0 && len2 != 0)
        {
            tlen = len1 + len2 - 1;
            if (len1 >= len2)
                _nmod_poly_mul(t, vec1[i].coeffs, len1, vec2[len - 1 - i].coeffs, len2, mod);
            else
                _nmod_poly_mul(t, vec2[len - 1 - i].coeffs, len2, vec1[i].coeffs, len1, mod);

            _nmod_poly_add(s, s, slen, t, tlen, mod);
            slen = FLINT_MAX(slen, tlen);
        }
    }

    if (initial == NULL)
    {
        if (subtract)
            _nmod_vec_neg(s, s, slen, mod);
    }
    else
    {
        len2 = initial->length;

        if (subtract)
            _nmod_poly_sub(s, initial->coeffs, len2, s, slen, mod);
        else
            _nmod_poly_add(s, initial->coeffs, len2, s, slen, mod);

        slen = FLINT_MAX(slen, len2);
    }

    while (slen > 0 && s[slen - 1] == 0)
        slen--;

    _fq_nmod_reduce(s, slen, FQ_CTX(ctx));
    slen = FLINT_MIN(slen, plen - 1);

    while (slen > 0 && s[slen - 1] == 0)
        slen--;

    nmod_poly_fit_length(res, slen);
    _nmod_vec_set(res->coeffs, s, slen);
    _nmod_poly_set_length(res, slen);

    GR_TMP_FREE(t, (4 * plen) * sizeof(mp_limb_t));

    return GR_SUCCESS;
}

/* todo: _fq_nmod_poly_mullow should do the right thing */
int
_gr_fq_nmod_poly_mullow(fq_nmod_struct * res,
    const fq_nmod_struct * poly1, slong len1,
    const fq_nmod_struct * poly2, slong len2, slong n, gr_ctx_t ctx)
{
    if (len1 + len2 - 1 == n)
    {
        if (poly1 == poly2 && len1 == len2)
            _fq_nmod_poly_sqr(res, poly1, len1, FQ_CTX(ctx));
        else if (len1 >= len2)
            _fq_nmod_poly_mul(res, poly1, len1, poly2, len2, FQ_CTX(ctx));
        else
            _fq_nmod_poly_mul(res, poly2, len2, poly1, len1, FQ_CTX(ctx));
    }
    else
    {
        if (len1 >= len2)
            _fq_nmod_poly_mullow(res, poly1, len1, poly2, len2, n, FQ_CTX(ctx));
        else
            _fq_nmod_poly_mullow(res, poly2, len2, poly1, len1, n, FQ_CTX(ctx));
    }

    return GR_SUCCESS;
}

/* todo: also need the _other version ... ? */
/* todo: implement generically */

int
_gr_fq_nmod_roots_gr_poly(gr_vec_t roots, gr_vec_t mult, const fq_nmod_poly_t poly, int flags, gr_ctx_t ctx)
{
    if (poly->length == 0)
        return GR_DOMAIN;

    {
        gr_ctx_t ZZ;
        fq_nmod_poly_factor_t fac;
        slong i, num;

        gr_ctx_init_fmpz(ZZ);
        fq_nmod_poly_factor_init(fac, FQ_CTX(ctx));
        fq_nmod_poly_roots(fac, poly, 1, FQ_CTX(ctx));

        num = fac->num;

        gr_vec_set_length(roots, num, ctx);
        gr_vec_set_length(mult, num, ZZ);

        for (i = 0; i < num; i++)
        {
            fq_nmod_neg(gr_vec_entry_ptr(roots, i, ctx), fac->poly[i].coeffs, FQ_CTX(ctx));

            /* work around flint bug: factors can be non-monic */
            if (!fq_nmod_is_one(fac->poly[i].coeffs + 1, FQ_CTX(ctx)))
                fq_nmod_div(gr_vec_entry_ptr(roots, i, ctx), gr_vec_entry_ptr(roots, i, ctx), fac->poly[i].coeffs + 1, FQ_CTX(ctx));

            fmpz_set_ui(((fmpz *) mult->entries) + i, fac->exp[i]);
        }

        fq_nmod_poly_factor_clear(fac, FQ_CTX(ctx));
        gr_ctx_clear(ZZ);
    }

    return GR_SUCCESS;
}

int
_gr_fq_nmod_mat_mul(fq_nmod_mat_t res, const fq_nmod_mat_t x, const fq_nmod_mat_t y, gr_ctx_t ctx)
{
    fq_nmod_mat_mul(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int _fq_nmod_methods_initialized = 0;

gr_static_method_table _fq_nmod_methods;

gr_method_tab_input _fq_nmod_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_fq_nmod_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_fq_nmod_ctx_write},
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
    {GR_METHOD_INIT,            (gr_funcptr) _gr_fq_nmod_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_fq_nmod_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_fq_nmod_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_fq_nmod_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_fq_nmod_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_fq_nmod_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_fq_nmod_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_fq_nmod_one},
    {GR_METHOD_GEN,             (gr_funcptr) _gr_fq_nmod_gen},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_fq_nmod_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_fq_nmod_is_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_fq_nmod_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_fq_nmod_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_fq_nmod_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_fq_nmod_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_fq_nmod_set_fmpz},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_fq_nmod_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_fq_nmod_add},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_fq_nmod_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_fq_nmod_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_fq_nmod_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_fq_nmod_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_fq_nmod_mul_fmpz},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fq_nmod_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fq_nmod_inv},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_fq_nmod_div},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_fq_nmod_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_fq_nmod_sqrt},

    {GR_METHOD_CTX_FQ_PRIME,            (gr_funcptr) _gr_ctx_fq_nmod_prime},
    {GR_METHOD_CTX_FQ_DEGREE,           (gr_funcptr) _gr_ctx_fq_nmod_degree},
    {GR_METHOD_CTX_FQ_ORDER,            (gr_funcptr) _gr_ctx_fq_nmod_order},
    {GR_METHOD_FQ_FROBENIUS,            (gr_funcptr) _gr_fq_nmod_frobenius},
    {GR_METHOD_FQ_MULTIPLICATIVE_ORDER, (gr_funcptr) _gr_fq_nmod_multiplicative_order},
    {GR_METHOD_FQ_NORM,                 (gr_funcptr) _gr_fq_nmod_norm},
    {GR_METHOD_FQ_TRACE,                (gr_funcptr) _gr_fq_nmod_trace},
    {GR_METHOD_FQ_IS_PRIMITIVE,         (gr_funcptr) _gr_fq_nmod_is_primitive},
    {GR_METHOD_FQ_PTH_ROOT,             (gr_funcptr) _gr_fq_nmod_pth_root},

    {GR_METHOD_VEC_DOT,         (gr_funcptr) __gr_fq_nmod_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) __gr_fq_nmod_vec_dot_rev},

    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_fq_nmod_poly_mullow},
    {GR_METHOD_POLY_ROOTS,      (gr_funcptr) _gr_fq_nmod_roots_gr_poly},

    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_fq_nmod_mat_mul},
    {0,                         (gr_funcptr) NULL},
};

void
_gr_ctx_init_fq_nmod_from_ref(gr_ctx_t ctx, const void * fq_nmod_ctx)
{
    ctx->which_ring = GR_CTX_FQ_NMOD;
    ctx->sizeof_elem = sizeof(fq_nmod_struct);
    GR_CTX_DATA_AS_PTR(ctx) = (fq_nmod_ctx_struct *) fq_nmod_ctx;
    ctx->size_limit = WORD_MAX;
    ctx->methods = _fq_nmod_methods;

    if (!_fq_nmod_methods_initialized)
    {
        gr_method_tab_init(_fq_nmod_methods, _fq_nmod_methods_input);
        _fq_nmod_methods_initialized = 1;
    }
}

void
gr_ctx_init_fq_nmod(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var)
{
    fq_nmod_ctx_struct * fq_nmod_ctx;

    if (fmpz_bits(p) > FLINT_BITS)
    {
        flint_throw(FLINT_ERROR, "gr_ctx_init_fq_nmod: expected a word-size p\n");
    }

    fq_nmod_ctx = flint_malloc(sizeof(fq_nmod_ctx_struct));
    fq_nmod_ctx_init(fq_nmod_ctx, p, d, var == NULL ? "a" : var);
    _gr_ctx_init_fq_nmod_from_ref(ctx, fq_nmod_ctx);
}
