/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "fmpz.h"
#include "fq.h"
#include "fq_nmod.h"
#include "fq_nmod_mat.h"
#include "fq_nmod_poly.h"
#include "fq_nmod_poly_factor.h"
#include "fmpz_mod_poly.h"
#include "gr.h"
#include "gr/impl.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"
#include "gr_generic.h"

#define FQ_CTX(ring_ctx) ((fq_nmod_ctx_struct *)(GR_CTX_DATA_AS_PTR(ring_ctx)))

static const char * default_var = "a";

static void
_gr_fq_nmod_ctx_clear(gr_ctx_t ctx)
{
    fq_nmod_ctx_clear(FQ_CTX(ctx));
    flint_free(GR_CTX_DATA_AS_PTR(ctx));
}

static int
_gr_fq_nmod_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_stream_write(out, "GF(");
    status |= gr_stream_write_ui(out, fq_nmod_ctx_prime(FQ_CTX(ctx)));
    status |= gr_stream_write(out, "^");
    status |= gr_stream_write_si(out, fq_nmod_ctx_degree(FQ_CTX(ctx)));
    status |= gr_stream_write(out, ") (fq_nmod)");
    return status;
}

static int _gr_fq_nmod_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
{
    slong len;
    len = strlen(s);

    FQ_CTX(ctx)->var = flint_realloc(FQ_CTX(ctx)->var, len + 1);
    memcpy(FQ_CTX(ctx)->var, s, len + 1);
    return GR_SUCCESS;
}

static int _gr_fq_nmod_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)
{
    return _gr_fq_nmod_ctx_set_gen_name(ctx, s[0]);
}

static int
_gr_fq_nmod_ctx_gen_name(char ** name, slong i, gr_ctx_t ctx)
{
    if (i != 0)
        return GR_DOMAIN;

    char * var = FQ_CTX(ctx)->var;
    size_t len = strlen(var);
    * name = flint_malloc(len + 1);
    if (* name == NULL)
        return GR_UNABLE;
    strncpy(* name, var, len + 1);

    return GR_SUCCESS;
}

static void
_gr_fq_nmod_init(fq_nmod_t x, const gr_ctx_t ctx)
{
    fq_nmod_init(x, FQ_CTX(ctx));
}

static void
_gr_fq_nmod_clear(fq_nmod_t x, const gr_ctx_t ctx)
{
    fq_nmod_clear(x, FQ_CTX(ctx));
}

static void
_gr_fq_nmod_swap(fq_nmod_t x, fq_nmod_t y, const gr_ctx_t ctx)
{
    fq_nmod_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

static void
_gr_fq_nmod_set_shallow(fq_nmod_t res, const fq_nmod_t x, const gr_ctx_t ctx)
{
    *res = *x;
}

static int
_gr_fq_nmod_randtest(fq_nmod_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    fq_nmod_randtest(res, state, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_write(gr_stream_t out, const fq_nmod_t x, const gr_ctx_t ctx)
{
    return gr_stream_write_free(out, fq_nmod_get_str_pretty(x, FQ_CTX(ctx)));
}

static int
_gr_fq_nmod_zero(fq_nmod_t x, const gr_ctx_t ctx)
{
    fq_nmod_zero(x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_one(fq_nmod_t x, const gr_ctx_t ctx)
{
    fq_nmod_one(x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_set_si(fq_nmod_t res, slong v, const gr_ctx_t ctx)
{
    fq_nmod_set_si(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_set_ui(fq_nmod_t res, ulong v, const gr_ctx_t ctx)
{
    fq_nmod_set_ui(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_set_fmpz(fq_nmod_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    fq_nmod_set_fmpz(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

/*
    An element of F_q is an integer exactly when it lies in the prime
    field, in which case it is its own constant coefficient.
*/
static int
_gr_fq_nmod_get_fmpz(fmpz_t res, const fq_nmod_t x, const gr_ctx_t ctx)
{
    /* an fq_nmod_t is an nmod_poly_t, reduced, so this is exact */
    if (nmod_poly_length(x) > 1)
        return GR_DOMAIN;

    fmpz_set_ui(res, nmod_poly_get_coeff_ui(x, 0));
    return GR_SUCCESS;
}

static truth_t
_gr_fq_nmod_is_zero(const fq_nmod_t x, const gr_ctx_t ctx)
{
    return fq_nmod_is_zero(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_fq_nmod_is_one(const fq_nmod_t x, const gr_ctx_t ctx)
{
    return fq_nmod_is_one(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_fq_nmod_equal(const fq_nmod_t x, const fq_nmod_t y, const gr_ctx_t ctx)
{
    return fq_nmod_equal(x, y, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

static int
_gr_fq_nmod_set(fq_nmod_t res, const fq_nmod_t x, const gr_ctx_t ctx)
{
    fq_nmod_set(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_neg(fq_nmod_t res, const fq_nmod_t x, const gr_ctx_t ctx)
{
    fq_nmod_neg(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_add(fq_nmod_t res, const fq_nmod_t x, const fq_nmod_t y, const gr_ctx_t ctx)
{
    fq_nmod_add(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_sub(fq_nmod_t res, const fq_nmod_t x, const fq_nmod_t y, const gr_ctx_t ctx)
{
    fq_nmod_sub(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_mul(fq_nmod_t res, const fq_nmod_t x, const fq_nmod_t y, const gr_ctx_t ctx)
{
    fq_nmod_mul(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_mul_si(fq_nmod_t res, const fq_nmod_t x, slong y, const gr_ctx_t ctx)
{
    fq_nmod_mul_si(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

/* fq_nmod_mul_ui without the normalisation (used for example by
   _gr_poly_derivative) */
static int
_gr_fq_nmod_mul_ui(fq_nmod_t res, const fq_nmod_t x, ulong y, const gr_ctx_t ctx)
{
    nmod_t mod = FQ_CTX(ctx)->mod;
    slong len = x->length;

    if (y >= mod.n)
        NMOD_RED(y, y, mod);

    if (y == 0 || len == 0)
    {
        res->length = 0;
    }
    else if (y == 1)
    {
        nmod_poly_set(res, x);
    }
    else
    {
        if (res->alloc < len)
            nmod_poly_fit_length(res, len);

        /* no cancellation as the characteristic is prime */
        _nmod_vec_scalar_mul_nmod(res->coeffs, x->coeffs, len, y, mod);
        res->length = len;
    }

    return GR_SUCCESS;
}

static int
_gr_fq_nmod_mul_fmpz(fq_nmod_t res, const fq_nmod_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fq_nmod_mul_fmpz(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_vec_mul_scalar_2exp_si(fq_nmod_struct * res, const fq_nmod_struct * vec, slong len, slong c, gr_ctx_t ctx)
{
    nmod_t mod = FQ_CTX(ctx)->mod;
    ulong t;
    slong i;

    if (c == 0)
    {
        for (i = 0; i < len; i++)
            fq_nmod_set(res + i, vec + i, FQ_CTX(ctx));
        return GR_SUCCESS;
    }

    if (mod.n == 2)
    {
        if (c < 0)
            return GR_DOMAIN;

        for (i = 0; i < len; i++)
            fq_nmod_zero(res + i, FQ_CTX(ctx));
        return GR_SUCCESS;
    }

    /* 1/2 = (p + 1) / 2 */
    if (c > 0)
        t = nmod_pow_ui(2, c, mod);
    else
        t = nmod_pow_ui((mod.n + 1) / 2, -(ulong) c, mod);

    for (i = 0; i < len; i++)
        nmod_poly_scalar_mul_nmod(res + i, vec + i, t);

    return GR_SUCCESS;
}

static int
_gr_fq_nmod_mul_2exp_si(fq_nmod_t res, const fq_nmod_t x, slong c, gr_ctx_t ctx)
{
    return _gr_fq_nmod_vec_mul_scalar_2exp_si(res, x, 1, c, ctx);
}

static int
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

static int
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

static int
_gr_fq_nmod_sqr(fq_nmod_t res, const fq_nmod_t x, const gr_ctx_t ctx)
{
    fq_nmod_sqr(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_pow_ui(fq_nmod_t res, const fq_nmod_t x, ulong y, const gr_ctx_t ctx)
{
    fq_nmod_pow_ui(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_pow_fmpz(fq_nmod_t res, const fq_nmod_t x, const fmpz_t y, gr_ctx_t ctx)
{
    if (fmpz_sgn(y) < 0)
    {
        return gr_generic_pow_fmpz(res, x, y, ctx);
    }
    else
    {
        fq_nmod_pow(res, x, y, FQ_CTX(ctx));
        return GR_SUCCESS;
    }
}

static truth_t
_gr_fq_nmod_is_invertible(const fq_nmod_t x, const gr_ctx_t ctx)
{
    return (!fq_nmod_is_zero(x, FQ_CTX(ctx))) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_fq_nmod_is_square(const fq_nmod_t x, const gr_ctx_t ctx)
{
    return fq_nmod_is_square(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

static int
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

static int
_gr_ctx_fq_nmod_prime(fmpz_t p, gr_ctx_t ctx)
{
    fmpz_set_ui(p, fq_nmod_ctx_prime(FQ_CTX(ctx)));
    return GR_SUCCESS;
}

static int
_gr_ctx_fq_nmod_degree(slong * deg, gr_ctx_t ctx)
{
    *deg = fq_nmod_ctx_degree(FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_ctx_fq_nmod_order(fmpz_t q, gr_ctx_t ctx)
{
    fq_nmod_ctx_order(q, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_gen(gr_ptr res, gr_ctx_t ctx)
{
    fq_nmod_gen(res, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_frobenius(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx)
{
    fq_nmod_frobenius(res, x, e, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_multiplicative_order(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    int ret;
    ret = fq_nmod_multiplicative_order(res, x, FQ_CTX(ctx));

    if (ret == 1)
        return GR_SUCCESS;

    /* todo: better solution? */
    return GR_DOMAIN;
}

static int
_gr_fq_nmod_norm(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_nmod_norm(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_trace(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_nmod_trace(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static truth_t
_gr_fq_nmod_is_primitive(gr_srcptr x, gr_ctx_t ctx)
{
    return fq_nmod_is_primitive(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

static int
_gr_fq_nmod_pth_root(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_nmod_pth_root(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

/* todo: basecase multiplication without reductions */
static int
__gr_fq_nmod_vec_dot(fq_nmod_struct * res, const fq_nmod_struct * initial, int subtract, const fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    nn_ptr s, t;
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

    t = GR_TMP_ALLOC((4 * plen) * sizeof(ulong));
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

    GR_TMP_FREE(t, (4 * plen) * sizeof(ulong));

    return GR_SUCCESS;
}

/* todo: basecase multiplication without reductions */
static int
__gr_fq_nmod_vec_dot_rev(fq_nmod_struct * res, const fq_nmod_struct * initial, int subtract, const fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    nn_ptr s, t;
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

    t = GR_TMP_ALLOC((4 * plen) * sizeof(ulong));
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

    GR_TMP_FREE(t, (4 * plen) * sizeof(ulong));

    return GR_SUCCESS;
}

/* todo: _fq_nmod_poly_mullow should do the right thing */
/* gcd and xgcd of the fq_nmod_poly module (Euclid or half-gcd with its
   cutoffs) */
static int
_gr_fq_nmod_poly_gcd(fq_nmod_struct * G, slong * lenG, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, gr_ctx_t ctx)
{
    *lenG = _fq_nmod_poly_gcd(G, A, lenA, B, lenB, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_poly_xgcd(slong * lenG, fq_nmod_struct * G, fq_nmod_struct * S, fq_nmod_struct * T, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, gr_ctx_t ctx)
{
    *lenG = _fq_nmod_poly_xgcd(G, S, T, A, lenA, B, lenB, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
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

/* The middle product is computed by the classical algorithm in the range
   where the classical full or low product would be used anyway. */
static int
_gr_fq_nmod_poly_mulmid(fq_nmod_struct * res,
    const fq_nmod_struct * poly1, slong len1,
    const fq_nmod_struct * poly2, slong len2, slong nlo, slong nhi, gr_ctx_t ctx)
{
    if (nlo != 0 && (nhi < FQ_NMOD_MULLOW_CLASSICAL_CUTOFF ||
                     FLINT_MAX(len1, len2) < FQ_NMOD_MUL_CLASSICAL_CUTOFF))
        return _gr_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, ctx);
    else
        return _gr_poly_mulmid_generic(res, poly1, len1, poly2, len2, nlo, nhi, ctx);
}

static int
_gr_fq_nmod_poly_divrem(fq_nmod_struct * Q, fq_nmod_struct * R,
    const fq_nmod_struct * A, slong lenA,
    const fq_nmod_struct * B, slong lenB, gr_ctx_t ctx)
{
    /* The Newton division only involves products of length about lenQ and
       lenB, which are fast (Kronecker substitution), and in this ring it
       beats the basecase division already for tiny lengths. */
    return _gr_poly_divrem_newton(Q, R, A, lenA, B, lenB, ctx);
}

/* Horner's rule */
static int
_gr_fq_nmod_poly_evaluate(fq_nmod_t res, const fq_nmod_struct * f, slong len,
    const fq_nmod_t x, gr_ctx_t ctx)
{
    const fq_nmod_ctx_struct * fctx = FQ_CTX(ctx);

    if (len == 0)
    {
        fq_nmod_zero(res, fctx);
    }
    else if (len == 1 || fq_nmod_is_zero(x, fctx))
    {
        fq_nmod_set(res, f, fctx);
    }
    else
    {
        fq_nmod_t t, u;
        fq_nmod_struct * s;
        slong i;

        fq_nmod_init(t, fctx);

        if (res == x)
        {
            fq_nmod_init(u, fctx);
            s = u;
        }
        else
        {
            s = res;
        }

        fq_nmod_set(s, f + len - 1, fctx);

        for (i = len - 2; i >= 0; i--)
        {
            fq_nmod_mul(t, s, x, fctx);
            fq_nmod_add(s, f + i, t, fctx);
        }

        if (res == x)
        {
            fq_nmod_swap(res, u, fctx);
            fq_nmod_clear(u, fctx);
        }

        fq_nmod_clear(t, fctx);
    }

    return GR_SUCCESS;
}

/* todo: also need the _other version ... ? */
/* todo: implement generically */

static int
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

static int
_gr_fq_nmod_mat_mul(fq_nmod_mat_t res, const fq_nmod_mat_t x, const fq_nmod_mat_t y, gr_ctx_t ctx)
{
    fq_nmod_mat_mul(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fq_nmod_mat_nonsingular_solve_tril(fq_nmod_mat_t X, const fq_nmod_mat_t L, const fq_nmod_mat_t B, int unit, gr_ctx_t ctx)
{
    if (B->r < 64 || B->c < 64)
        return gr_mat_nonsingular_solve_tril_classical((gr_mat_struct *) X, (const gr_mat_struct *) L, (const gr_mat_struct *) B, unit, ctx);
    else
        return gr_mat_nonsingular_solve_tril_recursive((gr_mat_struct *) X, (const gr_mat_struct *) L, (const gr_mat_struct *) B, unit, ctx);
}

static int
_gr_fq_nmod_mat_nonsingular_solve_triu(fq_nmod_mat_t X, const fq_nmod_mat_t U, const fq_nmod_mat_t B, int unit, gr_ctx_t ctx)
{
    if (B->r < 64 || B->c < 64)
        return gr_mat_nonsingular_solve_triu_classical((gr_mat_struct *) X, (const gr_mat_struct *) U, (const gr_mat_struct *) B, unit, ctx);
    else
        return gr_mat_nonsingular_solve_triu_recursive((gr_mat_struct *) X, (const gr_mat_struct *) U, (const gr_mat_struct *) B, unit, ctx);
}

static int
_gr_fq_nmod_mat_charpoly(fq_nmod_struct * res, const fq_nmod_mat_t mat, gr_ctx_t ctx)
{
    slong n = mat->r;

    if (n <= 12)
        return _gr_mat_charpoly_berkowitz(res, (const gr_mat_struct *) mat, ctx);
    else
        return _gr_mat_charpoly_danilevsky(res, (const gr_mat_struct *) mat, ctx);
}

static int
_gr_fq_nmod_mat_reduce_row(slong * column, fq_nmod_mat_t mat, slong * P, slong * L, slong n, gr_ctx_t ctx)
{
    *column = fq_nmod_mat_reduce_row(mat, P, L, n, FQ_CTX(ctx));
    return GR_SUCCESS;
}

/* Vector methods avoiding per-element dispatch through the method table */

static int
_gr_fq_nmod_vec_set(fq_nmod_struct * res, const fq_nmod_struct * vec, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        fq_nmod_set(res + i, vec + i, FQ_CTX(ctx));

    return GR_SUCCESS;
}

static int
_gr_fq_nmod_vec_zero(fq_nmod_struct * res, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        fq_nmod_zero(res + i, FQ_CTX(ctx));

    return GR_SUCCESS;
}

static int
_gr_fq_nmod_vec_add(fq_nmod_struct * res, const fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        fq_nmod_add(res + i, vec1 + i, vec2 + i, FQ_CTX(ctx));

    return GR_SUCCESS;
}

static int
_gr_fq_nmod_vec_sub(fq_nmod_struct * res, const fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        fq_nmod_sub(res + i, vec1 + i, vec2 + i, FQ_CTX(ctx));

    return GR_SUCCESS;
}

static int
_gr_fq_nmod_vec_add_scalar(fq_nmod_struct * res, const fq_nmod_struct * vec, slong len, const fq_nmod_struct * c, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        fq_nmod_add(res + i, vec + i, c, FQ_CTX(ctx));

    return GR_SUCCESS;
}

static truth_t
_gr_fq_nmod_vec_is_zero(const fq_nmod_struct * vec, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        if (!fq_nmod_is_zero(vec + i, FQ_CTX(ctx)))
            return T_FALSE;

    return T_TRUE;
}

static truth_t
_gr_fq_nmod_vec_equal(const fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        if (!fq_nmod_equal(vec1 + i, vec2 + i, FQ_CTX(ctx)))
            return T_FALSE;

    return T_TRUE;
}

int _fq_nmod_methods_initialized = 0;

gr_static_method_table _fq_nmod_methods;

gr_method_tab_input _fq_nmod_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_fq_nmod_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_fq_nmod_ctx_write},
    {GR_METHOD_CTX_SET_GEN_NAME,    (gr_funcptr) _gr_fq_nmod_ctx_set_gen_name},
    {GR_METHOD_CTX_SET_GEN_NAMES,   (gr_funcptr) _gr_fq_nmod_ctx_set_gen_names},
    {GR_METHOD_CTX_NGENS,       (gr_funcptr) gr_generic_ctx_ngens_1},
    {GR_METHOD_CTX_GEN_NAME,    (gr_funcptr) _gr_fq_nmod_ctx_gen_name},
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
    {GR_METHOD_GENS,            (gr_funcptr) gr_generic_gens_single},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_fq_nmod_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_fq_nmod_is_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_fq_nmod_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_fq_nmod_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_fq_nmod_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_fq_nmod_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_fq_nmod_set_fmpz},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_fq_nmod_get_fmpz},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_fq_nmod_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_fq_nmod_add},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_fq_nmod_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_fq_nmod_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_fq_nmod_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_fq_nmod_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_fq_nmod_mul_fmpz},
    {GR_METHOD_MUL_2EXP_SI,     (gr_funcptr) _gr_fq_nmod_mul_2exp_si},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fq_nmod_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fq_nmod_inv},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_fq_nmod_div},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_fq_nmod_sqr},
    {GR_METHOD_POW_UI,           (gr_funcptr) _gr_fq_nmod_pow_ui},
    {GR_METHOD_POW_FMPZ,         (gr_funcptr) _gr_fq_nmod_pow_fmpz},

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
    {GR_METHOD_VEC_MUL_SCALAR_2EXP_SI,  (gr_funcptr) _gr_fq_nmod_vec_mul_scalar_2exp_si},
    {GR_METHOD_VEC_SET,         (gr_funcptr) _gr_fq_nmod_vec_set},
    {GR_METHOD_VEC_ZERO,        (gr_funcptr) _gr_fq_nmod_vec_zero},
    {GR_METHOD_VEC_ADD,         (gr_funcptr) _gr_fq_nmod_vec_add},
    {GR_METHOD_VEC_SUB,         (gr_funcptr) _gr_fq_nmod_vec_sub},
    {GR_METHOD_VEC_ADD_SCALAR,  (gr_funcptr) _gr_fq_nmod_vec_add_scalar},
    {GR_METHOD_VEC_IS_ZERO,     (gr_funcptr) _gr_fq_nmod_vec_is_zero},
    {GR_METHOD_VEC_EQUAL,       (gr_funcptr) _gr_fq_nmod_vec_equal},

    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_fq_nmod_poly_mullow},
    {GR_METHOD_POLY_MULMID,     (gr_funcptr) _gr_fq_nmod_poly_mulmid},
    {GR_METHOD_POLY_DIVREM,     (gr_funcptr) _gr_fq_nmod_poly_divrem},
    {GR_METHOD_POLY_EVALUATE,   (gr_funcptr) _gr_fq_nmod_poly_evaluate},
    {GR_METHOD_POLY_GCD,        (gr_funcptr) _gr_fq_nmod_poly_gcd},
    {GR_METHOD_POLY_XGCD,       (gr_funcptr) _gr_fq_nmod_poly_xgcd},

    {GR_METHOD_POLY_FACTOR,     (gr_funcptr) _gr_poly_factor_finite_field_method},
    {GR_METHOD_POLY_ROOTS,      (gr_funcptr) _gr_fq_nmod_roots_gr_poly},

    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_fq_nmod_mat_mul},
    {GR_METHOD_MAT_NONSINGULAR_SOLVE_TRIL,      (gr_funcptr) _gr_fq_nmod_mat_nonsingular_solve_tril},
    {GR_METHOD_MAT_NONSINGULAR_SOLVE_TRIU,      (gr_funcptr) _gr_fq_nmod_mat_nonsingular_solve_triu},
    {GR_METHOD_MAT_CHARPOLY,    (gr_funcptr) _gr_fq_nmod_mat_charpoly},
    {GR_METHOD_MAT_REDUCE_ROW,  (gr_funcptr) _gr_fq_nmod_mat_reduce_row},
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
gr_ctx_init_fq_nmod(gr_ctx_t ctx, ulong p, slong d, const char * var)
{
    fq_nmod_ctx_struct * fq_nmod_ctx;

    fq_nmod_ctx = flint_malloc(sizeof(fq_nmod_ctx_struct));
    fq_nmod_ctx_init_ui(fq_nmod_ctx, p, d, var == NULL ? default_var : var);
    _gr_ctx_init_fq_nmod_from_ref(ctx, fq_nmod_ctx);
}

int gr_ctx_init_fq_nmod_modulus_nmod_poly(gr_ctx_t ctx, const nmod_poly_t modulus, const char * var)
{
    fq_nmod_ctx_struct * fq_nmod_ctx;
    fq_nmod_ctx = flint_malloc(sizeof(fq_nmod_ctx_struct));
    fq_nmod_ctx_init_modulus(fq_nmod_ctx, modulus, var == NULL ? default_var : var);
    _gr_ctx_init_fq_nmod_from_ref(ctx, fq_nmod_ctx);
    return GR_SUCCESS;
}

int
gr_ctx_init_fq_nmod_modulus_fmpz_mod_poly(gr_ctx_t ctx, const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var)
{
    nmod_poly_t nmodulus;
    int status;

    if (!fmpz_abs_fits_ui(mod_ctx->n))
        return GR_UNABLE;

    nmod_poly_init(nmodulus, fmpz_get_ui(mod_ctx->n));
    fmpz_mod_poly_get_nmod_poly(nmodulus, modulus);
    status = gr_ctx_init_fq_nmod_modulus_nmod_poly(ctx, nmodulus, var);
    nmod_poly_clear(nmodulus);
    return status;
}
