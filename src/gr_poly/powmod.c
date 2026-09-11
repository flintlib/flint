/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart
    Copyright (C) 2011, 2025 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "long_extras.h"
#include "fmpz.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"

/* todo: this might want to be a public method */
/*
    out = a * b mod f, using the actual lengths of a and b (which must
    be reduced, i.e. less than lenf) and skipping the reduction when the
    product already has degree less than deg(f). Writes the normalised
    length to *lenout. Requires space for lenf - 1 coefficients in out,
    which must not alias a or b, and scratch space T for 3 (lenf - 1) - 1
    coefficients (the product and the quotient).
*/
static int
_gr_poly_mulmod_len(gr_ptr out, slong * lenout,
    gr_srcptr a, slong lena, gr_srcptr b, slong lenb,
    const gr_poly_preinv_t P, gr_ptr T, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong lenf = P->lenf, lenT;

    if (lena == 0 || lenb == 0)
    {
        *lenout = 0;
        return GR_SUCCESS;
    }

    lenT = lena + lenb - 1;

    if (lenT < lenf)
    {
        if (lena >= lenb)
            status |= _gr_poly_mul(out, a, lena, b, lenb, ctx);
        else
            status |= _gr_poly_mul(out, b, lenb, a, lena, ctx);
        *lenout = lenT;
    }
    else
    {
        gr_ptr Q = GR_ENTRY(T, lenT, ctx->sizeof_elem);

        if (lena >= lenb)
            status |= _gr_poly_mul(T, a, lena, b, lenb, ctx);
        else
            status |= _gr_poly_mul(T, b, lenb, a, lena, ctx);

        status |= _gr_poly_preinv_divrem(Q, out, T, lenT, P, ctx);
        *lenout = lenf - 1;
    }

    GR_IGNORE(_gr_vec_normalise(lenout, out, *lenout, ctx));
    return status;
}

/* Pads res (holding len coefficients) with zeros to length lenf - 1,
   the output convention of the underscore powmod functions. */
static int
_gr_poly_powmod_pad(gr_ptr res, slong len, slong lenf, gr_ctx_t ctx)
{
    return _gr_vec_zero(GR_ENTRY(res, len, ctx->sizeof_elem), lenf - 1 - len, ctx);
}

/*
    Binary exponentiation with the actual lengths of the intermediate
    powers. For lenfinv = 0, no precomputed inverse is used.
    Requires len < lenf, lenf >= 2 and e >= 0. res may alias poly.
*/
int
_gr_poly_preinv_powmod_fmpz_binexp(gr_ptr res, gr_srcptr poly, slong len,
    const fmpz_t e, const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    gr_ptr W, a, b, t, T;
    slong lena, lenb, i;
    slong lenf = P->lenf;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    GR_IGNORE(_gr_vec_normalise(&len, poly, len, ctx));

    if (fmpz_is_zero(e))
    {
        status |= gr_one(res, ctx);
        return status | _gr_poly_powmod_pad(res, 1, lenf, ctx);
    }

    if (len == 0)
        return _gr_poly_powmod_pad(res, 0, lenf, ctx);

    if (fmpz_is_one(e))
    {
        if (res != poly)
            status |= _gr_vec_set(res, poly, len, ctx);
        return status | _gr_poly_powmod_pad(res, len, lenf, ctx);
    }

    GR_TMP_INIT_VEC(W, 5 * (lenf - 1) - 1, ctx);
    a = W;
    b = GR_ENTRY(W, lenf - 1, sz);
    T = GR_ENTRY(W, 2 * (lenf - 1), sz);

    /* a = poly^2 (top bit of e), then square-and-multiply */
    status |= _gr_poly_mulmod_len(a, &lena, poly, len, poly, len, P, T, ctx);

    for (i = fmpz_sizeinbase(e, 2) - 2; i >= 0 && status == GR_SUCCESS; i--)
    {
        if (fmpz_tstbit(e, i))
        {
            status |= _gr_poly_mulmod_len(b, &lenb, a, lena, poly, len, P, T, ctx);
            t = a; a = b; b = t; lena = lenb;
        }

        if (i > 0)
        {
            status |= _gr_poly_mulmod_len(b, &lenb, a, lena, a, lena, P, T, ctx);
            t = a; a = b; b = t; lena = lenb;
        }
    }

    status |= _gr_vec_set(res, a, lena, ctx);
    status |= _gr_poly_powmod_pad(res, lena, lenf, ctx);

    GR_TMP_CLEAR_VEC(W, 5 * (lenf - 1) - 1, ctx);
    return status;
}

int
_gr_poly_powmod_fmpz_binexp_preinv(gr_ptr res, gr_srcptr poly, slong len,
    const fmpz_t e, gr_srcptr f, slong lenf,
    gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)
{
    gr_poly_preinv_t P;
    _gr_poly_preinv_init_newton_shallow(P, f, lenf, finv, lenfinv, ctx);
    return _gr_poly_preinv_powmod_fmpz_binexp(res, poly, len, e, P, ctx);
}

int
_gr_poly_powmod_fmpz_binexp(gr_ptr res, gr_srcptr poly, slong len,
    const fmpz_t e, gr_srcptr f, slong lenf, gr_ctx_t ctx)
{
    return _gr_poly_powmod_fmpz_binexp_preinv(res, poly, len, e, f, lenf, NULL, 0, ctx);
}

int
_gr_poly_powmod_ui_binexp_preinv(gr_ptr res, gr_srcptr poly, slong len,
    ulong e, gr_srcptr f, slong lenf,
    gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)
{
    fmpz_t e2;
    int status;
    fmpz_init_set_ui(e2, e);
    status = _gr_poly_powmod_fmpz_binexp_preinv(res, poly, len, e2, f, lenf, finv, lenfinv, ctx);
    fmpz_clear(e2);
    return status;
}

int
_gr_poly_powmod_ui_binexp(gr_ptr res, gr_srcptr poly, slong len,
    ulong e, gr_srcptr f, slong lenf, gr_ctx_t ctx)
{
    return _gr_poly_powmod_ui_binexp_preinv(res, poly, len, e, f, lenf, NULL, 0, ctx);
}

/*
    Sliding window exponentiation with window size k (k = 0 selects the
    window size automatically) and the actual lengths of the intermediate
    powers. Requires len < lenf, lenf >= 2 and e >= 0. res may alias poly.
*/
int
_gr_poly_preinv_powmod_fmpz_sliding(gr_ptr res, gr_srcptr poly, slong len,
    const fmpz_t e, ulong k, const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    gr_ptr W, a, b, t, sq, precomp, T;
    slong * lenprecomp;
    slong lena, lenb, lensq, i, j, l, np;
    slong lenf = P->lenf;
    slong bits, sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    GR_IGNORE(_gr_vec_normalise(&len, poly, len, ctx));

    if (fmpz_is_zero(e))
    {
        status |= gr_one(res, ctx);
        return status | _gr_poly_powmod_pad(res, 1, lenf, ctx);
    }

    if (len == 0)
        return _gr_poly_powmod_pad(res, 0, lenf, ctx);

    if (fmpz_is_one(e))
    {
        if (res != poly)
            status |= _gr_vec_set(res, poly, len, ctx);
        return status | _gr_poly_powmod_pad(res, len, lenf, ctx);
    }

    bits = fmpz_bits(e);

    if (k == 0)
    {
        /* Window size k costs 2^(k-1) precomputed products and saves
           about bits/2 - bits/(k+1) multiplications compared to binary
           exponentiation; crossovers measured including the overhead. */
        if (bits <= 24)
            k = 1;
        else if (bits <= 64)
            k = 2;
        else if (bits <= 160)
            k = 3;
        else if (bits <= 480)
            k = 4;
        else if (bits <= 1300)
            k = 5;
        else if (bits <= 7000)
            k = 6;
        else
            k = 7;
    }

    if (k == 1)
        return _gr_poly_preinv_powmod_fmpz_binexp(res, poly, len, e, P, ctx);

    /* precomp[j] = poly^(2j+1) for 0 <= j < np */
    np = WORD(1) << (k - 1);

    GR_TMP_INIT_VEC(W, (np + 6) * (lenf - 1) - 1, ctx);
    lenprecomp = flint_malloc(np * sizeof(slong));
    precomp = W;
    sq = GR_ENTRY(W, np * (lenf - 1), sz);
    a = GR_ENTRY(sq, lenf - 1, sz);
    b = GR_ENTRY(a, lenf - 1, sz);
    T = GR_ENTRY(b, lenf - 1, sz);

    status |= _gr_vec_set(precomp, poly, len, ctx);
    lenprecomp[0] = len;
    status |= _gr_poly_mulmod_len(sq, &lensq, poly, len, poly, len, P, T, ctx);
    for (j = 1; j < np && status == GR_SUCCESS; j++)
        status |= _gr_poly_mulmod_len(GR_ENTRY(precomp, j * (lenf - 1), sz), lenprecomp + j,
            GR_ENTRY(precomp, (j - 1) * (lenf - 1), sz), lenprecomp[j - 1], sq, lensq, P, T, ctx);

    /* a = 1 */
    status |= gr_one(a, ctx);
    lena = 1;

    i = bits - 1;
    while (i >= 0 && status == GR_SUCCESS)
    {
        if (!fmpz_tstbit(e, i))
        {
            status |= _gr_poly_mulmod_len(b, &lenb, a, lena, a, lena, P, T, ctx);
            t = a; a = b; b = t; lena = lenb;
            i--;
        }
        else
        {
            /* find the largest window e[i..l] (with l >= i - k + 1) ending in a 1 */
            l = FLINT_MAX(i - (slong) k + 1, 0);
            while (!fmpz_tstbit(e, l))
                l++;

            for (j = 0; j < i - l + 1 && status == GR_SUCCESS; j++)
            {
                status |= _gr_poly_mulmod_len(b, &lenb, a, lena, a, lena, P, T, ctx);
                t = a; a = b; b = t; lena = lenb;
            }

            /* multiply by poly^(e[i..l]) which is odd */
            {
                ulong w = 0;
                for (j = i; j >= l; j--)
                    w = 2 * w + fmpz_tstbit(e, j);
                j = (w - 1) / 2;
                status |= _gr_poly_mulmod_len(b, &lenb, a, lena,
                    GR_ENTRY(precomp, j * (lenf - 1), sz), lenprecomp[j], P, T, ctx);
                t = a; a = b; b = t; lena = lenb;
            }

            i = l - 1;
        }
    }

    status |= _gr_vec_set(res, a, lena, ctx);
    status |= _gr_poly_powmod_pad(res, lena, lenf, ctx);

    flint_free(lenprecomp);
    GR_TMP_CLEAR_VEC(W, (np + 6) * (lenf - 1) - 1, ctx);
    return status;
}

int
_gr_poly_powmod_fmpz_sliding_preinv(gr_ptr res, gr_srcptr poly, slong len,
    const fmpz_t e, ulong k, gr_srcptr f, slong lenf,
    gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)
{
    gr_poly_preinv_t P;
    _gr_poly_preinv_init_newton_shallow(P, f, lenf, finv, lenfinv, ctx);
    return _gr_poly_preinv_powmod_fmpz_sliding(res, poly, len, e, k, P, ctx);
}

/*
    Common wrapper: reduces poly if necessary, handles trivial cases and
    calls the underscore function (which supports aliasing).
    algorithm: 0 = binexp, 1 = sliding (with window size k).
*/
static int
_gr_poly_preinv_powmod_wrapper(gr_poly_t res, const gr_poly_t poly, const fmpz_t e,
    const gr_poly_preinv_t P, int algorithm, ulong k, gr_ctx_t ctx)
{
    slong len = poly->length;
    slong lenf = P->lenf;
    int status = GR_SUCCESS;

    if (lenf == 0)
        return GR_DOMAIN;

    /* Not implemented. */
    if (fmpz_sgn(e) < 0)
        return GR_UNABLE;

    if (len >= lenf)
    {
        gr_poly_t r;
        gr_poly_init(r, ctx);
        status = gr_poly_preinv_rem(r, poly, P, ctx);
        if (status == GR_SUCCESS)
            status |= _gr_poly_preinv_powmod_wrapper(res, r, e, P, algorithm, k, ctx);
        gr_poly_clear(r, ctx);
        return status;
    }

    if (fmpz_is_zero(e))
    {
        if (lenf == 1)
            return gr_poly_zero(res, ctx);
        else
            return gr_poly_one(res, ctx);
    }

    if (lenf == 1 || len == 0)
        return gr_poly_zero(res, ctx);

    if (fmpz_is_one(e))
        return gr_poly_set(res, poly, ctx);

    gr_poly_fit_length(res, lenf - 1, ctx);

    if (algorithm == 0)
        status |= _gr_poly_preinv_powmod_fmpz_binexp(res->coeffs, poly->coeffs, len, e, P, ctx);
    else
        status |= _gr_poly_preinv_powmod_fmpz_sliding(res->coeffs, poly->coeffs, len, e, k, P, ctx);

    _gr_poly_set_length_normalise(res, lenf - 1, ctx);
    return status;
}

/* Wrapper for the functions taking f and (optionally) its Newton inverse
   as polynomials. */
static int
_gr_poly_powmod_wrapper(gr_poly_t res, const gr_poly_t poly, const fmpz_t e,
    const gr_poly_t f, const gr_poly_t finv, int algorithm, ulong k, gr_ctx_t ctx)
{
    gr_poly_preinv_t P;
    int status;

    if (f->length == 0)
        return GR_DOMAIN;

    if (res == f || res == finv)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status = _gr_poly_powmod_wrapper(t, poly, e, f, finv, algorithm, k, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    _gr_poly_preinv_init_newton_shallow(P, f->coeffs, f->length,
        finv == NULL ? NULL : finv->coeffs, finv == NULL ? 0 : finv->length, ctx);

    return _gr_poly_preinv_powmod_wrapper(res, poly, e, P, algorithm, k, ctx);
}

int
gr_poly_preinv_powmod_fmpz_binexp(gr_poly_t res, const gr_poly_t poly, const fmpz_t e,
    const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    return _gr_poly_preinv_powmod_wrapper(res, poly, e, P, 0, 0, ctx);
}

int
gr_poly_preinv_powmod_fmpz_sliding(gr_poly_t res, const gr_poly_t poly, const fmpz_t e,
    ulong k, const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    return _gr_poly_preinv_powmod_wrapper(res, poly, e, P, 1, k, ctx);
}

int
gr_poly_powmod_fmpz_binexp(gr_poly_t res, const gr_poly_t poly, const fmpz_t e,
    const gr_poly_t f, gr_ctx_t ctx)
{
    return _gr_poly_powmod_wrapper(res, poly, e, f, NULL, 0, 0, ctx);
}

int
gr_poly_powmod_fmpz_binexp_preinv(gr_poly_t res, const gr_poly_t poly, const fmpz_t e,
    const gr_poly_t f, const gr_poly_t finv, gr_ctx_t ctx)
{
    return _gr_poly_powmod_wrapper(res, poly, e, f, finv, 0, 0, ctx);
}

int
gr_poly_powmod_ui_binexp(gr_poly_t res, const gr_poly_t poly, ulong e,
    const gr_poly_t f, gr_ctx_t ctx)
{
    fmpz_t e2;
    int status;
    fmpz_init_set_ui(e2, e);
    status = _gr_poly_powmod_wrapper(res, poly, e2, f, NULL, 0, 0, ctx);
    fmpz_clear(e2);
    return status;
}

int
gr_poly_powmod_ui_binexp_preinv(gr_poly_t res, const gr_poly_t poly, ulong e,
    const gr_poly_t f, const gr_poly_t finv, gr_ctx_t ctx)
{
    fmpz_t e2;
    int status;
    fmpz_init_set_ui(e2, e);
    status = _gr_poly_powmod_wrapper(res, poly, e2, f, finv, 0, 0, ctx);
    fmpz_clear(e2);
    return status;
}

int
gr_poly_powmod_fmpz_sliding_preinv(gr_poly_t res, const gr_poly_t poly, const fmpz_t e,
    ulong k, const gr_poly_t f, const gr_poly_t finv, gr_ctx_t ctx)
{
    return _gr_poly_powmod_wrapper(res, poly, e, f, finv, 1, k, ctx);
}

int
_gr_poly_preinv_powmod_x_fmpz(gr_ptr res, const fmpz_t e,
    const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    gr_ptr W, a, b, t, T, Q;
    slong lena, lenb, lenT;
    slong lenf = P->lenf;
    slong i, window, l, c;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (fmpz_is_zero(e))
    {
        status |= gr_one(res, ctx);
        return status | _gr_poly_powmod_pad(res, 1, lenf, ctx);
    }

    /* a, b: current power and scratch (lenf - 1 each); T: product or
       shifted power (at most 2 lenf - 3) followed by the quotient */
    GR_TMP_INIT_VEC(W, 5 * (lenf - 1) - 1, ctx);
    a = W;
    b = GR_ENTRY(W, lenf - 1, sz);
    T = GR_ENTRY(W, 2 * (lenf - 1), sz);

    status |= gr_one(a, ctx);
    lena = 1;

    /* Windows of l + 1 bits, l = bits(deg f) - 2, so that x^window
       has degree less than deg f (window < 2^(l+1) <= deg f). */
    l = z_sizeinbase(lenf - 1, 2) - 2;
    window = (WORD(1) << l);
    c = l;
    i = fmpz_sizeinbase(e, 2) - 2;

    if (i <= l)
    {
        window = (WORD(1) << i);
        c = i;
        l = i;
    }

#define SHIFT_STEP() \
    do { \
        if (lena + window < lenf) \
        { \
            status |= _gr_poly_shift_left(b, a, lena, window, ctx); \
            lenb = lena + window; \
        } \
        else \
        { \
            lenT = lena + window; \
            Q = GR_ENTRY(T, lenT, sz); \
            status |= _gr_poly_shift_left(T, a, lena, window, ctx); \
            status |= _gr_poly_preinv_divrem(Q, b, T, lenT, P, ctx); \
            lenb = lenf - 1; \
            GR_IGNORE(_gr_vec_normalise(&lenb, b, lenb, ctx)); \
        } \
        t = a; a = b; b = t; lena = lenb; \
    } while (0)

    if (c == 0)
    {
        SHIFT_STEP();
        c = l + 1;
        window = 0;
    }

    for (; i >= 0 && status == GR_SUCCESS; i--)
    {
        status |= _gr_poly_mulmod_len(b, &lenb, a, lena, a, lena, P, T, ctx);
        t = a; a = b; b = t; lena = lenb;

        c--;
        if (fmpz_tstbit(e, i))
        {
            if (window == 0 && i <= l - 1)
                c = i;
            if (c >= 0)
                window = window | (WORD(1) << c);
        }
        else if (window == 0)
        {
            c = l + 1;
        }

        if (c == 0)
        {
            SHIFT_STEP();
            c = l + 1;
            window = 0;
        }
    }

#undef SHIFT_STEP

    status |= _gr_vec_set(res, a, lena, ctx);
    status |= _gr_poly_powmod_pad(res, lena, lenf, ctx);

    GR_TMP_CLEAR_VEC(W, 5 * (lenf - 1) - 1, ctx);

    return status;
}

int
_gr_poly_powmod_x_fmpz_preinv(gr_ptr res, const fmpz_t e, gr_srcptr f, slong lenf,
    gr_srcptr finv, slong lenfinv, gr_ctx_t ctx)
{
    gr_poly_preinv_t P;
    _gr_poly_preinv_init_newton_shallow(P, f, lenf, finv, lenfinv, ctx);
    return _gr_poly_preinv_powmod_x_fmpz(res, e, P, ctx);
}

int
gr_poly_preinv_powmod_x_fmpz(gr_poly_t res, const fmpz_t e, const gr_poly_preinv_t P, gr_ctx_t ctx)
{
    slong lenf = P->lenf;
    int status = GR_SUCCESS;
    gr_poly_t tmp;

    if (lenf == 0)
        return GR_DOMAIN;

    if (lenf == 1)
        return gr_poly_zero(res, ctx);

    if (fmpz_sgn(e) < 0)
        return GR_UNABLE;

    if (lenf == 2 || fmpz_cmp_ui(e, 2) <= 0)
    {
        gr_poly_init(tmp, ctx);
        status |= gr_poly_gen(tmp, ctx);
        status |= gr_poly_preinv_powmod_fmpz_binexp(res, tmp, e, P, ctx);
        gr_poly_clear(tmp, ctx);
        return status;
    }

    gr_poly_fit_length(res, lenf - 1, ctx);
    status |= _gr_poly_preinv_powmod_x_fmpz(res->coeffs, e, P, ctx);
    _gr_poly_set_length_normalise(res, lenf - 1, ctx);
    return status;
}

int
gr_poly_powmod_x_fmpz_preinv(gr_poly_t res,
                                        const fmpz_t e,
                                        const gr_poly_t f,
                                        const gr_poly_t finv,
                                        gr_ctx_t ctx)
{
    slong lenf = f->length;
    slong trunc = lenf - 1;
    gr_poly_t tmp;
    int status = GR_SUCCESS;

    if (lenf == 0)
        return GR_DOMAIN;

    if (fmpz_sgn(e) < 0)
        return GR_UNABLE;

    if (lenf == 1)
        return gr_poly_zero(res, ctx);

    if (lenf == 2)
    {
        gr_poly_t r, poly;
        gr_poly_init(tmp, ctx);
        gr_poly_init(r, ctx);
        gr_poly_init2(poly, 2, ctx);
        status |= gr_poly_gen(poly, ctx);
        status |= gr_poly_divrem(tmp, r, poly, f, ctx);
        status |= gr_poly_powmod_fmpz_binexp_preinv(res, r, e, f, finv, ctx);
        gr_poly_clear(tmp, ctx);
        gr_poly_clear(r, ctx);
        gr_poly_clear(poly, ctx);
        return status;
    }

    if (fmpz_is_zero(e))
    {
        return gr_poly_one(res, ctx);
    }
    else if (fmpz_is_one(e))
    {
        gr_poly_t r;
        gr_poly_init2(r, 2, ctx);
        gr_poly_init(tmp, ctx);
        status |= gr_poly_gen(r, ctx);
        status |= gr_poly_divrem(tmp, res, r, f, ctx);
        gr_poly_clear(tmp, ctx);
        gr_poly_clear(r, ctx);
        return status;
    }
    else if (*e == WORD(2))
    {
        gr_poly_init2(tmp, 2, ctx);
        status |= gr_poly_gen(tmp, ctx);
        status |= gr_poly_mulmod(res, tmp, tmp, f, ctx);
        gr_poly_clear(tmp, ctx);
        return status;
    }

    if ((res == f) || (res == finv))
    {
        gr_poly_init2(tmp, trunc, ctx);
        status |= _gr_poly_powmod_x_fmpz_preinv(tmp->coeffs, e, f->coeffs,
                                                 lenf, finv->coeffs,
                                                 finv->length, ctx);
        gr_poly_swap(res, tmp, ctx);
        gr_poly_clear(tmp, ctx);
    }
    else
    {
        gr_poly_fit_length(res, trunc, ctx);
        status |= _gr_poly_powmod_x_fmpz_preinv(res->coeffs, e, f->coeffs,
                                                 lenf, finv->coeffs,
                                                 finv->length, ctx);
    }

    _gr_poly_set_length_normalise(res, trunc, ctx);

    return status;
}
