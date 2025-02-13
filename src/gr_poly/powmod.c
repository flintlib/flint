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
static int _gr_poly_sqr(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
{
    return _gr_poly_mul(res, poly, len, poly, len, ctx);
}

/* todo: all of this could be implemented more simply by constructing
   a fixed-length polmod ring and using generic powering functions */

/* todo: preinversions should be automatic */

int
_gr_poly_powmod_fmpz_binexp(
    gr_ptr res,
    gr_srcptr poly,
    const fmpz_t e,
    gr_srcptr f, slong lenf,
    gr_ctx_t ctx)
{
    gr_ptr T, Q;
    slong lenT, lenQ;
    slong i;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (lenf == 2)
        return gr_pow_fmpz(res, poly, e, ctx);

    lenT = 2 * lenf - 3;
    lenQ = lenT - lenf + 1;

    GR_TMP_INIT_VEC(T, lenT + lenQ, ctx);
    Q = GR_ENTRY(T, lenT, sz);

    /* todo: should use divrem_preinv1 */
    /* status |= gr_inv(invf, f + lenf - 1, ctx); */

    status |= _gr_vec_set(res, poly, lenf - 1, ctx);

    for (i = fmpz_sizeinbase(e, 2) - 2; i >= 0; i--)
    {
        status |= _gr_poly_sqr(T, res, lenf - 1, ctx);
        status |= _gr_poly_divrem(Q, res, T, 2 * lenf - 3, f, lenf, ctx);

        if (fmpz_tstbit(e, i))
        {
            status |= _gr_poly_mul(T, res, lenf - 1, poly, lenf - 1, ctx);
            status |= _gr_poly_divrem(Q, res, T, 2 * lenf - 3, f, lenf, ctx);
        }
    }

    GR_TMP_CLEAR_VEC(T, lenT + lenQ, ctx);

    return status;
}

int
gr_poly_powmod_fmpz_binexp(gr_poly_t res,
                                      const gr_poly_t poly,
                                      const fmpz_t e,
                                      const gr_poly_t f,
                                      gr_ctx_t ctx)
{
    gr_ptr q;
    slong len = poly->length;
    slong lenf = f->length;
    slong trunc = lenf - 1;
    int qcopy = 0;
    int status = GR_SUCCESS;

    if (lenf == 0)
        return GR_DOMAIN;

    /* Not implemented. */
    if (fmpz_sgn(e) < 0)
        return GR_UNABLE;

    if (len >= lenf)
    {
        gr_poly_t t, r;
        gr_poly_init(t, ctx);
        gr_poly_init(r, ctx);
        status = gr_poly_divrem(t, r, poly, f, ctx);
        if (status == GR_SUCCESS)
            status |= gr_poly_powmod_fmpz_binexp(res, r, e, f, ctx);
        gr_poly_clear(t, ctx);
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
    else if (*e == WORD(2))
        return gr_poly_mulmod(res, poly, poly, f, ctx);

    if (poly->length < trunc)
    {
        GR_TMP_INIT_VEC(q, trunc, ctx);
        status |= _gr_vec_set(q, poly->coeffs, len, ctx);
        qcopy = 1;
    }
    else
    {
        q = poly->coeffs;
    }

    if ((res == poly && !qcopy) || (res == f))
    {
        gr_poly_t t;
        gr_poly_init2(t, 2 * lenf - 3, ctx);
        status |= _gr_poly_powmod_fmpz_binexp(t->coeffs, q, e, f->coeffs,
                                               lenf, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, 2 * lenf - 3, ctx);
        status |= _gr_poly_powmod_fmpz_binexp(res->coeffs, q, e, f->coeffs,
                                               lenf, ctx);
    }

    if (qcopy)
        GR_TMP_CLEAR_VEC(q, trunc, ctx);

    _gr_poly_set_length(res, trunc, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
_gr_poly_powmod_fmpz_binexp_preinv(
    gr_ptr res,
    gr_srcptr poly,
    const fmpz_t e,
    gr_srcptr f, slong lenf,
    gr_srcptr finv, slong lenfinv,
    gr_ctx_t ctx)
{
    gr_ptr T, Q;
    slong lenT, lenQ;
    slong i;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (lenf == 2)
        return gr_pow_fmpz(res, poly, e, ctx);

    lenT = 2 * lenf - 3;
    lenQ = lenT - lenf + 1;

    GR_TMP_INIT_VEC(T, lenT + lenQ, ctx);
    Q = GR_ENTRY(T, lenT, sz);

    status |= _gr_vec_set(res, poly, lenf - 1, ctx);

    for (i = fmpz_sizeinbase(e, 2) - 2; i >= 0; i--)
    {
        status |= _gr_poly_sqr(T, res, lenf - 1, ctx);
        status |= _gr_poly_divrem_newton_n_preinv(Q, res, T, 2 * lenf - 3, f,
                                                   lenf, finv, lenfinv, ctx);

        if (fmpz_tstbit(e, i))
        {
            status |= _gr_poly_mul(T, res, lenf - 1, poly, lenf - 1, ctx);
            status |= _gr_poly_divrem_newton_n_preinv(Q, res, T, 2 * lenf - 3,
                                                       f, lenf, finv, lenfinv,
                                                       ctx);
        }
    }

    GR_TMP_CLEAR_VEC(T, lenT + lenQ, ctx);

    return status;
}


int
gr_poly_powmod_fmpz_binexp_preinv(gr_poly_t res,
                                             const gr_poly_t poly,
                                             const fmpz_t e,
                                             const gr_poly_t f,
                                             const gr_poly_t finv,
                                             gr_ctx_t ctx)
{
    gr_ptr q;
    slong len = poly->length;
    slong lenf = f->length;
    slong trunc = lenf - 1;
    int qcopy = 0;
    int status = GR_SUCCESS;

    if (lenf == 0)
        return GR_DOMAIN;

    /* Not implemented. */
    if (fmpz_sgn(e) < 0)
        return GR_UNABLE;

    if (len >= lenf)
    {
        gr_poly_t t, r;
        gr_poly_init(t, ctx);
        gr_poly_init(r, ctx);
        status = gr_poly_divrem(t, r, poly, f, ctx);  /* todo: use preinv */
        if (status == GR_SUCCESS)
            status |= gr_poly_powmod_fmpz_binexp_preinv(res, r, e, f, finv, ctx);
        gr_poly_clear(t, ctx);
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
    else if (*e == WORD(2))
        return gr_poly_mulmod(res, poly, poly, f, ctx);

    if (poly->length < trunc)
    {
        GR_TMP_INIT_VEC(q, trunc, ctx);
        status |= _gr_vec_set(q, poly->coeffs, len, ctx);
        qcopy = 1;
    }
    else
    {
        q = poly->coeffs;
    }

    if ((res == poly && !qcopy) || (res == f) || (res == finv))
    {
        gr_poly_t t;
        gr_poly_init2(t, 2 * lenf - 3, ctx);
        status |= _gr_poly_powmod_fmpz_binexp_preinv(t->coeffs, q, e, f->coeffs,
                                               lenf, finv->coeffs, finv->length, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, 2 * lenf - 3, ctx);
        status |= _gr_poly_powmod_fmpz_binexp_preinv(res->coeffs, q, e, f->coeffs,
                                               lenf, finv->coeffs, finv->length, ctx);
    }

    if (qcopy)
        GR_TMP_CLEAR_VEC(q, trunc, ctx);

    _gr_poly_set_length(res, trunc, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}


int
_gr_poly_powmod_ui_binexp(gr_ptr res,
                                     gr_srcptr poly, ulong e,
                                     gr_srcptr f, slong lenf,
                                     gr_ctx_t ctx)
{
    fmpz_t e2;
    int status;
    fmpz_init_set_ui(e2, e);
    status = _gr_poly_powmod_fmpz_binexp(res, poly, e2, f, lenf, ctx);
    fmpz_clear(e2);
    return status;
}

int
gr_poly_powmod_ui_binexp(gr_poly_t res,
                                    const gr_poly_t poly, ulong e,
                                    const gr_poly_t f,
                                    gr_ctx_t ctx)
{
    fmpz_t e2;
    int status;
    fmpz_init_set_ui(e2, e);
    status = gr_poly_powmod_fmpz_binexp(res, poly, e2, f, ctx);
    fmpz_clear(e2);
    return status;
}

int
_gr_poly_powmod_ui_binexp_preinv(
    gr_ptr res,
    gr_srcptr poly,
    ulong e,
    gr_srcptr f, slong lenf,
    gr_srcptr finv, slong lenfinv,
    gr_ctx_t ctx)
{
    fmpz_t e2;
    int status;
    fmpz_init_set_ui(e2, e);
    status = _gr_poly_powmod_fmpz_binexp_preinv(res, poly, e2, f, lenf, finv, lenfinv, ctx);
    fmpz_clear(e2);
    return status;
}

int
gr_poly_powmod_ui_binexp_preinv(gr_poly_t res,
                                           const gr_poly_t poly,
                                           ulong e,
                                           const gr_poly_t f,
                                           const gr_poly_t finv,
                                           gr_ctx_t ctx)
{
    fmpz_t e2;
    int status;
    fmpz_init_set_ui(e2, e);
    status = gr_poly_powmod_fmpz_binexp_preinv(res, poly, e2, f, finv, ctx);
    fmpz_clear(e2);
    return status;
}

int
_gr_poly_powmod_fmpz_sliding_preinv(
    gr_ptr res,
    gr_srcptr poly,
    const fmpz_t e, ulong k,
    gr_srcptr f, slong lenf,
    gr_srcptr finv, slong lenfinv,
    gr_ctx_t ctx)
{
    gr_ptr T, Q;
    gr_poly_struct * precomp;
    gr_poly_t poly_squared;
    ulong twokm1;
    slong lenT, lenQ;
    slong i, l, j;
    int index;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (lenf == 2)
        return gr_pow_fmpz(res, poly, e, ctx);

    lenT = 2 * lenf - 3;
    lenQ = lenT - lenf + 1;

    GR_TMP_INIT_VEC(T, lenT + lenQ, ctx);
    Q = GR_ENTRY(T, lenT, sz);

    /* Precomputation */
    twokm1 = n_pow(2, k - 1);
    precomp = flint_malloc(twokm1 * sizeof(gr_poly_struct));
    gr_poly_init(precomp, ctx);
    gr_poly_fit_length(precomp, lenf - 1, ctx);
    status |= _gr_vec_set(precomp->coeffs, poly, lenf - 1, ctx);

    gr_poly_init(poly_squared, ctx);
    if (k > 1)
    {
        gr_poly_fit_length(poly_squared, lenf - 1, ctx);
        status |= _gr_poly_mul(T, poly, lenf - 1, poly, lenf - 1, ctx);
        status |= _gr_poly_divrem_newton_n_preinv(Q, poly_squared->coeffs, T,
                                                   2 * lenf - 3, f, lenf, finv,
                                                   lenfinv, ctx);
    }
    for (i = 1; (ulong) i < twokm1; i++)
    {
        gr_poly_init(precomp + i, ctx);
        gr_poly_fit_length(precomp + i, lenf - 1, ctx);
        status |= _gr_poly_mul(T, (precomp + i - 1)->coeffs, lenf - 1,
                                poly_squared->coeffs, lenf - 1, ctx);
        status |= _gr_poly_divrem_newton_n_preinv(Q, (precomp + i)->coeffs, T,
                                                   2 * lenf - 3, f, lenf, finv,
                                                   lenfinv, ctx);
    }

    status |= _gr_vec_set(res, poly, lenf - 1, ctx);

    i = fmpz_sizeinbase(e, 2) - 2;
    while (i >= 0)
    {
        if (fmpz_tstbit(e, i) == 0)
        {
            status |= _gr_poly_sqr(T, res, lenf - 1, ctx);
            status |= _gr_poly_divrem_newton_n_preinv(Q, res, T, 2 * lenf - 3,
                                                       f, lenf, finv, lenfinv,
                                                       ctx);
            i -= 1;
        }
        else
        {
            l = FLINT_MAX(i - k + 1, 0);
            while (fmpz_tstbit(e, l) == 0)
            {
                l += 1;
            }
            for (j = 0; j < i - l + 1; j++)
            {
                status |= _gr_poly_sqr(T, res, lenf - 1, ctx);
                status |= _gr_poly_divrem_newton_n_preinv(Q, res, T,
                                                           2 * lenf - 3, f,
                                                           lenf, finv, lenfinv,
                                                           ctx);
            }

            index = fmpz_tstbit(e, i);
            for (j = i - 1; j >= l; j--)
            {
                index = index << 1;
                index += fmpz_tstbit(e, j);
            }
            index = (index - 1) / 2;

            status |= _gr_poly_mul(T, res, lenf - 1,
                                    (precomp + index)->coeffs, lenf - 1, ctx);
            status |= _gr_poly_divrem_newton_n_preinv(Q, res, T, 2 * lenf - 3,
                                                       f, lenf, finv, lenfinv,
                                                       ctx);
            i = l - 1;
        }
    }

    for (j = 0; (ulong) j < twokm1; j++)
    {
        gr_poly_clear(precomp + j, ctx);
    }
    flint_free(precomp);
    gr_poly_clear(poly_squared, ctx);

    GR_TMP_CLEAR_VEC(T, lenT + lenQ, ctx);

    return status;
}


int
gr_poly_powmod_fmpz_sliding_preinv(gr_poly_t res,
                                              const gr_poly_t poly,
                                              const fmpz_t e, ulong k,
                                              const gr_poly_t f,
                                              const gr_poly_t finv,
                                              gr_ctx_t ctx)
{
    gr_ptr q;
    slong len = poly->length;
    slong lenf = f->length;
    slong trunc = lenf - 1;
    flint_bitcnt_t bits;
    int qcopy = 0;
    int status = GR_SUCCESS;

    if (lenf == 0)
        return GR_DOMAIN;

    /* Not implemented. */
    if (fmpz_sgn(e) < 0)
        return GR_UNABLE;

    if (len >= lenf)
    {
        gr_poly_t t, r;
        gr_poly_init(t, ctx);
        gr_poly_init(r, ctx);
        status = gr_poly_divrem(t, r, poly, f, ctx);
        if (status == GR_SUCCESS)
            status |= gr_poly_powmod_fmpz_sliding_preinv(res, r, e, k, f, finv, ctx);
        gr_poly_clear(t, ctx);
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
    else if (*e == WORD(2))
        return gr_poly_mulmod(res, poly, poly, f, ctx);

    if (poly->length < trunc)
    {
        GR_TMP_INIT_VEC(q, trunc, ctx);
        status |= _gr_vec_set(q, poly->coeffs, len, ctx);
        qcopy = 1;
    }
    else
    {
        q = poly->coeffs;
    }

    /* Determine "optimum" sliding window size */
    if (k == 0)
    {
        bits = fmpz_bits(e);
        if (bits < 9)
            k = 1;
        else if (bits < 15)
            k = 2;
        else if (bits < 62)
            k = 3;
        else if (bits < 203)
            k = 4;
        else if (bits < 587)
            k = 5;
        else if (bits < 1560)
            k = 6;
        else
            k = 7;
    }

    if ((res == poly && !qcopy) || (res == f) || (res == finv))
    {
        gr_poly_t t;
        gr_poly_init2(t, 2 * lenf - 3, ctx);
        status |= _gr_poly_powmod_fmpz_sliding_preinv(t->coeffs, q, e, k, f->coeffs,
                                               lenf, finv->coeffs, finv->length, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, 2 * lenf - 3, ctx);
        status |= _gr_poly_powmod_fmpz_sliding_preinv(res->coeffs, q, e, k, f->coeffs,
                                               lenf, finv->coeffs, finv->length, ctx);
    }

    if (qcopy)
        GR_TMP_CLEAR_VEC(q, trunc, ctx);

    _gr_poly_set_length(res, trunc, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}


int
_gr_poly_powmod_x_fmpz_preinv(
    gr_ptr res,
    const fmpz_t e,
    gr_srcptr f, slong lenf,
    gr_srcptr finv, slong lenfinv,
    gr_ctx_t ctx)
{
    gr_ptr T, Q;
    slong lenT, lenQ;
    slong i, window, l, c;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    lenT = 2 * lenf - 3;
    lenQ = lenT - lenf + 1;

    GR_TMP_INIT_VEC(T, lenT + lenQ, ctx);
    Q = GR_ENTRY(T, lenT, sz);

    status |= gr_one(res, ctx);
    status |= _gr_vec_zero(GR_ENTRY(res, 1, sz), lenf - 2, ctx);

    l = z_sizeinbase(lenf - 1, 2) - 2;
    window = 0;
    window = (WORD(1) << l);
    c = l;
    i = fmpz_sizeinbase(e, 2) - 2;

    if (i <= l)
    {
        window = 0;
        window = (1 << i);
        c = i;
        l = i;
    }

    if (c == 0)
    {
        status |= _gr_poly_shift_left(T, res, lenf - 1, window, ctx);
        status |= _gr_poly_divrem_newton_n_preinv(Q, res, T,
                                                   lenf - 1 + window, f, lenf,
                                                   finv, lenfinv, ctx);
        c = l + 1;
        window = 0;
    }

    for (; i >= 0; i--)
    {
        status |= _gr_poly_sqr(T, res, lenf - 1, ctx);
        status |= _gr_poly_divrem_newton_n_preinv(Q, res, T, 2 * lenf - 3, f,
                                                   lenf, finv, lenfinv, ctx);

        c--;
        if (fmpz_tstbit(e, i))
        {
            if (window == 0 && i <= l - 1)
                c = i;
            if (c >= 0)
                window = window | (1 << c);
        }
        else if (window == 0)
        {
            c = l + 1;
        }

        if (c == 0)
        {
            status |= _gr_poly_shift_left(T, res, lenf - 1, window, ctx);
            status |= _gr_poly_divrem_newton_n_preinv(Q, res, T,
                                                       lenf - 1 + window, f,
                                                       lenf, finv, lenfinv,
                                                       ctx);
            c = l + 1;
            window = 0;
        }
    }

    GR_TMP_CLEAR_VEC(T, lenT + lenQ, ctx);

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

    _gr_poly_set_length(res, trunc, ctx);
    _gr_poly_normalise(res, ctx);

    return status;
}
