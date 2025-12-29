/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"
#include "longlong.h"

/* compose by poly2 = a*x^n + c, no aliasing; n >= 1 */
static int
_gr_poly_compose_axnc(gr_ptr res, gr_srcptr poly1, slong len1,
        gr_srcptr c, gr_srcptr a, slong n, gr_ctx_t ctx)
{
    slong i, sz = ctx->sizeof_elem;
    int status;

    /* shift by c (c = 0 case will be fast) */
    status = _gr_poly_taylor_shift(res, poly1, len1, c, ctx);

    /* multiply by powers of a */
    if (gr_is_one(a, ctx) != T_TRUE)
    {
        if (gr_is_neg_one(a, ctx) == T_TRUE)
        {
            for (i = 1; i < len1; i += 2)
                status |= gr_neg(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), ctx);
        }
        else if (len1 == 2)
        {
            status |= gr_mul(GR_ENTRY(res, 1, sz), GR_ENTRY(res, 1, sz), a, ctx);
        }
        else
        {
            int maxbit = FLINT_CLOG2(len1);
            gr_ptr t;

            /* Prefer squaring for powers? cf. _gr_vec_set_powers */
            if (gr_ctx_is_finite(ctx) == T_TRUE || gr_ctx_has_real_prec(ctx) == T_TRUE)
            {
                GR_TMP_INIT_VEC(t, maxbit, ctx);

                status |= gr_set(GR_ENTRY(t, 0, sz), a, ctx);
                status |= gr_mul(GR_ENTRY(res, 1, sz), GR_ENTRY(res, 1, sz), a, ctx);
                for (int j = 1; j < maxbit; ++j)
                {
                    status |= gr_sqr(GR_ENTRY(t, j, sz), GR_ENTRY(t, j-1, sz), ctx);
                    status |= gr_mul(GR_ENTRY(res, 1<<j, sz), GR_ENTRY(res, 1<<j, sz), GR_ENTRY(t, j, sz), ctx);
                }
                for (i = ((slong) 1 << (maxbit-1)) + 1; i < len1; i++)
                {
                    int bit = flint_ctz(i);
                    status |= gr_mul(GR_ENTRY(t, maxbit-1-bit, sz), GR_ENTRY(t, maxbit-1-bit, sz), a, ctx);
                    status |= gr_mul(GR_ENTRY(res, i>>bit, sz), GR_ENTRY(res, i>>bit, sz), GR_ENTRY(t, maxbit-1-bit, sz), ctx);

                    for (int j = bit; j > 0; --j)
                    {
                        status |= gr_sqr(GR_ENTRY(t, maxbit-j, sz), GR_ENTRY(t, maxbit-j-1, sz), ctx);
                        status |= gr_mul(GR_ENTRY(res, i>>(j-1), sz), GR_ENTRY(res, i>>(j-1), sz), GR_ENTRY(t, maxbit-j, sz), ctx);
                    }
                }

                for (i = (len1 + 1) >> 1; i < ((slong)1 << (maxbit-1)); i++)
                {
                    int bit = flint_ctz(i);
                    status |= gr_mul(GR_ENTRY(t, maxbit-2-bit, sz), GR_ENTRY(t, maxbit-2-bit, sz), a, ctx);
                    status |= gr_mul(GR_ENTRY(res, i>>bit, sz), GR_ENTRY(res, i>>bit, sz), GR_ENTRY(t, maxbit-2-bit, sz), ctx);

                    for (int j = bit; j > 0; --j)
                    {
                        status |= gr_sqr(GR_ENTRY(t, maxbit-j-1, sz), GR_ENTRY(t, maxbit-j-2, sz), ctx);
                        status |= gr_mul(GR_ENTRY(res, i>>(j-1), sz), GR_ENTRY(res, i>>(j-1), sz), GR_ENTRY(t, maxbit-j-1, sz), ctx);
                    }
                }

                GR_TMP_CLEAR_VEC(t, maxbit, ctx);
            }
            else
            {
                GR_TMP_INIT(t, ctx);

                status |= gr_set(t, a, ctx);

                for (i = 1; i < len1; i++)
                {
                    status |= gr_mul(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), t, ctx);
                    if (i + 1 < len1)
                        status |= gr_mul(t, t, a, ctx);
                }

                GR_TMP_CLEAR(t, ctx);
            }
        }
    }

    /* stretch */
    for (i = len1 - 1; i >= 1 && n > 1; i--)
    {
        gr_swap(GR_ENTRY(res, i * n, sz), GR_ENTRY(res, i, sz), ctx);
        status |= _gr_vec_zero(GR_ENTRY(res, (i - 1) * n + 1, sz), n - 1, ctx);
    }

    return status;
}

int
_gr_poly_compose(gr_ptr res,
    gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, gr_ctx_t ctx)
{
    slong sz;

    if (len1 == 1)
        return gr_set(res, poly1, ctx);

    if (len2 == 1)
        return _gr_poly_evaluate(res, poly1, len1, poly2, ctx);

    sz = ctx->sizeof_elem;

    if (_gr_vec_is_zero(GR_ENTRY(poly2, 1, sz), len2 - 2, ctx) == T_TRUE)
        return _gr_poly_compose_axnc(res, poly1, len1, poly2, GR_ENTRY(poly2, len2 - 1, sz), len2 - 1, ctx);

    if (len1 <= 7)
        return _gr_poly_compose_horner(res, poly1, len1, poly2, len2, ctx);

    return _gr_poly_compose_divconquer(res, poly1, len1, poly2, len2, ctx);
}

int gr_poly_compose(gr_poly_t res,
    const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;

    if (len1 == 0)
    {
        return gr_poly_zero(res, ctx);
    }
    else if (len1 == 1 || len2 == 0)
    {
        return gr_poly_set_scalar(res, poly1->coeffs, ctx);
    }
    else
    {
        const slong lenr = (len1 - 1) * (len2 - 1) + 1;
        int status;

        if (res != poly1 && res != poly2)
        {
            gr_poly_fit_length(res, lenr, ctx);
            status = _gr_poly_compose(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, ctx);
        }
        else
        {
            gr_poly_t t;
            gr_poly_init2(t, lenr, ctx);
            status = _gr_poly_compose(t->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, ctx);
            gr_poly_swap(res, t, ctx);
            gr_poly_clear(t, ctx);
        }

        _gr_poly_set_length(res, lenr, ctx);
        _gr_poly_normalise(res, ctx);
        return status;
    }
}
