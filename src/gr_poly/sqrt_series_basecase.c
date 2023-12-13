/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* todo: with a slight rewrite, could allow aliasing */
/* todo: use _gr_poly_pow_series_fmpq_recurrence if possible, when flen is short */
int
_gr_poly_sqrt_series_basecase(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    flen = FLINT_MIN(flen, len);

    status |= gr_sqrt(res, f, ctx);

    if (status != GR_SUCCESS)
        return status;

    if (flen == 1)
    {
        status |= _gr_vec_zero(GR_ENTRY(res, 1, sz), len - 1, ctx);
    }
    else if (len == 2)
    {
        status |= gr_mul(GR_ENTRY(res, 1, sz), res, GR_ENTRY(f, 1, sz), ctx);
        status |= gr_div(GR_ENTRY(res, 1, sz), GR_ENTRY(res, 1, sz), GR_ENTRY(f, 0, sz), ctx);
        status |= gr_mul_2exp_si(GR_ENTRY(res, 1, sz), GR_ENTRY(res, 1, sz), -1, ctx);
    }
    else
    {
        int is_one, have_inv;
        gr_ptr rinv = NULL;
        slong i, l;

        is_one = (gr_is_one(res, ctx) == T_TRUE);

        if (!is_one)
        {
            GR_TMP_INIT(rinv, ctx);
            have_inv = (gr_inv(rinv, res, ctx) == GR_SUCCESS);
        }

        for (i = 1; i < len; i++)
        {
            gr_srcptr initial = GR_ENTRY(res, i, sz);

            l = (i - 1) / 2;

            if (i % 2 == 1)
            {
                if (i < flen)
                    status |= gr_mul_2exp_si(GR_ENTRY(res, i, sz), GR_ENTRY(f, i, sz), -1, ctx);
                else
                    initial = NULL;
            }
            else
            {
                status |= gr_sqr(GR_ENTRY(res, i, sz), GR_ENTRY(res, i / 2, sz), ctx);

                if (i < flen)
                    status |= gr_sub(GR_ENTRY(res, i, sz), GR_ENTRY(f, i, sz), GR_ENTRY(res, i, sz), ctx);
                else
                    status |= gr_neg(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), ctx);

                status |= gr_mul_2exp_si(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), -1, ctx);
            }

            if (status != GR_SUCCESS)
                break;

            status |= _gr_vec_dot_rev(GR_ENTRY(res, i, sz), initial, 1, GR_ENTRY(res, 1, sz), GR_ENTRY(res, i - l, sz), l, ctx);

            if (!is_one)
            {
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
                if (have_inv)
                {
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC diagnostic pop
#endif
                    status |= gr_mul(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), rinv, ctx);
                }
                else
                {
                    status |= gr_div(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), res, ctx);
                    if (status != GR_SUCCESS)
                        break;
                }
            }
        }

        if (!is_one)
        {
            GR_TMP_CLEAR(rinv, ctx);
        }
    }

    return status;
}

int
gr_poly_sqrt_series_basecase(gr_poly_t res, const gr_poly_t h, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong hlen;

    hlen = h->length;

    if (hlen == 0 || len == 0)
        return gr_poly_zero(res, ctx);

    if (hlen == 1)
        len = 1;

    if (res == h)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status = gr_poly_sqrt_series_basecase(t, h, len, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_sqrt_series_basecase(res->coeffs, h->coeffs, h->length, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
