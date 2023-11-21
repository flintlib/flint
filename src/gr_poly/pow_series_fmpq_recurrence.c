/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "fmpq.h"
#include "gr_vec.h"
#include "gr_poly.h"

/*
Implements the recurrence
https://en.wikipedia.org/wiki/Formal_power_series#Operations_on_formal_power_series
*/

int
_gr_poly_pow_series_fmpq_recurrence(gr_ptr h, gr_srcptr f, slong flen, const fmpq_t g, slong len, int precomp, gr_ctx_t ctx)
{
    gr_ptr a, b, s, t;
    slong i, l, alloc;
    slong sz = ctx->sizeof_elem;
    int use_divexact;
    int precomp_constant_term;
    int precomp_reciprocals;
    int status = GR_SUCCESS;

    flen = FLINT_MIN(flen, len);

    precomp_constant_term = ((precomp & 1) != 0);
    precomp_reciprocals   = ((precomp & 2) != 0);

    if (!precomp_constant_term)
    {
        status |= gr_pow_fmpq(h, f, g, ctx);

        if (status != GR_SUCCESS)
            return status;
    }

    if (flen == 1)
        return _gr_vec_zero(GR_ENTRY(h, 1, sz), len - 1, ctx);

    use_divexact = (fmpz_is_one(fmpq_denref(g)) &&
                (gr_ctx_is_integral_domain(ctx) == T_TRUE) &&
                (gr_ctx_is_finite_characteristic(ctx) == T_FALSE));

    alloc = 2 * flen + 2;

    GR_TMP_INIT_VEC(a, alloc, ctx);
    b = GR_ENTRY(a, flen, sz);
    s = GR_ENTRY(b, flen, sz);
    t = GR_ENTRY(s, 1, sz);

    /* b[k] = f[k] * q */
    if (fmpz_is_one(fmpq_denref(g)))
        status |= _gr_vec_set(b, f, flen, ctx);
    /* else if (*fmpq_denref(g) == 2)        --  todo: implement
        status |= _gr_vec_mul_two(b, f, flen, ctx); */
    else
        status |= _gr_vec_mul_scalar_fmpz(b, f, flen, fmpq_denref(g), ctx);

    /* a[k-1] = f[k] * (k * p) */
    status |= _gr_poly_derivative(a, f, flen, ctx);

    if (!fmpz_is_one(fmpq_numref(g)))
    {
        if (*fmpq_numref(g) == -1)
            status |= _gr_vec_neg(a, a, flen - 1, ctx);
        else
            status |= _gr_vec_mul_scalar_fmpz(a, a, flen - 1, fmpq_numref(g), ctx);
    }

    if (precomp_reciprocals)
        status |= gr_inv(b, b, ctx);

    for (i = 1; i < len && status == GR_SUCCESS; i++)
    {
        l = FLINT_MIN(i, flen - 1);

        status |= _gr_vec_sub(GR_ENTRY(a, 0, sz), GR_ENTRY(a, 0, sz), GR_ENTRY(b, 1, sz), FLINT_MIN(i, flen) - 1, ctx);
        status |= _gr_vec_dot_rev(s, NULL, 0, GR_ENTRY(a, 0, sz), GR_ENTRY(h, i - l, sz), l, ctx);

        /* h[i] = s / (i * (q * f[0])) */
        if (!precomp_reciprocals)
        {
            status |= gr_mul_ui(t, b, i, ctx);

            if (use_divexact)
                status |= gr_divexact(GR_ENTRY(h, i, sz), s, t, ctx);
            else
                status |= gr_div(GR_ENTRY(h, i, sz), s, t, ctx);
        }
        else
        {
            /* h[i] = s * (i^(-1) * (q * f[0])^(-1)) */
            status |= gr_mul(t, b, GR_ENTRY(h, i, sz), ctx);
            status |= gr_mul(GR_ENTRY(h, i, sz), s, t, ctx);
        }
    }

    GR_TMP_CLEAR_VEC(a, alloc, ctx);

    return status;
}

int
gr_poly_pow_series_fmpq_recurrence(gr_poly_t res,
    const gr_poly_t poly, const fmpq_t exp, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong flen;

    flen = poly->length;
    len = FLINT_MAX(len, 0);
    flen = FLINT_MIN(flen, len);

    if (fmpq_is_zero(exp))
    {
        if (len == 0)
            return gr_poly_zero(res, ctx);
        else
            return gr_poly_one(res, ctx);
    }
    else if (flen == 0)
    {
        if (fmpz_sgn(fmpq_numref(exp)) < 0)
            return GR_DOMAIN;
        else
            return gr_poly_zero(res, ctx);
    }
    else
    {
        if (flen == 1)
            len = 1;
        else if (fmpz_is_one(fmpq_denref(exp)) && !COEFF_IS_MPZ(*fmpq_numref(exp)) && *fmpq_numref(exp) > 0)
            len = poly_pow_length(flen, *fmpq_numref(exp), len);

        if (res != poly)
        {
            gr_poly_fit_length(res, len, ctx);
            status |= _gr_poly_pow_series_fmpq_recurrence(res->coeffs, poly->coeffs, flen, exp, len, 0, ctx);
            _gr_poly_set_length(res, len, ctx);
            _gr_poly_normalise(res, ctx);
        }
        else
        {
            gr_poly_t t;
            gr_poly_init2(t, len, ctx);
            status |= _gr_poly_pow_series_fmpq_recurrence(t->coeffs, poly->coeffs, flen, exp, len, 0, ctx);
            _gr_poly_set_length(t, len, ctx);
            _gr_poly_normalise(t, ctx);
            gr_poly_swap(res, t, ctx);
            gr_poly_clear(t, ctx);
        }

        return status;
    }
}
