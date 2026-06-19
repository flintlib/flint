/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_special.h"

/* Simultaneous Newton iteration for S = sin(h) and C = cos(h), derived
   from the Newton iteration applied to F = exp(ih) = C + iS,
   G = exp(-ih) = C - iS, FG = 1. */
int
_gr_poly_sin_cos_series_newton(gr_ptr S, gr_ptr C,
    gr_srcptr h, slong hlen, slong len, slong cutoff,
    int times_pi, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong i, j, m, n, l, r, nm, Alen;
    gr_ptr hprime, A, B, P, Q, tmp, recip = NULL;
    int use_reciprocals = 0;
    gr_ptr S_tmp = NULL, C_tmp = NULL;
    slong a[FLINT_BITS];
    int want_S = (S != NULL);
    int want_C = (C != NULL);

    hlen = FLINT_MIN(hlen, len);

    /* Allocate scratch for the unwanted output if NULL was passed. */
    if (!want_S)
    {
        GR_TMP_INIT_VEC(S_tmp, len, ctx);
        S = S_tmp;
    }
    if (!want_C)
    {
        GR_TMP_INIT_VEC(C_tmp, len, ctx);
        C = C_tmp;
    }

    if (len < cutoff || hlen <= 1)
    {
        status |= _gr_poly_sin_cos_series_basecase(S, C, h, hlen, len, times_pi, ctx);
        goto cleanup;
    }

    /* Must be computed before basecase in case of aliasing */
    GR_TMP_INIT_VEC(hprime, hlen - 1, ctx);
    status |= _gr_poly_derivative(hprime, h, hlen, ctx);

    cutoff = FLINT_MAX(cutoff, 2);

    a[i = 0] = n = len;
    while (n >= cutoff)
        a[++i] = (n = (n + 1) / 2);

    status |= _gr_poly_sin_cos_series_basecase(S, C, h, hlen, n, times_pi, ctx);

    slong half = len / 2 + 1;

    GR_TMP_INIT_VEC(A,      half,     ctx);
    GR_TMP_INIT_VEC(B,      half,     ctx);
    GR_TMP_INIT_VEC(P,      half,     ctx);
    GR_TMP_INIT_VEC(Q,      half,     ctx);
    GR_TMP_INIT_VEC(tmp,    half,     ctx);

    use_reciprocals = (gr_ctx_is_finite_characteristic(ctx) == T_TRUE);
    if (use_reciprocals)
    {
        GR_TMP_INIT_VEC(recip, len - 1, ctx);
        status |= _gr_vec_reciprocals(recip, len - 1, ctx);
    }

    if (times_pi)
    {
        gr_ptr pi;
        GR_TMP_INIT(pi, ctx);
        status |= gr_pi(pi, ctx);
        status |= _gr_vec_mul_scalar(hprime, hprime, hlen - 1, pi, ctx);
        GR_TMP_CLEAR(pi, ctx);
    }

    /* Always use Karatsuba multiplication? Should generally be faster;
       drawback is somewhat increased precision loss with numerical
       coefficients. May want disabling for low precision, huge len. */
    int use_karatsuba = 1;

    for (i--; i >= 0; i--)
    {
        m = n;
        n = a[i];
        nm = n - m;

        l = FLINT_MIN(hlen - 1, n);
        r = FLINT_MIN(l + m - 1, n - 1); /* mulmid extracts positions m-1..r-1 */
        Alen = r - m + 1;                /* number of terms in A and B */

        /* Compute A = [h'C_m]_{m-1..r-1}  and  B = [h'S_m]_{m-1..r-1} */
        status |= _gr_poly_mulmid(A, hprime, l, C, m, m - 1, r, ctx);
        status |= _gr_poly_mulmid(B, hprime, l, S, m, m - 1, r, ctx);

        /* Compute P = A*S_m - B*C_m  and  Q = A*C_m + B*S_m */
        if (!use_karatsuba)
        {
            status |= _gr_poly_mullow(P,   A, Alen, S, nm, nm, ctx);
            status |= _gr_poly_mullow(tmp, B, Alen, C, nm, nm, ctx);
            status |= _gr_vec_sub(P, P, tmp, nm, ctx);
            status |= _gr_poly_mullow(Q,   A, Alen, C, nm, nm, ctx);
            status |= _gr_poly_mullow(tmp, B, Alen, S, nm, nm, ctx);
            status |= _gr_vec_add(Q, Q, tmp, nm, ctx);
        }
        else
        {
            status |= _gr_poly_mullow(P,   B, Alen, C, nm, nm, ctx);
            status |= _gr_vec_neg(P, P, nm, ctx);
            status |= _gr_poly_mullow(Q, A, Alen, S, nm, nm, ctx);
            status |= _gr_vec_neg(Q, Q, nm, ctx);
            status |= _gr_vec_sub(A, A, B, Alen, ctx);
            status |= _gr_vec_sub(tmp, C, S, nm, ctx);
            status |= _gr_poly_mullow(B, A, Alen, tmp, nm, nm, ctx);
            status |= _gr_vec_sub(B, B, P, nm, ctx);
            status |= _gr_vec_sub(B, B, Q, nm, ctx);
            status |= _gr_vec_sub(P, P, Q, nm, ctx);
            status |= _gr_vec_set(Q, B, nm, ctx);
        }

        /* Integrate P and Q */
        if (use_reciprocals)
        {
            status |= _gr_vec_mul(P, P, GR_ENTRY(recip, m - 1, sz), nm, ctx);
            status |= _gr_vec_mul(Q, Q, GR_ENTRY(recip, m - 1, sz), nm, ctx);
        }
        else
        {
            for (j = 0; j < nm; j++)
            {
                status |= gr_div_ui(GR_ENTRY(P, j, sz), GR_ENTRY(P, j, sz), m + j, ctx);
                status |= gr_div_ui(GR_ENTRY(Q, j, sz), GR_ENTRY(Q, j, sz), m + j, ctx);
            }
        }

        /* Update S and C */
        /* On the last Newton step (i == 0) when only one output is wanted,
           use the standard branch regardless of use_karatsuba so we can skip
           the two mullows for the unwanted output. */
        if (use_karatsuba && (i != 0 || (want_S && want_C)))
        {
            status |= _gr_poly_mullow(GR_ENTRY(C, m, sz), C, nm, P, nm, nm, ctx);
            status |= _gr_poly_mullow(GR_ENTRY(S, m, sz), S, nm, Q, nm, nm, ctx);
            status |= _gr_vec_add(tmp, C, S, nm, ctx);
            status |= _gr_vec_add(P, P, Q, nm, ctx);
            status |= _gr_poly_mullow(A, tmp, nm, P, nm, nm, ctx);
            status |= _gr_vec_sub(A, A, GR_ENTRY(C, m, sz), nm, ctx);
            status |= _gr_vec_sub(A, A, GR_ENTRY(S, m, sz), nm, ctx);
            status |= _gr_vec_sub(GR_ENTRY(C, m, sz),
                                  GR_ENTRY(C, m, sz), GR_ENTRY(S, m, sz), nm, ctx);
            status |= _gr_vec_set(GR_ENTRY(S, m, sz), A, nm, ctx);
        }
        else
        {
            if (want_C || i != 0)
            {
                status |= _gr_poly_mullow(GR_ENTRY(C, m, sz), C, nm, P, nm, nm, ctx);
                status |= _gr_poly_mullow(tmp, S, nm, Q, nm, nm, ctx);
                status |= _gr_vec_sub(GR_ENTRY(C, m, sz), GR_ENTRY(C, m, sz), tmp, nm, ctx);
            }

            if (want_S || i != 0)
            {
                status |= _gr_poly_mullow(GR_ENTRY(S, m, sz), S, nm, P, nm, nm, ctx);
                status |= _gr_poly_mullow(tmp, C, nm, Q, nm, nm, ctx);
                status |= _gr_vec_add(GR_ENTRY(S, m, sz), GR_ENTRY(S, m, sz), tmp, nm, ctx);
            }
        }
    }

    GR_TMP_CLEAR_VEC(hprime, hlen - 1, ctx);
    GR_TMP_CLEAR_VEC(A,      half,     ctx);
    GR_TMP_CLEAR_VEC(B,      half,     ctx);
    GR_TMP_CLEAR_VEC(P,      half,     ctx);
    GR_TMP_CLEAR_VEC(Q,      half,     ctx);
    GR_TMP_CLEAR_VEC(tmp,    half,     ctx);
    if (use_reciprocals)
        GR_TMP_CLEAR_VEC(recip, len - 1, ctx);

cleanup:
    if (S_tmp != NULL) { GR_TMP_CLEAR_VEC(S_tmp, len, ctx); }
    if (C_tmp != NULL) { GR_TMP_CLEAR_VEC(C_tmp, len, ctx); }

    return status;
}

int
gr_poly_sin_cos_series_newton(gr_poly_t s, gr_poly_t c,
                                    const gr_poly_t h, slong n, slong cutoff, int times_pi, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong hlen = h->length;

    if (n == 0)
    {
        status |= gr_poly_zero(s, ctx);
        status |= gr_poly_zero(c, ctx);
        return status;
    }

    if (hlen == 0)
    {
        status |= gr_poly_zero(s, ctx);
        status |= gr_poly_one(c, ctx);
        return status;
    }

    gr_poly_fit_length(s, n, ctx);
    gr_poly_fit_length(c, n, ctx);
    status |= _gr_poly_sin_cos_series_newton(s->coeffs, c->coeffs,
        h->coeffs, hlen, n, cutoff, times_pi, ctx);
    _gr_poly_set_length(s, n, ctx);
    _gr_poly_normalise(s, ctx);
    _gr_poly_set_length(c, n, ctx);
    _gr_poly_normalise(c, ctx);

    return status;
}

