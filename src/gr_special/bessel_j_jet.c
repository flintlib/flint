/*
    Copyright (C) 2026 Joel Dahne

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "fmpz.h"

/*
    For z != 0, uses the recurrence from https://fungrim.org/entry/9b2f38/

    For z = 0 and integer nu = n, uses the series expansion

        J_n(z) = sum_{k=0}^{inf} (-1)^k / (k! * (k+n)! * 2^{2k+n}) * z^{2k+n}

    for non-negative n and J_{n}(z) = (-1)^n J_{-n}(z) for negative n.
*/
int
gr_bessel_j_jet(gr_ptr res, gr_srcptr nu, gr_srcptr z, slong len, gr_ctx_t ctx)
{
    gr_ptr nu_shifted, t, z2, z2_minus_nu2;
    fmpz_t tmp_fmpz;
    slong sz = ctx->sizeof_elem;
    slong k, n;
    int status = GR_SUCCESS;
    int is_integer;

    if (len <= 0)
        return status;

    /* TODO: This should use gr_is_integer */
    is_integer = gr_get_si(&n, nu, ctx) == GR_SUCCESS;

    /* TODO: Handle gr_is_zero(z, ctx) == T_UNKNOWN by using
     * https://fungrim.org/entry/2488BB/ and gr_poly_binomial_transform (which
     * is yet to be implemented) */
    if (gr_is_zero(z, ctx) == T_TRUE && is_integer)
    {
        slong n;
        if (gr_get_si(&n, nu, ctx) == GR_SUCCESS)
        {
            slong abs_n;

            abs_n = FLINT_ABS(n);

            /* Many of the terms are zero. We zero the entire vector
             * and adjust the non-zero ones. */
            status |= _gr_vec_zero(res, len, ctx);

            if (abs_n < len)
            {
                fmpz_init(tmp_fmpz);

                /* Start: res[abs_n] = 1 / (abs_n! * 2^abs_n) */
                gr_ptr c_start = GR_ENTRY(res, abs_n, sz);
                status |= gr_fac_ui(c_start, abs_n, ctx);
                status |= gr_mul_2exp_si(c_start, c_start, abs_n, ctx);
                status |= gr_inv(c_start, c_start, ctx);

                for (k = 1; 2 * k + abs_n < len; k++)
                {
                    fmpz_set_si(tmp_fmpz, -4 * k);
                    fmpz_mul_si(tmp_fmpz, tmp_fmpz, k + abs_n);

                    status |= gr_div_fmpz(
                        GR_ENTRY(res, abs_n + 2 * k, sz),
                        GR_ENTRY(res, abs_n + 2 * (k - 1), sz),
                        tmp_fmpz,
                        ctx
                        );
                }

                if ((n < 0) && (abs_n % 2 != 0))
                    status |= _gr_vec_neg(res, res, len, ctx);

                fmpz_clear(tmp_fmpz);
            }

            return status;
        }
        else
        {
            /* In this case nu is an integer, but doesn't fit in a
             * slong. It follows that len < abs(nu). Since the lowest
             * order term is of degree abs(nu) all terms up to len are
             * zero. */
            status |= _gr_vec_zero(res, len, ctx);
            return status;
        }
    }

    /* res[0] = J_{nu}(z) */
    status |= gr_bessel_j(GR_ENTRY(res, 0, sz), nu, z, ctx);

    if (len == 1)
        return status;

    /* res[1] = (J_{nu-1}(z) - J_{nu+1}(z)) / 2 */
    GR_TMP_INIT2(nu_shifted, t, ctx);

    status |= gr_sub_si(nu_shifted, nu, 1, ctx);
    status |= gr_bessel_j(GR_ENTRY(res, 1, sz), nu_shifted, z, ctx);
    status |= gr_add_si(nu_shifted, nu, 1, ctx);
    status |= gr_bessel_j(t, nu_shifted, z, ctx);
    status |= gr_sub(GR_ENTRY(res, 1, sz), GR_ENTRY(res, 1, sz), t, ctx);
    status |= gr_mul_2exp_si(GR_ENTRY(res, 1, sz), GR_ENTRY(res, 1, sz), -1, ctx);

    if (len == 2)
    {
        GR_TMP_CLEAR2(nu_shifted, t, ctx);
        return status;
    }

    GR_TMP_INIT2(z2, z2_minus_nu2, ctx);

    status |= gr_sqr(z2, z, ctx);
    status |= gr_sqr(z2_minus_nu2, nu, ctx);
    status |= gr_sub(z2_minus_nu2, z2, z2_minus_nu2, ctx);

    /* res[2] = -(z * res[1] + (z^2 - nu^2) * res[0]) / (2 * z^2) */
    gr_ptr c2 = GR_ENTRY(res, 2, sz);
    status |= gr_mul(c2, z, GR_ENTRY(res, 1, sz), ctx);
    status |= gr_addmul(c2, z2_minus_nu2, GR_ENTRY(res, 0, sz), ctx);
    status |= gr_div(c2, c2, z2, ctx);
    status |= gr_mul_2exp_si(c2, c2, -1, ctx);
    status |= gr_neg(c2, c2, ctx);

    if (len == 3)
    {
        GR_TMP_CLEAR2(nu_shifted, t, ctx);
        GR_TMP_CLEAR2(z2, z2_minus_nu2, ctx);
        return status;
    }

    fmpz_init(tmp_fmpz);

    /* res[3] = -(6 * z * res[2] + (1 + z^2 - nu^2) * res[1] + 2 * z * res[0]) / (6 * z^2) */
    gr_ptr c3 = GR_ENTRY(res, 3, sz);
    status |= gr_mul(c3, z, GR_ENTRY(res, 2, sz), ctx);
    status |= gr_mul_si(c3, c3, 6, ctx);
    status |= gr_add_si(t, z2_minus_nu2, 1, ctx);
    status |= gr_addmul(c3, t, GR_ENTRY(res, 1, sz), ctx);
    status |= gr_mul_2exp_si(t, z, 1, ctx);
    status |= gr_addmul(c3, t, GR_ENTRY(res, 0, sz), ctx);
    status |= gr_mul_si(t, z2, 6, ctx);
    status |= gr_div(c3, c3, t, ctx);
    status |= gr_neg(c3, c3, ctx);

    for (k = 2; k < len - 2; k++)
    {
        gr_ptr ck2 = GR_ENTRY(res, k+2, sz);
        gr_ptr ck1 = GR_ENTRY(res, k+1, sz);
        gr_ptr ck  = GR_ENTRY(res, k,   sz);
        gr_ptr ckm1 = GR_ENTRY(res, k-1, sz);
        gr_ptr ckm2 = GR_ENTRY(res, k-2, sz);

        /* res[k+2] = (k + 1) * (2k + 1) * z * res[k+1] */
        status |= gr_mul(ck2, z, ck1, ctx);
        fmpz_set_ui(tmp_fmpz, k + 1);
        fmpz_mul_ui(tmp_fmpz, tmp_fmpz, 2 * k + 1);
        status |= gr_mul_fmpz(ck2, ck2, tmp_fmpz, ctx);

        /* res[k+2] += (k^2 + z^2 - nu^2) * res[k] */
        status |= gr_add_si(t, z2_minus_nu2, k * k, ctx);
        status |= gr_addmul(ck2, t, ck, ctx);

        /* res[k+2] += 2 * z * res[k-1] */
        status |= gr_mul_2exp_si(t, z, 1, ctx);
        status |= gr_addmul(ck2, t, ckm1, ctx);

        /* res[k+2] += res[k-2] */
        status |= gr_add(ck2, ck2, ckm2, ctx);

        /* res[k+2] /= -((k+1) * (k+2) * z^2) */
        fmpz_set_ui(tmp_fmpz, k + 1);
        fmpz_mul_ui(tmp_fmpz, tmp_fmpz, k + 2);
        status |= gr_mul_fmpz(t, z2, tmp_fmpz, ctx);
        status |= gr_div(ck2, ck2, t, ctx);
        status |= gr_neg(ck2, ck2, ctx);
    }

    GR_TMP_CLEAR2(nu_shifted, t, ctx);
    GR_TMP_CLEAR2(z2, z2_minus_nu2, ctx);
    fmpz_clear(tmp_fmpz);

    return status;
}
