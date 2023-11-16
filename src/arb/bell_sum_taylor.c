/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_poly.h"

/* tuning parameter */
#define RADIUS_BITS 3

void
_arb_bell_sum_taylor(arb_t res, const fmpz_t n,
        const fmpz_t a, const fmpz_t b, const fmpz_t mmag, slong tol)
{
    fmpz_t m, r, R, tmp;
    mag_t B, C, D, bound;
    arb_t t, u;
    slong wp, k, N;

    if (_fmpz_sub_small(b, a) < 5)
    {
        arb_bell_sum_bsplit(res, n, a, b, mmag, tol);
        return;
    }

    fmpz_init(m);
    fmpz_init(r);
    fmpz_init(R);
    fmpz_init(tmp);

    /* r = max(m - a, b - m) */
    /* m = a + (b - a) / 2 */
    fmpz_sub(r, b, a);
    fmpz_cdiv_q_2exp(r, r, 1);
    fmpz_add(m, a, r);

    fmpz_mul_2exp(R, r, RADIUS_BITS);

    mag_init(B);
    mag_init(C);
    mag_init(D);
    mag_init(bound);

    arb_init(t);
    arb_init(u);

    if (fmpz_cmp(R, m) >= 0)
    {
        mag_inf(C);
        mag_inf(D);
    }
    else
    {
        /* C = exp(R * |F'(m)| + (1/2) R^2 * (n/(m-R)^2 + 1/(m-R))) */
        /* C = exp(R * (|F'(m)| + (1/2) R * (n/(m-R) + 1)/(m-R))) */
        /* D = (1/2) R * (n/(m-R) + 1)/(m-R) */
        fmpz_sub(tmp, m, R);
        mag_set_fmpz(D, n);
        mag_div_fmpz(D, D, tmp);
        mag_one(C);
        mag_add(D, D, C);
        mag_div_fmpz(D, D, tmp);
        mag_mul_fmpz(D, D, R);
        mag_mul_2exp_si(D, D, -1);

        /* C = |F'(m)| */
        wp = 20 + 1.05 * fmpz_bits(n);
        arb_set_fmpz(t, n);
        arb_div_fmpz(t, t, m, wp);
        fmpz_add_ui(tmp, m, 1);
        arb_set_fmpz(u, tmp);
        arb_digamma(u, u, wp);
        arb_sub(t, t, u, wp);
        arb_get_mag(C, t);

        /* C = exp(R * (C + D)) */
        mag_add(C, C, D);
        mag_mul_fmpz(C, C, R);
        mag_exp(C, C);
    }

    if (mag_cmp_2exp_si(C, tol / 4 + 2) > 0)
    {
        _arb_bell_sum_taylor(res, n, a, m, mmag, tol);
        _arb_bell_sum_taylor(t, n, m, b, mmag, tol);
        arb_add(res, res, t, 2 * tol);
    }
    else
    {
        arb_ptr mx, ser1, ser2, ser3;

        /* D = T(m) */
        wp = 20 + 1.05 * fmpz_bits(n);
        arb_set_fmpz(t, m);
        arb_pow_fmpz(t, t, n, wp);
        fmpz_add_ui(tmp, m, 1);
        arb_gamma_fmpz(u, tmp, wp);
        arb_div(t, t, u, wp);
        arb_get_mag(D, t);

        /* error bound: (b-a) * C * D * B^N / (1 - B), B = r/R */
        /*              ((b-a) * C * D * 2) * 2^(-N*RADIUS_BITS) */

        /* ((b-a) * C * D * 2) */
        mag_mul(bound, C, D);
        mag_mul_2exp_si(bound, bound, 1);
        fmpz_sub(tmp, b, a);
        mag_mul_fmpz(bound, bound, tmp);

        /* N = (tol + log2((b-a)*C*D*2) - mmag) / RADIUS_BITS */
        if (mmag == NULL)
        {
            /* estimate D ~= 2^mmag */
            fmpz_add_ui(tmp, MAG_EXPREF(C), tol);
            fmpz_cdiv_q_ui(tmp, tmp, RADIUS_BITS);
        }
        else
        {
            fmpz_sub(tmp, MAG_EXPREF(bound), mmag);
            fmpz_add_ui(tmp, tmp, tol);
            fmpz_cdiv_q_ui(tmp, tmp, RADIUS_BITS);
        }

        if (fmpz_cmp_ui(tmp, 5 * tol / 4) > 0)
            N = 5 * tol / 4;
        else if (fmpz_cmp_ui(tmp, 2) < 0)
            N = 2;
        else
            N = fmpz_get_ui(tmp);

        /* multiply by 2^(-N*RADIUS_BITS) */
        mag_mul_2exp_si(bound, bound, -N * RADIUS_BITS);

        mx = _arb_vec_init(2);
        ser1 = _arb_vec_init(N);
        ser2 = _arb_vec_init(N);
        ser3 = _arb_vec_init(N);

        /* estimate (this should work for moderate n and tol) */
        wp = 1.1 * tol + 1.05 * fmpz_bits(n) + 5;

        /* increase precision until convergence */
        while (1)
        {
            /* (m+x)^n / gamma(m+1+x) */
            arb_set_fmpz(mx, m);
            arb_one(mx + 1);
            _arb_poly_log_series(ser1, mx, 2, N, wp);
            for (k = 0; k < N; k++)
                arb_mul_fmpz(ser1 + k, ser1 + k, n, wp);
            arb_add_ui(mx, mx, 1, wp);
            _arb_poly_lgamma_series(ser2, mx, 2, N, wp);
            _arb_vec_sub(ser1, ser1, ser2, N, wp);
            _arb_poly_exp_series(ser3, ser1, N, N, wp);

            /* t = a - m, u = b - m */
            arb_set_fmpz(t, a);
            arb_sub_fmpz(t, t, m, wp);
            arb_set_fmpz(u, b);
            arb_sub_fmpz(u, u, m, wp);
            arb_power_sum_vec(ser1, t, u, N, wp);

            arb_zero(res);
            for (k = 0; k < N; k++)
                arb_addmul(res, ser3 + k, ser1 + k, wp);

            if (mmag != NULL)
            {
                if (_fmpz_sub_small(MAG_EXPREF(arb_radref(res)), mmag) <= -tol)
                    break;
            }
            else
            {
                if (arb_rel_accuracy_bits(res) >= tol)
                    break;
            }

            wp = 2 * wp;
        }

        /* add the series truncation bound */
        arb_add_error_mag(res, bound);

        _arb_vec_clear(mx, 2);
        _arb_vec_clear(ser1, N);
        _arb_vec_clear(ser2, N);
        _arb_vec_clear(ser3, N);
    }

    mag_clear(B);
    mag_clear(C);
    mag_clear(D);
    mag_clear(bound);
    arb_clear(t);
    arb_clear(u);

    fmpz_clear(m);
    fmpz_clear(r);
    fmpz_clear(R);
    fmpz_clear(tmp);
}

void
arb_bell_sum_taylor(arb_t res, const fmpz_t n,
        const fmpz_t a, const fmpz_t b, const fmpz_t mmag, slong prec)
{
    _arb_bell_sum_taylor(res, n, a, b, mmag, prec + 5);
}

