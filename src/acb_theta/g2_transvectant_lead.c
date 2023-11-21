/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_theta.h"

void
acb_theta_g2_transvectant_lead(acb_t res, const acb_poly_t g, const acb_poly_t h,
    slong m, slong n, slong k, slong prec)
{
    acb_ptr s, t;
    fmpz_t num, f;
    slong j;

    s = _acb_vec_init(k + 1);
    t = _acb_vec_init(k + 1);
    fmpz_init(num);
    fmpz_init(f);

    /* Set i = m - k (resp. n - k) in g2_transvectant and use acb_dot */
    for (j = 0; j <= k; j++)
    {
        acb_poly_get_coeff_acb(&s[j], g, m - j);
        acb_poly_get_coeff_acb(&t[j], h, n - k + j);
        /* Put all factorials in s */
        fmpz_fac_ui(num, m - j);
        fmpz_fac_ui(f, n - k + j);
        fmpz_mul(num, num, f);
        if ((k - j) % 2 == 1)
        {
            fmpz_neg(num, num);
        }
        acb_mul_fmpz(&s[j], &s[j], num, prec);
    }
    acb_dot(res, NULL, 0, s, 1, t, 1, k + 1, prec);

    fmpz_fac_ui(num, k);
    acb_set_fmpz(t, num);
    fmpz_fac_ui(num, m);
    fmpz_fac_ui(f, n);
    fmpz_mul(num, num, f);
    acb_div_fmpz(t, t, num, prec);
    acb_mul(res, res, t, prec);

    _acb_vec_clear(s, k + 1);
    _acb_vec_clear(t, k + 1);
    fmpz_clear(num);
    fmpz_clear(f);
}
