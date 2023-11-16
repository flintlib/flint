/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"
#include "bernoulli.h"

/* todo: use log(1-z) when this is better? would also need to
   adjust strategy in the main function */
void
acb_hypgeom_dilog_bernoulli(acb_t res, const acb_t z, slong prec)
{
    acb_t s, w, w2;
    slong n, k;
    fmpz_t c, d;
    mag_t m, err;
    double lm;
    int real;

    acb_init(s);
    acb_init(w);
    acb_init(w2);
    fmpz_init(c);
    fmpz_init(d);
    mag_init(m);
    mag_init(err);

    real = 0;
    if (acb_is_real(z))
    {
        arb_sub_ui(acb_realref(w), acb_realref(z), 1, 30);
        real = arb_is_nonpositive(acb_realref(w));
    }

    acb_log(w, z, prec);
    acb_get_mag(m, w);

    /* for k >= 4, the terms are bounded by  (|w| / (2 pi))^k */
    mag_set_ui_2exp_si(err, 2670177, -24);  /* upper bound for 1/(2pi) */
    mag_mul(err, err, m);
    lm = mag_get_d_log2_approx(err);

    if (lm < -0.25)
    {
        n = prec / (-lm) + 1;
        n = FLINT_MAX(n, 4);
        mag_geom_series(err, err, n);

        BERNOULLI_ENSURE_CACHED(n)

        acb_mul(w2, w, w, prec);

        for (k = n - (n % 2 == 0); k >= 3; k -= 2)
        {
            fmpz_mul_ui(c, fmpq_denref(bernoulli_cache + k - 1), k - 1);
            fmpz_mul_ui(d, c, (k + 1) * (k + 2));
            acb_mul(s, s, w2, prec);
            acb_mul_fmpz(s, s, c, prec);
            fmpz_mul_ui(c, fmpq_numref(bernoulli_cache + k - 1), (k + 1) * (k + 2));
            acb_sub_fmpz(s, s, c, prec);
            acb_div_fmpz(s, s, d, prec);
        }

        acb_mul(s, s, w, prec);
        acb_mul_2exp_si(s, s, 1);
        acb_sub_ui(s, s, 3, prec);
        acb_mul(s, s, w2, prec);
        acb_mul_2exp_si(s, s, -1);
        acb_const_pi(w2, prec);
        acb_addmul(s, w2, w2, prec);
        acb_div_ui(s, s, 6, prec);

        acb_neg(w2, w);
        acb_log(w2, w2, prec);
        acb_submul(s, w2, w, prec);
        acb_add(res, s, w, prec);

        acb_add_error_mag(res, err);
        if (real)
            arb_zero(acb_imagref(res));
    }
    else
    {
        acb_indeterminate(res);
    }

    acb_clear(s);
    acb_clear(w);
    acb_clear(w2);
    fmpz_clear(c);
    fmpz_clear(d);
    mag_clear(m);
    mag_clear(err);
}

