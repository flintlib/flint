/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_fmpz_poly.h"
#include "ca_qqbar.h"

int
_ca_qqbar_validate_enclosure(acb_t res, const fmpz_poly_t poly, const acb_t z, slong prec)
{
    if (!acb_is_finite(z) || fmpz_poly_degree(poly) < 1)
    {
        if (res != NULL)
            acb_set(res, z);
        return 0;
    }
    else if (acb_is_exact(z) || fmpz_poly_degree(poly) == 1)
    {
        if (res != NULL)
            acb_set(res, z);
        return 1;
    }
    else
    {
        slong acc;
        acb_t z2, zmid, t, u;
        fmpz_poly_t deriv;
        mag_t zmag, eps;
        int pure_real, ans;

        pure_real = acb_is_real(z);

        if (prec == 0)
        {
            acc = acb_rel_accuracy_bits(z);
            /* todo: maybe use the bits * degree of poly instead as bound */
            acc = FLINT_MIN(acc, ARF_PREC_EXACT / 4);
            prec = FLINT_MAX(acc, 32) * 2 + 10;
        }

        acb_init(z2);
        acb_init(zmid);
        acb_init(t);
        acb_init(u);
        mag_init(eps);
        mag_init(zmag);
        fmpz_poly_init(deriv);

        /* Slightly inflate enclosure -- needed e.g. in case of complex
           interval with one very narrow real/imaginary part */
        acb_get_mag(zmag, z);
#if 0
        mag_mul_2exp_si(zmag, zmag, -prec / 2);
#else
        mag_mul_2exp_si(zmag, zmag, -3 * prec / 4);
#endif

        mag_hypot(eps, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
        mag_mul_2exp_si(eps, eps, -4);

        mag_max(eps, eps, zmag);

        acb_set(z2, z);
        if (pure_real)
            arb_add_error_mag(acb_realref(z2), eps);
        else
            acb_add_error_mag(z2, eps);

        acb_get_mid(zmid, z2);
        fmpz_poly_derivative(deriv, poly);

        /* Do one interval Newton step, t = zmid - poly(zmid) / poly'(z2) */
        arb_fmpz_poly_evaluate_acb(t, poly, zmid, prec);
        arb_fmpz_poly_evaluate_acb(u, deriv, z2, prec);
        acb_div(t, t, u, prec);
        acb_sub(t, zmid, t, prec);

        if (pure_real)
        {
            ans = arb_contains_interior(acb_realref(z2), acb_realref(t));
            arb_zero(acb_imagref(t));
        }
        else
        {
            ans = acb_contains_interior(z2, t);
        }

        if (res != NULL)
        {
            if (ans)
                acb_set(res, t);
            else
                acb_set(res, z);
        }

        acb_clear(z2);
        acb_clear(zmid);
        acb_clear(t);
        acb_clear(u);
        mag_clear(eps);
        mag_clear(zmag);
        fmpz_poly_clear(deriv);

        return ans;
    }
}

