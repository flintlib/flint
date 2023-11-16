/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

#ifdef __GNUC__
# define log __builtin_log
#else
# include <math.h>
#endif

void
acb_hypgeom_erfc(acb_t res, const acb_t z, slong prec)
{
    double x, y, abs_z2, log_z, log_erfc_z_asymp;

    if (!acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_is_real(z))
    {
        arb_hypgeom_erfc(acb_realref(res), acb_realref(z), prec);
        arb_zero(acb_imagref(res));
        return;
    }

    if (acb_is_zero(z))
    {
        acb_one(res);
        return;
    }

    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), -64) < 0 &&
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), -64) < 0) ||
        arf_sgn(arb_midref(acb_realref(z))) < 0)
    {
        acb_hypgeom_erf(res, z, prec);
        acb_sub_ui(res, res, 1, prec);
        acb_neg(res, res);
        return;
    }

    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 64) > 0 ||
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 64) > 0))
    {
        acb_hypgeom_erf_asymp(res, z, 1, prec, prec);
        return;
    }

    x = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
    y = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);

    abs_z2 = x * x + y * y;

    if (abs_z2 > (prec + 8) * 0.69314718055994530942)
    {
        acb_hypgeom_erf_asymp(res, z, 1, prec, prec);
    }
    else
    {
        slong wp;

        log_z = 0.5 * log(abs_z2);
        log_erfc_z_asymp = y * y - x * x - log_z;

        wp = prec + 2;
        if (log_erfc_z_asymp < 0.0)
            wp += (-log_erfc_z_asymp) * 1.4426950408889634074;

        if (acb_rel_accuracy_bits(z) >= wp)
        {
            acb_hypgeom_erf(res, z, wp);
        }
        else
        {
            acb_t zmid;
            mag_t re_err, im_err;

            acb_init(zmid);
            mag_init(re_err);
            mag_init(im_err);

            acb_hypgeom_erf_propagated_error(re_err, im_err, z);
            arf_set(arb_midref(acb_realref(zmid)), arb_midref(acb_realref(z)));
            arf_set(arb_midref(acb_imagref(zmid)), arb_midref(acb_imagref(z)));

            acb_hypgeom_erf(res, zmid, wp);

            arb_add_error_mag(acb_realref(res), re_err);
            arb_add_error_mag(acb_imagref(res), im_err);

            acb_clear(zmid);
            mag_clear(re_err);
            mag_clear(im_err);
        }

        acb_sub_ui(res, res, 1, prec);
        acb_neg(res, res);
    }
}

