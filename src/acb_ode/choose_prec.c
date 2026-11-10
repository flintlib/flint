/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>

#include "acb_types.h"
#include "acb_poly.h"
#include "acb_ode.h"
#include "ulong_extras.h"

/* XXX
   - Include overestimation of tail bound? how?
   - Cap to input accuracy?
   - Better overestimation estimate based on the roots of the characteristic
     equation (cf. _acb_ode_solution_growth). Idea (implicit in a debugging
     session with claude): if the dominant root (or root cluster) α has
     multiplicity m, the coefficient sequence grows roughly like the solutions
     of (1 - α S^{-1})^m, while the radii grow roughly like the solutions of
     (1 + |α| S^{-1})^m - 2 (absolute values everywhere, flip sign of lc). The
     magnitude ratio of the dominant roots is 1/(2^{1/m}-1) <= m/ln(2), so we
     expect to lose about log2(m/ln(2)) bits per term. */
slong
acb_ode_choose_prec(slong * rec_prec, const acb_poly_struct * dop, slong dop_len,
                    mag_srcptr rad, mag_srcptr cvrad, slong tgt_prec)
{
    slong sum_prec, prec0, prec;
    double lgmag, nterms;
    mag_t ratio;

    mag_init(ratio);

    prec0 = tgt_prec + 2*dop_len;

    /* rough estimate of bits lost to overestimation when unrolling a recurrence
       in interval arithmetic
       todo: linear error analysis */
    slong dop_clen = _acb_poly_vec_length(dop, dop_len);
    prec = 8 + prec0 * (1 + n_clog(FLINT_MAX(1, dop_clen - 3), 2));

    lgmag = 0.;
    if (mag_is_zero(rad))
    {
        nterms = 1; /* todo: better choice? */
    }
    else if (mag_is_finite(cvrad))
    {
        mag_div_lower(ratio, cvrad, rad);
        nterms = prec / mag_get_d_log2_approx(ratio);
    }
    else
    {
        /* estimate hump height based on worst-case asymptotic behavior */
        mag_t order, base;

        mag_init(order);
        mag_init(base);

        _acb_ode_solution_growth(order, base, dop, dop_len);

        /* mag so we don't have to worry about overflows */
        mag_mul(base, base, rad);

        double base_d = mag_get_d(base);
        double order_d = mag_get_d(order);

        double hump = exp(log(base_d) * order_d + 1.);
        /* cap the cancellation we'll attempt to absorb */
        hump = FLINT_MIN(hump, 100. + prec * log2(prec));
        double den = log2(prec) / order_d;
        if (base_d < 1.)
            den = FLINT_MAX(den, - log2(base_d));
        nterms = hump + prec / den;
        lgmag = FLINT_MAX(0., order_d * hump * log2(hump));

        mag_clear(order);
        mag_clear(base);
    }

    double prec1 = 2 * lgmag + log2(nterms);
    if (!(fabs(prec1) <= prec0))
        prec1 = prec0;

    sum_prec = 8 + 1.125 * (prec0 + prec1);
    prec = FLINT_MAX(prec, sum_prec);

    mag_clear(ratio);

    if (rec_prec != NULL)
        * rec_prec = prec;
    return sum_prec;
}
