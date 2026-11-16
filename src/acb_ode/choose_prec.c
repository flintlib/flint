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
#include "acb_ode.h"
#include "mag.h"

#define LOG2 0.69314718055994531

/* XXX
   - Include overestimation of tail bound? how?
   - Cap to input accuracy? */
slong
acb_ode_choose_prec(slong * rec_prec,
                    const acb_poly_struct * dop, slong dop_len,
                    mag_srcptr rad, mag_srcptr cvrad, slong tgt_prec)
{
    slong sum_prec, prec0, wp;
    double lgmag, nterms;
    mag_t order, base, maj_base, maj_ratio;

    prec0 = tgt_prec + 2*dop_len;

    if (mag_is_zero(rad))
    {
        /* this assumes that we only care about the sum */
        if (rec_prec != NULL)
            * rec_prec = prec0;
        return prec0;
    }

    mag_init(order);
    mag_init(base);
    mag_init(maj_base);
    mag_init(maj_ratio);

    _acb_ode_solution_growth(order, base, maj_base, dop, dop_len);
    if (!mag_is_finite(order))  /* cvrad may be more accurate or > generic */
    {
        if (mag_is_finite(cvrad))
            mag_inv(base, cvrad);
        else
            mag_div_ui(base, base, 10);  /* kludge */
    }
    if (mag_cmp(maj_base, base) > 0)
        mag_div(maj_ratio, maj_base, base);
    else
        /* this should only happen due to inaccurate root-finding */
        mag_one(maj_ratio);

    mag_mul(base, base, rad);
    double log2_base_d = mag_get_d_log2_approx(base);

    if (!mag_is_finite(order))
    {
        lgmag = 0.;
        nterms = prec0 / (- log2_base_d);
    }
    else
    {
        double order_d = mag_get_d(order);

        /* estimate hump height based on worst-case asymptotic behavior */
        double hump = exp(LOG2 * log2_base_d * order_d + 1.);
        /* cap the cancellation that we will attempt to absorb */
        hump = FLINT_MIN(hump, 100. + prec0 * log2(prec0));
        double den = log2(prec0) / order_d;
        if (log2_base_d < 0.)
            den = FLINT_MAX(den, - log2_base_d);
        nterms = hump + prec0 / den;
        lgmag = FLINT_MAX(0., order_d * hump * log2(hump));
    }

    /* accumulation precision */

    double prec1 = 2 * lgmag + log2(nterms);
    if (prec1 < 0.)
        prec1 = 0.;
    else if (!(prec1 <= prec0))
        prec1 = prec0;

    sum_prec = 8 + 1.125 * (prec0 + prec1);
    sum_prec = FLINT_MAX(FLINT_MIN(sum_prec, 1000 + 10 * tgt_prec), tgt_prec);

    /* If the coefficients are computed iteratively in ball arithmetic, their
       radii grow roughly like a generic solution of the majorant recurrence, so
       we expect to lose about log2(maj_ratio) bits per term.

       (If the characteristic roots are well separated,
       maj_ratio ≈ 1/(2^{1/m}-1) <= m/ln(2)
       where m is the total multiplicity of characteristic roots of minimal
       modulus, i.e., when the leading coefficient is non-constant, the total
       multiplicity of minimal singularities.)

       todo: Linear error analysis. */

    wp = sum_prec + nterms * mag_get_d_log2_approx(maj_ratio);
    wp = FLINT_MAX(FLINT_MIN(wp, 10000 + 20 * tgt_prec), sum_prec);

    mag_clear(order);
    mag_clear(base);
    mag_clear(maj_base);
    mag_clear(maj_ratio);

    if (rec_prec != NULL)
        * rec_prec = wp;
    return sum_prec;
}
