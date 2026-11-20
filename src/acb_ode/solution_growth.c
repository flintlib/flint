/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_types.h"
#include "acb.h"
#include "acb_ode.h"
#include "acb_poly.h"

static void
l1_maj(acb_poly_t eqn)
{
    acb_t tmp;
    acb_init(tmp);

    for (slong i = 0; i < eqn->length; i++)
    {
        acb_swap(tmp, eqn->coeffs + i);
        /* l1 norm because of the wrapping effect */
        arb_abs(acb_realref(tmp), acb_realref(tmp));
        arb_abs(acb_imagref(tmp), acb_imagref(tmp));
        arb_add(acb_realref(eqn->coeffs + i),
                acb_realref(tmp), acb_imagref(tmp),
                MAG_BITS);
        arb_zero(acb_imagref(eqn->coeffs + i));
    }

    acb_ptr lc = eqn->coeffs + eqn->length - 1;
    FLINT_ASSERT(acb_is_one(lc));
    acb_set_si(lc, -1);

    acb_clear(tmp);
}

static void
max_root(mag_t base, const acb_poly_t eqn, slong p, slong q)
{
    acb_ptr rts = _acb_vec_init(eqn->length - 1);
    mag_t mag;
    mag_init(mag);

    /* change to a cheap but rigorous lower bound? (not needed at the moment but
       may be useful for other applications) */
    double eps = _acb_poly_find_roots_double(rts, eqn->coeffs, NULL,
                                             eqn->length, 32, 16);
    if (eps > 0x1p-8)
        _acb_poly_find_roots(rts, eqn->coeffs, rts, eqn->length, 0, 16);

    mag_zero(base);
    for (slong i = 0; i < eqn->length - 1; i++)
    {
        acb_get_mag(mag, rts + i);
        mag_max(base, base, mag);
    }

    _acb_vec_clear(rts, eqn->length - 1);
    mag_clear(mag);
}

/* max growth of
   - coefficient sequence ≈ base^n / n!^(1/order)
   - function at infinity (if entire) ≈ exp((1/order)*(base*x)^order) */

void
_acb_ode_solution_growth(mag_t order, mag_t base, mag_t maj_base,
                         const acb_poly_struct * dop, slong dop_len)
{
    slong p = -1, q = 0;

    FLINT_ASSERT (!acb_is_zero(dop[dop_len - 1].coeffs + 0));

    /* (j0, i0) = origin of the main edge of the Newton polygon under the
       above assumption */
    slong j0 = 0;
    slong i0 = dop_len - 1;

    slong lcd = acb_poly_length(dop + dop_len - 1) - 1;

    if (lcd >= 1)  /* expected radius of convergence is finite */
    {
        p = 0;
        q = lcd;
    }
    else  /* infinite radius of convergence */
    {
        /* Newton polygon:
           If the general term of f(x) is ≈ base^n/n!^(1/order)·x^n,
           applying x^j·(x·d/dx)^i to f(x) multiplies the general term by
           ≈ base^(-j)·n^(i+j/order). So the possible values for order are those
           such that i+j/order is reached twice or more, and we want the largest
           of these. */

        for (slong i = 0; i < dop_len; i++)
            for (slong j = 0; j < dop[i].length; j++)
            {
                if (j > j0 && i < i0 && !acb_is_zero(dop[i].coeffs + j))
                {
                    /* (j, i) belongs to the Newton polygon;
                       slope = Δi/Δj (sic) = p/q */
                    slong p1 = i - i0, q1 = j - j0;
                    if (p1 * q > p * q1)
                        p = p1, q = q1;
                }
            }
    }

    FLINT_ASSERT(p <= 0);
    FLINT_ASSERT(q >= 0);

    mag_set_ui(order, q);
    if (q == 0)  /* poly(theta), rec solns have finite support */
    {
        mag_zero(base);
        mag_inf(order);
        return;
    }
    mag_div_ui(order, order, -p);  /* order = -1/slope */

    /* estimate the roots of the characteristic equations of the recurrence and
       the majorant recurrence */

    acb_poly_t eqn;
    acb_poly_init(eqn);

    for (slong i = i0 + p; i <= i0; i++)
        for (slong j = j0; j <= FLINT_MIN(j0 + q, dop[i].length - 1); j++)
            if ((i - i0) * q == p * (j - j0))
                acb_poly_set_coeff_acb(eqn, j0 + q - j, dop[i].coeffs + j);

    acb_ptr lc = eqn->coeffs + eqn->length - 1;
    acb_poly_scalar_div(eqn, eqn, lc, MAG_BITS);
    acb_one(lc);

    max_root(base, eqn, p, q);

    l1_maj(eqn);
    max_root(maj_base, eqn, p, q);

    acb_poly_clear(eqn);
}
