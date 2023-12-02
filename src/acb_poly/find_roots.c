/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpq.h"
#include "acb_poly.h"

slong
_acb_get_mid_mag(const acb_t z)
{
    slong rm, im;

    rm = arf_abs_bound_lt_2exp_si(arb_midref(acb_realref(z)));
    im = arf_abs_bound_lt_2exp_si(arb_midref(acb_imagref(z)));

    return FLINT_MAX(rm, im);
}

slong
_acb_get_rad_mag(const acb_t z)
{
    slong rm, im;

    /* TODO: write mag function */
    arf_t t;
    arf_init(t);

    arf_set_mag(t, arb_radref(acb_realref(z)));
    rm = arf_abs_bound_lt_2exp_si(t);

    arf_set_mag(t, arb_radref(acb_imagref(z)));
    im = arf_abs_bound_lt_2exp_si(t);

    arf_clear(t);

    return FLINT_MAX(rm, im);
}

void
_acb_poly_roots_initial_values(acb_ptr roots, slong deg, slong prec)
{
    slong i;

    fmpq_t q;
    fmpq_init(q);

    fmpq_set_si(q, 4, 10);
    arb_set_fmpq(acb_realref(roots + 0), q, prec);
    mag_zero(arb_radref(acb_realref(roots + 0)));

    fmpq_set_si(q, 9, 10);
    arb_set_fmpq(acb_imagref(roots + 0), q, prec);
    mag_zero(arb_radref(acb_imagref(roots + 0)));
    fmpq_clear(q);

    for (i = 1; i < deg; i++)
    {
        acb_mul(roots + i, roots + i - 1, roots + 0, prec);
        mag_zero(arb_radref(acb_realref(roots + i)));
        mag_zero(arb_radref(acb_imagref(roots + i)));
    }
}

slong
_acb_poly_find_roots(acb_ptr roots,
    acb_srcptr poly,
    acb_srcptr initial, slong len, slong maxiter, slong prec)
{
    slong iter, i, deg;
    slong rootmag, max_rootmag, correction, max_correction;

    deg = len - 1;

    if (deg == 0)
    {
        return 0;
    }
    else if (acb_contains_zero(poly + len - 1))
    {
        /* if the leading coefficient contains zero, roots can be anywhere */
        for (i = 0; i < deg; i++)
        {
            arb_zero_pm_inf(acb_realref(roots + i));
            arb_zero_pm_inf(acb_imagref(roots + i));
        }
        return 0;
    }
    else if (deg == 1)
    {
        acb_inv(roots + 0, poly + 1, prec);
        acb_mul(roots + 0, roots + 0, poly + 0, prec);
        acb_neg(roots + 0, roots + 0);
        return 1;
    }

    if (initial == NULL)
        _acb_poly_roots_initial_values(roots, deg, prec);
    else
        _acb_vec_set(roots, initial, deg);

    if (maxiter == 0)
        maxiter = 2 * deg + n_sqrt(prec);

    for (iter = 0; iter < maxiter; iter++)
    {
        max_rootmag = -ARF_PREC_EXACT;
        for (i = 0; i < deg; i++)
        {
            rootmag = _acb_get_mid_mag(roots + i);
            max_rootmag = FLINT_MAX(rootmag, max_rootmag);
        }

        _acb_poly_refine_roots_durand_kerner(roots, poly, len, prec);

        max_correction = -ARF_PREC_EXACT;
        for (i = 0; i < deg; i++)
        {
            correction = _acb_get_rad_mag(roots + i);
            max_correction = FLINT_MAX(correction, max_correction);
        }

        /* estimate the correction relative to the whole set of roots */
        max_correction -= max_rootmag;

        /* flint_printf("ITER %wd MAX CORRECTION: %wd\n", iter, max_correction); */

        if (max_correction < -prec / 2)
            maxiter = FLINT_MIN(maxiter, iter + 2);
        else if (max_correction < -prec / 3)
            maxiter = FLINT_MIN(maxiter, iter + 3);
        else if (max_correction < -prec / 4)
            maxiter = FLINT_MIN(maxiter, iter + 4);
    }

    return _acb_poly_validate_roots(roots, poly, len, prec);
}


slong
acb_poly_find_roots(acb_ptr roots,
    const acb_poly_t poly, acb_srcptr initial,
    slong maxiter, slong prec)
{
    slong len = poly->length;

    if (len == 0)
    {
        flint_throw(FLINT_ERROR, "find_roots: expected a nonzero polynomial");
    }

    return _acb_poly_find_roots(roots, poly->coeffs, initial,
                len, maxiter, prec);
}

