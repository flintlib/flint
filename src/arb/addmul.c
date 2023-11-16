/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_addmul_arf(arb_t z, const arb_t x, const arf_t y, slong prec)
{
    mag_t ym;
    int inexact;

    if (arb_is_exact(x))
    {
        inexact = arf_addmul(arb_midref(z), arb_midref(x), y, prec, ARB_RND);

        if (inexact)
            arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
    }
    else if (ARB_IS_LAGOM(x) && ARF_IS_LAGOM(y) && ARB_IS_LAGOM(z))
    {
        mag_fast_init_set_arf(ym, y);
        mag_fast_addmul(arb_radref(z), ym, arb_radref(x));
        inexact = arf_addmul(arb_midref(z), arb_midref(x), y, prec, ARB_RND);

        if (inexact)
            arf_mag_fast_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
    }
    else if (arf_is_inf(y) && arb_is_nonzero(x))
    {
        if (arf_sgn(arb_midref(x)) > 0)
            arb_add_arf(z, z, y, prec);
        else
            arb_sub_arf(z, z, y, prec);
    }
    else
    {
        mag_init_set_arf(ym, y);
        mag_addmul(arb_radref(z), ym, arb_radref(x));

        inexact = arf_addmul(arb_midref(z), arb_midref(x), y, prec, ARB_RND);
        if (inexact)
            arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);

        mag_clear(ym);
    }
}

void
arb_addmul(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    mag_t zr, xm, ym;
    int inexact;

    if (arb_is_exact(y))
    {
        arb_addmul_arf(z, x, arb_midref(y), prec);
    }
    else if (arb_is_exact(x))
    {
        arb_addmul_arf(z, y, arb_midref(x), prec);
    }
    else if (ARB_IS_LAGOM(x) && ARB_IS_LAGOM(y) && ARB_IS_LAGOM(z))
    {
        mag_fast_init_set_arf(xm, arb_midref(x));
        mag_fast_init_set_arf(ym, arb_midref(y));

        mag_fast_init_set(zr, arb_radref(z));
        mag_fast_addmul(zr, xm, arb_radref(y));
        mag_fast_addmul(zr, ym, arb_radref(x));
        mag_fast_addmul(zr, arb_radref(x), arb_radref(y));

        inexact = arf_addmul(arb_midref(z), arb_midref(x), arb_midref(y),
            prec, ARF_RND_DOWN);

        if (inexact)
            arf_mag_fast_add_ulp(zr, zr, arb_midref(z), prec);

        *arb_radref(z) = *zr;
    }
    else if (arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)) && arb_is_nonzero(y))
    {
        if (arf_sgn(arb_midref(y)) > 0)
            arb_add_arf(z, z, arb_midref(x), prec);
        else
            arb_sub_arf(z, z, arb_midref(x), prec);
    }
    else if (arf_is_inf(arb_midref(y)) && mag_is_finite(arb_radref(y)) && arb_is_nonzero(x))
    {
        if (arf_sgn(arb_midref(x)) > 0)
            arb_add_arf(z, z, arb_midref(y), prec);
        else
            arb_sub_arf(z, z, arb_midref(y), prec);
    }
    else
    {
        mag_init_set_arf(xm, arb_midref(x));
        mag_init_set_arf(ym, arb_midref(y));

        mag_init_set(zr, arb_radref(z));
        mag_addmul(zr, xm, arb_radref(y));
        mag_addmul(zr, ym, arb_radref(x));
        mag_addmul(zr, arb_radref(x), arb_radref(y));

        inexact = arf_addmul(arb_midref(z), arb_midref(x), arb_midref(y),
            prec, ARF_RND_DOWN);

        if (inexact)
            arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
        else
            mag_set(arb_radref(z), zr);

        mag_clear(zr);
        mag_clear(xm);
        mag_clear(ym);
    }
}

void
arb_addmul_si(arb_t z, const arb_t x, slong y, slong prec)
{
    arf_t t;
    arf_init_set_si(t, y); /* no need to free */
    arb_addmul_arf(z, x, t, prec);
}

void
arb_addmul_ui(arb_t z, const arb_t x, ulong y, slong prec)
{
    arf_t t;
    arf_init_set_ui(t, y); /* no need to free */
    arb_addmul_arf(z, x, t, prec);
}

void
arb_addmul_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)
{
    arf_t t;

    if (!COEFF_IS_MPZ(*y))
    {
        arf_init_set_si(t, *y); /* no need to free */
        arb_addmul_arf(z, x, t, prec);
    }
    else
    {
        arf_init(t);
        arf_set_fmpz(t, y);
        arb_addmul_arf(z, x, t, prec);
        arf_clear(t);
    }
}
