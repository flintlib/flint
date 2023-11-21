/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_fma_arf(arb_t res, const arb_t x, const arf_t y, const arb_t z, slong prec)
{
    mag_t ym;
    int inexact;

    if (arb_is_exact(x))
    {
        inexact = arf_fma(arb_midref(res), arb_midref(x), y, arb_midref(z), prec, ARB_RND);

        if (inexact)
            arf_mag_add_ulp(arb_radref(res), arb_radref(z), arb_midref(res), prec);
        else
            mag_set(arb_radref(res), arb_radref(z));
    }
    else if (ARB_IS_LAGOM(res) && ARB_IS_LAGOM(x) && ARF_IS_LAGOM(y) && ARB_IS_LAGOM(z))
    {
        mag_t tm;

        mag_fast_init_set_arf(ym, y);
        *tm = *arb_radref(z);
        mag_fast_addmul(tm, ym, arb_radref(x));
        *arb_radref(res) = *tm;

        inexact = arf_fma(arb_midref(res), arb_midref(x), y, arb_midref(z), prec, ARB_RND);
        if (inexact)
            arf_mag_fast_add_ulp(arb_radref(res), arb_radref(res), arb_midref(res), prec);
    }
    else if (arf_is_inf(y) && arb_is_nonzero(x))
    {
        if (arf_sgn(arb_midref(x)) > 0)
            arb_add_arf(res, z, y, prec);
        else
            arb_sub_arf(res, z, y, prec);
    }
    else
    {
        mag_t tm;
        mag_init(tm);

        mag_init_set_arf(ym, y);
        mag_set(tm, arb_radref(z));

        mag_addmul(tm, ym, arb_radref(x));
        mag_set(arb_radref(res), tm);

        inexact = arf_fma(arb_midref(res), arb_midref(x), y, arb_midref(z), prec, ARB_RND);
        if (inexact)
            arf_mag_add_ulp(arb_radref(res), arb_radref(res), arb_midref(res), prec);

        mag_clear(tm);
        mag_clear(ym);
    }
}

void
arb_fma(arb_t res, const arb_t x, const arb_t y, const arb_t z, slong prec)
{
    mag_t zr, xm, ym;
    int inexact;

    if (arb_is_exact(y))
    {
        arb_fma_arf(res, x, arb_midref(y), z, prec);
    }
    else if (arb_is_exact(x))
    {
        arb_fma_arf(res, y, arb_midref(x), z, prec);
    }
    else if (ARB_IS_LAGOM(res) && ARB_IS_LAGOM(x) && ARB_IS_LAGOM(y) && ARB_IS_LAGOM(z))
    {
        mag_fast_init_set_arf(xm, arb_midref(x));
        mag_fast_init_set_arf(ym, arb_midref(y));

        mag_fast_init_set(zr, arb_radref(z));
        mag_fast_addmul(zr, xm, arb_radref(y));
        mag_fast_addmul(zr, ym, arb_radref(x));
        mag_fast_addmul(zr, arb_radref(x), arb_radref(y));

        inexact = arf_fma(arb_midref(res), arb_midref(x), arb_midref(y), arb_midref(z),
            prec, ARF_RND_DOWN);

        if (inexact)
            arf_mag_fast_add_ulp(zr, zr, arb_midref(res), prec);

        *arb_radref(res) = *zr;
    }
    else if (arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)) && arb_is_nonzero(y))
    {
        if (arf_sgn(arb_midref(y)) > 0)
            arb_add_arf(res, z, arb_midref(x), prec);
        else
            arb_sub_arf(res, z, arb_midref(x), prec);
    }
    else if (arf_is_inf(arb_midref(y)) && mag_is_finite(arb_radref(y)) && arb_is_nonzero(x))
    {
        if (arf_sgn(arb_midref(x)) > 0)
            arb_add_arf(res, z, arb_midref(y), prec);
        else
            arb_sub_arf(res, z, arb_midref(y), prec);
    }
    else
    {
        mag_init_set_arf(xm, arb_midref(x));
        mag_init_set_arf(ym, arb_midref(y));

        mag_init_set(zr, arb_radref(z));
        mag_addmul(zr, xm, arb_radref(y));
        mag_addmul(zr, ym, arb_radref(x));
        mag_addmul(zr, arb_radref(x), arb_radref(y));

        inexact = arf_fma(arb_midref(res), arb_midref(x), arb_midref(y), arb_midref(z),
            prec, ARF_RND_DOWN);

        if (inexact)
            arf_mag_add_ulp(arb_radref(res), zr, arb_midref(res), prec);
        else
            mag_set(arb_radref(res), zr);

        mag_clear(zr);
        mag_clear(xm);
        mag_clear(ym);
    }
}

void
arb_fma_ui(arb_t res, const arb_t x, ulong y, const arb_t z, slong prec)
{
    arf_t t;
    arf_init_set_ui(t, y); /* no need to free */
    arb_fma_arf(res, x, t, z, prec);
}

void
arb_fma_si(arb_t res, const arb_t x, slong y, const arb_t z, slong prec)
{
    arf_t t;
    arf_init_set_si(t, y); /* no need to free */
    arb_fma_arf(res, x, t, z, prec);
}

void
arb_fma_fmpz(arb_t res, const arb_t x, const fmpz_t y, const arb_t z, slong prec)
{
    arf_t t;

    if (!COEFF_IS_MPZ(*y))
    {
        arf_init_set_si(t, *y); /* no need to free */
        arb_fma_arf(res, x, t, z, prec);
    }
    else
    {
        arf_init(t);
        arf_set_fmpz(t, y);
        arb_fma_arf(res, x, t, z, prec);
        arf_clear(t);
    }
}

