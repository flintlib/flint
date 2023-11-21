/*
    Copyright (C) 2012, 2015 Fredrik Johansson
    Copyright (C) 2015 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "acb.h"

/* printing *******************************************************************/

void
acb_fprint(FILE * file, const acb_t x)
{
    flint_fprintf(file, "(");
    arb_fprint(file, acb_realref(x));
    flint_fprintf(file, ", ");
    arb_fprint(file, acb_imagref(x));
    flint_fprintf(file, ")");
}

void
acb_fprintd(FILE * file, const acb_t z, slong digits)
{
    flint_fprintf(file, "(");
    arf_fprintd(file, arb_midref(acb_realref(z)), digits);

    if (arf_sgn(arb_midref(acb_imagref(z))) < 0)
    {
        arf_t t;
        arf_init_neg_shallow(t, arb_midref(acb_imagref(z)));
        flint_fprintf(file, " - ");
        arf_fprintd(file, t, digits);
    }
    else
    {
        flint_fprintf(file, " + ");
        arf_fprintd(file, arb_midref(acb_imagref(z)), digits);
    }
    flint_fprintf(file, "j)");

    flint_fprintf(file, "  +/-  ");

    flint_fprintf(file, "(");
    mag_fprintd(file, arb_radref(acb_realref(z)), 3);
    flint_fprintf(file, ", ");
    mag_fprintd(file, arb_radref(acb_imagref(z)), 3);
    flint_fprintf(file, "j)");
}

void
acb_fprintn(FILE * file, const acb_t z, slong digits, ulong flags)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        arb_fprintn(file, acb_realref(z), digits, flags);
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        arb_fprintn(file, acb_imagref(z), digits, flags);
        flint_fprintf(file, "*I");
    }
    else
    {
        arb_fprintn(file, acb_realref(z), digits, flags);

        if ((arb_is_exact(acb_imagref(z)) || (flags & ARB_STR_NO_RADIUS))
                && arf_sgn(arb_midref(acb_imagref(z))) < 0)
        {
            arb_t t;
            arb_init(t);
            arb_neg(t, acb_imagref(z));
            flint_fprintf(file, " - ");
            arb_fprintn(file, t, digits, flags);
            arb_clear(t);
        }
        else
        {
            flint_fprintf(file, " + ");
            arb_fprintn(file, acb_imagref(z), digits, flags);
        }

        flint_fprintf(file, "*I");
    }
}

void acb_print(const acb_t x) { acb_fprint(stdout, x); }
void acb_printd(const acb_t z, slong digits) { acb_fprintd(stdout, z, digits); }
void acb_printn(const acb_t x, slong digits, ulong flags) { acb_fprintn(stdout, x, digits, flags); }

void
_acb_vec_printn(acb_srcptr vec, slong len, slong ndigits, ulong flags)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        acb_printn(vec + i, ndigits, flags);
        if (i < len - 1)
            flint_printf(", ");
    }
}

void
_acb_vec_printd(acb_srcptr vec, slong len, slong ndigits)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        acb_printd(vec + i, ndigits);
        if (i < len - 1)
        {
            flint_printf(", ");
        }
    }
    flint_printf("\n");
}
