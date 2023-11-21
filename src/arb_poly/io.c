/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2015 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "arb_poly.h"

/* printing *******************************************************************/

void
arb_poly_fprintd(FILE * file, const arb_poly_t poly, slong digits)
{
    slong i;

    flint_fprintf(file, "[");

    for (i = 0; i < poly->length; i++)
    {
        flint_fprintf(file, "(");
        arb_fprintd(file, poly->coeffs + i, digits);
        flint_fprintf(file, ")");

        if (i + 1 < poly->length)
            flint_fprintf(file, ", ");
    }

    flint_fprintf(file, "]");
}

void arb_poly_printd(const arb_poly_t poly, slong digits) { arb_poly_fprintd(stdout, poly, digits); }
