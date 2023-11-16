/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
ca_poly_print(const ca_poly_t poly, ca_ctx_t ctx)
{
    slong i, len;

    len = poly->length;

    flint_printf("ca_poly of length %wd:\n", len);

    for (i = 0; i < len; i++)
    {
        flint_printf("    ");
        ca_print(poly->coeffs + i, ctx);
        flint_printf("\n");
    }

    flint_printf("\n");
}

void
ca_poly_printn(const ca_poly_t poly, slong digits, ca_ctx_t ctx)
{
    slong len, i;

    len = poly->length;

    flint_printf("[");

    for (i = 0; i < len; i++)
    {
        ca_printn(poly->coeffs + i, digits, ctx);
        if (i < len - 1)
            flint_printf(", ");
    }
    flint_printf("]\n");
}
