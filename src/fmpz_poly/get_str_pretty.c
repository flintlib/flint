/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

char *
_fmpz_poly_get_str_pretty(const fmpz * poly, slong len, const char *x)
{
    char *str;
    size_t off;
    slong i, bound, nz;

    if (len == 0)
    {
        str = flint_malloc(2);
        str[0] = '0';
        str[1] = '\0';
        return str;
    }

    if (len == 1)
    {
        str = fmpz_get_str(NULL, 10, poly);
        return str;
    }

    nz = 0;
    bound = 1;
    for (i = 0; i < len; i++)
        if (!fmpz_is_zero(poly + i))
        {
            bound += fmpz_sizeinbase(poly + i, 10) + 1;
            nz++;
        }
    bound += nz * (3 + strlen(x) + (slong) (ceil(log10(len))));

    str = flint_malloc(bound);
    off = 0;
    i = len - 1;

    if (poly[i] == WORD(1))
    {
    }
    else if (poly[i] == WORD(-1))
        str[off++] = '-';
    else if (!COEFF_IS_MPZ(poly[i]))
        off += flint_sprintf(str + off, "%wd*", poly[i]);
    else
        off += gmp_sprintf(str + off, "%Zd*", COEFF_TO_PTR(poly[i]));
    if (i > 1)
        off += flint_sprintf(str + off, "%s^%wd", x, i);
    else
        off += flint_sprintf(str + off, "%s", x);

    for (--i; i > 0; --i)
    {
        if (poly[i] == WORD(0))
            continue;
        if (fmpz_sgn(poly + i) > 0)
            str[off++] = '+';
        if (poly[i] == WORD(-1))
            str[off++] = '-';
        if (poly[i] != WORD(1) && poly[i] != WORD(-1))
        {
            if (!COEFF_IS_MPZ(poly[i]))
                off += flint_sprintf(str + off, "%wd*", poly[i]);
            else
                off += gmp_sprintf(str + off, "%Zd*", COEFF_TO_PTR(poly[i]));
        }
        if (i > 1)
            off += flint_sprintf(str + off, "%s^%wd", x, i);
        else
            off += flint_sprintf(str + off, "%s", x);
    }

    if (poly[i] != WORD(0))
    {
        if (fmpz_sgn(poly + i) > 0)
            str[off++] = '+';
        if (!COEFF_IS_MPZ(poly[i]))
            off += flint_sprintf(str + off, "%wd", poly[i]);
        else
            off += gmp_sprintf(str + off, "%Zd", COEFF_TO_PTR(poly[i]));
    }

    return str;
}

char *
fmpz_poly_get_str_pretty(const fmpz_poly_t poly, const char *x)
{
    return _fmpz_poly_get_str_pretty(poly->coeffs, poly->length, x);
}
