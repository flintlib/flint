/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010, 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

char *
_fmpz_poly_get_str_pretty(const fmpz * poly, len_t len, const char *x)
{
    char *str;
    size_t off;
    len_t i, bound, nz;

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
    bound += nz * (3 + strlen(x) + (len_t) (ceil(log10(len))));

    str = flint_malloc(bound);
    off = 0;
    i = len - 1;

    if (poly[i] == 1L)
    {
    }
    else if (poly[i] == -1L)
        str[off++] = '-';
    else if (!COEFF_IS_MPZ(poly[i]))
        off += sprintf(str + off, "%ld*", poly[i]);
    else
        off += gmp_sprintf(str + off, "%Zd*", COEFF_TO_PTR(poly[i]));
    if (i > 1)
        off += sprintf(str + off, "%s^%ld", x, i);
    else
        off += sprintf(str + off, "%s", x);

    for (--i; i > 0; --i)
    {
        if (poly[i] == 0L)
            continue;
        if (fmpz_sgn(poly + i) > 0)
            str[off++] = '+';
        if (poly[i] == -1L)
            str[off++] = '-';
        if (poly[i] != 1L && poly[i] != -1L)
        {
            if (!COEFF_IS_MPZ(poly[i]))
                off += sprintf(str + off, "%ld*", poly[i]);
            else
                off += gmp_sprintf(str + off, "%Zd*", COEFF_TO_PTR(poly[i]));
        }
        if (i > 1)
            off += sprintf(str + off, "%s^%ld", x, i);
        else
            off += sprintf(str + off, "%s", x);
    }

    if (poly[i] != 0L)
    {
        if (fmpz_sgn(poly + i) > 0)
            str[off++] = '+';
        if (!COEFF_IS_MPZ(poly[i]))
            off += sprintf(str + off, "%ld", poly[i]);
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
