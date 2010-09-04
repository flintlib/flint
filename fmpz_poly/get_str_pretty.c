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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

char *
_fmpz_poly_get_str_pretty(const fmpz * poly, long len, const char *x)
{
    long i, bound, nz;
    fmpz *top;
    char *str, *strbase;

    if (len == 0)
    {
        str = (char *) malloc(2 * sizeof(char));
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
    bound += nz * (3 + strlen(x) + (long) (ceil(log10((double) len))));

    strbase = (char *) malloc(bound * sizeof(char));
    str = strbase;
    top = (fmpz *) poly + (--len);

    if (*top == 1L)
    {
    }
    else if (*top == -1L)
        *str++ = '-';
    else if (!COEFF_IS_MPZ(*top))
        str += sprintf(str, "%li*", *top);
    else
        str += gmp_sprintf(str, "%Zd*", COEFF_TO_PTR(*top));
    str += sprintf(str, "%s^%li", x, len);

    while (--len, --top != poly)
    {
        if (*top == 0L)
            continue;
        if (fmpz_sgn(top) > 0)
            *str++ = '+';
        if (*top == -1L)
            *str++ = '-';
        if (*top != 1L && *top != -1L)
        {
            if (!COEFF_IS_MPZ(*top))
                str += sprintf(str, "%li*", *top);
            else
                str += gmp_sprintf(str, "%Zd*", COEFF_TO_PTR(*top));
        }
        str += sprintf(str, "%s^%li", x, len);
    }

    if (*top != 0L)
    {
        if (fmpz_sgn(top) > 0)
            *str++ = '+';
        if (!COEFF_IS_MPZ(*top))
            str += sprintf(str, "%li", *top);
        else
            str += gmp_sprintf(str, "%Zd", COEFF_TO_PTR(*top));
    }

    return strbase;
}

char *
fmpz_poly_get_str_pretty(const fmpz_poly_t poly, const char *x)
{
    return _fmpz_poly_get_str_pretty(poly->coeffs, poly->length, x);
}
