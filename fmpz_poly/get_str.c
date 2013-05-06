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
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

char *
_fmpz_poly_get_str(const fmpz * poly, long len)
{
    long i, bound;
    char *str, *strbase;

    if (len == 0)
    {
        str = (char *) flint_malloc(2 * sizeof(char));
        str[0] = '0';
        str[1] = '\0';
        return str;
    }

    bound = (long) (ceil(log10((double) (len + 1))));
    for (i = 0; i < len; i++)
        bound += fmpz_sizeinbase(poly + i, 10) + 1;
    bound += len + 2;

    strbase = (char *) flint_malloc(bound * sizeof(char));
    str = strbase;

    str += sprintf(str, "%li ", len);
    do
    {
        if (!COEFF_IS_MPZ(*poly))
            str += sprintf(str, " %li", *poly);
        else
            str += gmp_sprintf(str, " %Zd", COEFF_TO_PTR(*poly));
    } while (poly++, --len);

    return strbase;
}

char *
fmpz_poly_get_str(const fmpz_poly_t poly)
{
    return _fmpz_poly_get_str(poly->coeffs, poly->length);
}
