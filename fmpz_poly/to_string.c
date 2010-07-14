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

char * _fmpz_poly_to_string(const fmpz * poly, ulong len)
{
    ulong i, bound;
    char * str, * strbase;
    
    if (len == 0UL)
    {
        str = malloc(2 * sizeof(char));
        str[0] = '0';
        str[1] = '\0';
        return str;
    }
    
    bound  = (ulong) (ceil(log10((double) (len + 1UL))));
    for (i = 0; i < len; i++)
        bound += fmpz_sizeinbase(poly + i, 10) + 1UL;
    bound += len + 1UL;
    
    strbase = malloc(bound * sizeof(char));
    str = strbase;
    
    str += sprintf(str, "%lu ", len);
    do {
        if (!COEFF_IS_MPZ(*poly))
            str += sprintf(str, " %li", *poly);
        else
            str += gmp_sprintf(str, " %Zd", COEFF_TO_PTR(*poly));
    } while (poly++, --len);
    
    return strbase;
}

char * fmpz_poly_to_string(const fmpz_poly_t poly)
{
    return _fmpz_poly_to_string(poly->coeffs, poly->length);
}

