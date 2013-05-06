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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <gmp.h>

#include "fmpz.h"
#include "fmpz_mod_poly.h"

int _fmpz_mod_poly_fprint(FILE * file, const fmpz *poly, long len, 
                          const fmpz_t p)
{
    int r;
    long i;

    r = fprintf(file, "%ld ", len);
    if (r <= 0)
        return r;

    r = fmpz_fprint(file, p);
    if (r <= 0)
        return r;

    if (len == 0)
        return r;

    r = fprintf(file, " ");
    if (r <= 0)
        return r;

    for (i = 0; (r > 0) && (i < len); i++)
    {
        r = fprintf(file, " ");
        if (r <= 0)
            return r;
        r = fmpz_fprint(file, poly + i);
        if (r <= 0)
            return r;
    }

    return r;
}

int fmpz_mod_poly_fprint(FILE * file, const fmpz_mod_poly_t poly)
{
    return _fmpz_mod_poly_fprint(file, poly->coeffs, poly->length, &(poly->p));
}

