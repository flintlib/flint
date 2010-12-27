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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

/*
    Recall the return value conventions for fputc (of type int) 

    ``If there are no errors, the same character that has been written is 
    returned.  If an error occurs, EOF is returned and the error indicator 
    is set''

    where the EOF macro expands to a negative int, and fprintf (of type int)

    ``On success, the total number of characters written is returned.
    On failure, a negative number is returned.''
 */

int _fmpz_poly_fprint(FILE * file, const fmpz * poly, long len)
{
    long i;
    int r;

    r = fprintf(file, "%li", len);
    if ((len > 0) && (r > 0))
    {
        r = fputc(' ', file);
        for (i = 0; (i < len) && (r > 0); i++)
        {
            r = fputc(' ', file);
            if (r > 0)
                r = fmpz_fprint(file, poly + i);
        }
    }

    return r;
}

int fmpz_poly_fprint(FILE * file, const fmpz_poly_t poly)
{
    return _fmpz_poly_fprint(file, poly->coeffs, poly->length);
}

