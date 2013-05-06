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
#include <stdio.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

/*
    Recall the return value conventions for fputc (of type int) 

    ``If there are no errors, the same character that has been written is 
    returned.  If an error occurs, EOF is returned and the error indicator 
    is set''

    where the EOF macro expands to a negative int, and fprintf (of type int)

    ``On success, the total number of characters written is returned.
    On failure, a negative number is returned.''
 */

int 
_fmpq_poly_fprint(FILE * file, const fmpz * poly, const fmpz_t den, long len)
{
    int r;
    long i;
    fmpz_t n, d, g;

    fmpz_init(n);
    fmpz_init(d);
    fmpz_init(g);

    r = fprintf(file, "%li", len);
    if ((len > 0) && (r > 0))
    {
        r = fputc(' ', file);
        for (i = 0; (i < len) && (r > 0); i++)
        {
            r = fputc(' ', file);
            if (r > 0)
            {
                fmpz_gcd(g, poly + i, den);
                fmpz_divexact(n, poly + i, g);
                fmpz_divexact(d, den, g);
                if (*d == 1L)
                    r = fmpz_fprint(file, n);
                else
                {
                    r = fmpz_fprint(file, n);
                    if (r > 0)
                        r = fputc('/', file);
                    if (r > 0)
                        r = fmpz_fprint(file, d);
                }
            }
        }
    }

    fmpz_clear(n);
    fmpz_clear(d);
    fmpz_clear(g);

    return r;
}

int fmpq_poly_fprint(FILE * file, const fmpq_poly_t poly)
{
    return _fmpq_poly_fprint(file, poly->coeffs, poly->den, poly->length);
}

