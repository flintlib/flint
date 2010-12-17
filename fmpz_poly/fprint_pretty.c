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
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void _fmpz_poly_fprint_pretty(FILE * file, 
                              const fmpz * poly, long len, const char * x)
{
    long i;

    if (len == 0)
    {
        fputc('0', file);
        return;
    }
    else if (len == 1)
    {
        fmpz_fprint(file, poly);
        return;
    }

    i = len - 1;

    if (*(poly + i) == 1)
        fprintf(file, "%s^%ld", x, i);
    else if (*(poly + i) == -1)
        fprintf(file, "-%s^%ld", x, i);
    else
    {
        fmpz_fprint(file, poly + i);
        fprintf(file, "*%s^%ld", x, i);
    }

    for ( ; i > 1; --i)
    {
        if (*(poly + i) == 0)
            continue;

        if (*(poly + i) == 1)
            fprintf(file, "+%s^%ld", x, i);
        else if (*(poly + i) == -1)
            fprintf(file, "-%s^%ld", x, i);
        else
        {
            if (fmpz_sgn(poly + i) > 0)
                fputc('+', file);
            fmpz_fprint(file, poly + i);
            fprintf(file, "*%s^%ld", x, i);
        }
    }

    if (*(poly + 1))
    {
        if (*(poly + 1) == 1)
        {
            fputc('+', file);
            fputs(x, file);
        }
        else if (*(poly + 1) == -1)
        {
            fputc('-', file);
            fputs(x, file);
        }
        else
        {
            if (fmpz_sgn(poly + 1) > 0)
                fputc('+', file);
            fmpz_fprint(file, poly + 1);
            fputc('*', file);
            fputs(x, file);
        }
    }
    if (*(poly))
    {
        if (fmpz_sgn(poly) > 0)
            fputc('+', file);
        fmpz_fprint(file, poly);
    }
}

void fmpz_poly_fprint_pretty(FILE * file, 
                             const fmpz_poly_t poly, const char * x)
{
    _fmpz_poly_fprint_pretty(file, poly->coeffs, poly->length, x);
}

