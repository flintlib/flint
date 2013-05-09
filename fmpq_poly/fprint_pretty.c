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

    Copyright (C) 2010, 2012 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "fmpq_poly.h"

/*
    Macro wrapping _fmpq_fprint(file, x, y), ensuring that the printed 
    rational is in lowest terms.  Assumes that y > 0.
 */

#define __fmpq_fprint(x,y)         \
do {                               \
    fmpz_gcd(g, x, y);             \
    if (fmpz_is_one(g))            \
    {                              \
        _fmpq_fprint(file, x, y);  \
    }                              \
    else                           \
    {                              \
        fmpz_divexact(n, x, g);    \
        fmpz_divexact(d, y, g);    \
        _fmpq_fprint(file, n, d);  \
    }                              \
} while (0)

int _fmpq_poly_fprint_pretty(FILE * file, 
                             const fmpz *poly, const fmpz_t den, len_t len, 
                             const char * x)
{
    fmpz_t n, d, g;

    fmpz_init(n);
    fmpz_init(d);
    fmpz_init(g);

    if (len == 0)
    {
        fputc('0', file);
    }
    else if (len == 1)
    {
        _fmpq_fprint(file, poly + 0, den);
    }
    else if (len == 2)
    {
        if (poly[1] == 1L)
        {
            fprintf(file, "%s", x);
        }
        else if (poly[1] == -1L)
        {
            fprintf(file, "-%s", x);
        }
        else
        {
            __fmpq_fprint(poly + 1, den);
            fprintf(file, "*%s", x);
        }
        
        if (fmpz_sgn(poly + 0) > 0)
        {
            fprintf(file, "+");
            __fmpq_fprint(poly + 0, den);
        }
        else if (fmpz_sgn(poly + 0) < 0)
        {
            __fmpq_fprint(poly + 0, den);
        }
    }
    else  /* len >= 3 */
    {
        len_t i = len - 1;  /* i >= 2 */
        {
            if (poly[i] == 1L)
               fprintf(file, "%s^%ld", x, i);
            else if (poly[i] == -1L)
               fprintf(file, "-%s^%ld", x, i);
            else
            {
               __fmpq_fprint(poly + i, den);
               fprintf(file, "*%s^%ld", x, i);
            }
            --i;
        }

        for (; i > 1; --i)
        {
            if (poly[i] == 0)
                continue;

            if (poly[i] == 1L)
                fprintf(file, "+%s^%ld", x, i);
            else if (poly[i] == -1L)
                fprintf(file, "-%s^%ld", x, i);
            else
            {
                if (fmpz_sgn(poly + i) > 0)
                {
                    fputc('+', file);
                }
                __fmpq_fprint(poly + i, den);
                fprintf(file, "*%s^%ld", x, i);
            }
        }

        if (poly[1])
        {
            if (poly[1] == 1L)
            {
                fputc('+', file);
                fputs(x, file);
            }
            else if (poly[1] == -1L)
            {
                fputc('-', file);
                fputs(x, file);
            }
            else
            {
                if (fmpz_sgn(poly + 1) > 0)
                {
                    fputc('+', file);
                }
                __fmpq_fprint(poly + 1, den);
                fputc('*', file);
                fputs(x, file);
            }
        }
        if (*(poly))
        {
            if (fmpz_sgn(poly) > 0)
            {
                fputc('+', file);
            }
            __fmpq_fprint(poly + 0, den);
        }
    }

    fmpz_clear(n);
    fmpz_clear(d);
    fmpz_clear(g);

    return 1;
}

#undef __fmpq_fprint

int fmpq_poly_fprint_pretty(FILE * file, 
                            const fmpq_poly_t poly, const char * var)
{
    return _fmpq_poly_fprint_pretty(file, poly->coeffs, poly->den, poly->length, var);
}

