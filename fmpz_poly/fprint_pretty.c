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

#include "fmpz_poly.h"

int _fmpz_poly_fprint_pretty(FILE * file, 
                             const fmpz * poly, len_t len, const char * x)
{
    int r;
    len_t i;

    FMPZ_VEC_NORM(poly, len);

    if (len == 0)
    {
        r = fputc('0', file);
        r = (r != EOF) ? 1 : EOF;
        return r;
    }
    else if (len == 1)
    {
        r = fmpz_fprint(file, poly + 0);
        return r;
    }
    else if (len == 2)
    {
        if (*(poly + 1) == 1L)
        {
            r = fprintf(file, "%s", x);
        }
        else if (*(poly + 1) == -1L)
        {
            r = fprintf(file, "-%s", x);
        }
        else
        {
            r = fmpz_fprint(file, poly + 1);
            if (r > 0)
                r = fprintf(file, "*%s", x);
        }
        
        if (r > 0)
        {
            if (fmpz_sgn(poly + 0) > 0)
            {
                r = fprintf(file, "+");
                if (r > 0)
                    r = fmpz_fprint(file, poly + 0);
            }
            else if (fmpz_sgn(poly + 0) < 0)
            {
                r = fmpz_fprint(file, poly + 0);
            }
        }
        return r;
    }

    i = len - 1;  /* i >= 2 */
    r = 1;
    {
        if (*(poly + i) == 1)
           r = fprintf(file, "%s^%ld", x, i);
        else if (*(poly + i) == -1)
           r = fprintf(file, "-%s^%ld", x, i);
        else
        {
           r = fmpz_fprint(file, poly + i);
           if (r > 0)
              r = fprintf(file, "*%s^%ld", x, i);
        }
        --i;
    }

    for (; (r > 0) && (i > 1); --i)
    {
        if (*(poly + i) == 0)
            continue;

        if (*(poly + i) == 1)
            r = fprintf(file, "+%s^%ld", x, i);
        else if (*(poly + i) == -1)
            r = fprintf(file, "-%s^%ld", x, i);
        else
        {
            if (fmpz_sgn(poly + i) > 0)
            {
                r = fputc('+', file);
                r = (r != EOF) ? 1 : EOF;
            }
            if (r > 0)
                r = fmpz_fprint(file, poly + i);
            if (r > 0)
                r = fprintf(file, "*%s^%ld", x, i);
        }
    }

    if ((r > 0) && *(poly + 1))
    {
        if (*(poly + 1) == 1)
        {
            r = fputc('+', file);
            r = (r != EOF) ? 1 : EOF;
            if (r > 0)
            {
                r = fputs(x, file);
                r = (r >= 0) ? 1 : -1;
            }
        }
        else if (*(poly + 1) == -1)
        {
            r = fputc('-', file);
            r = (r != EOF) ? 1 : EOF;
            if (r > 0)
            {
                r = fputs(x, file);
                r = (r >= 0) ? 1 : -1;
            }
        }
        else
        {
            if (fmpz_sgn(poly + 1) > 0)
            {
                r = fputc('+', file);
                r = (r != EOF) ? 1 : EOF;
            }
            if (r > 0)
                r = fmpz_fprint(file, poly + 1);
            if (r > 0)
            {
                r = fputc('*', file);
                r = (r != EOF) ? 1 : EOF;
            }
            if (r > 0)
            {
                r = fputs(x, file);
                r = (r >= 0) ? 1 : -1;
            }
        }
    }
    if ((r > 0) && *(poly))
    {
        if (fmpz_sgn(poly) > 0)
        {
            r = fputc('+', file);
            r = (r != EOF) ? 1 : EOF;
        }
        if (r > 0)
            r = fmpz_fprint(file, poly);
    }

    return r;
}

int fmpz_poly_fprint_pretty(FILE * file, 
                            const fmpz_poly_t poly, const char * x)
{
    return _fmpz_poly_fprint_pretty(file, poly->coeffs, poly->length, x);
}

