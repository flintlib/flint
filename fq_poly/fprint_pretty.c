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

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include "fq_poly.h"

/*
    TODO:  When op has only one non-zero term, this 
    function would omit the parentheses.
 */

static void __fq_print(FILE *file, const fq_struct *op, const fq_ctx_t ctx)
{
    fputc('(', file);
    fq_fprint_pretty(file, op, ctx);
    fputc(')', file);
}

int _fq_poly_fprint_pretty(FILE *file, const fq_struct *poly, long len, 
                           const char *x, const fq_ctx_t ctx)
{
    if (len == 0)
    {
        fputc('0', file);
    }
    else if (len == 1)
    {
        fq_fprint_pretty(file, poly + 0, ctx);
    }
    else if (len == 2)
    {
        if (fq_is_one(poly + 1))
            fprintf(file, "%s", x);
        else
        {
            __fq_print(file, poly + 1, ctx);
            fprintf(file, "*%s", x);
        }
        if (!fq_is_zero(poly + 0))
        {
            fputc('+', file);
            __fq_print(file, poly + 0, ctx);
        }
    }
    else
    {
        long i = len - 1;

        {
            if (fq_is_one(poly + i))
                fprintf(file, "%s^%ld", x, i);
            else
            {
                __fq_print(file, poly + i, ctx);
                fprintf(file, "*%s^%ld", x, i);
            }
            --i;
        }

        for (; i > 1; --i)
        {
            if (fq_is_zero(poly + i))
                continue;

            if (fq_is_one(poly + i))
                fprintf(file, "+%s^%ld", x, i);
            else
            {
                fputc('+', file);
                __fq_print(file, poly + i, ctx);
                fprintf(file, "*%s^%ld", x, i);
            }
        }

        if (!fq_is_zero(poly + 1))
        {
            if (fq_is_one(poly + 1))
            {
                fputc('+', file);
                fputs(x, file);
            }
            else
            {
                fputc('+', file);
                __fq_print(file, poly + 1, ctx);
                fputc('*', file);
                fputs(x, file);
            }
        }
        if (!fq_is_zero(poly + 0))
        {
            fputc('+', file);
            __fq_print(file, poly + 0, ctx);
        }
    }

    return 1;
}

int fq_poly_fprint_pretty(FILE *file, const fq_poly_t poly, const char *x, 
                          const fq_ctx_t ctx)
{
    return _fq_poly_fprint_pretty(file, poly->coeffs, poly->length, x, ctx);
}

