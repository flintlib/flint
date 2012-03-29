/*============================================================================

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

#include "padic_poly.h"

int _padic_poly_fprint_pretty(FILE *file, 
                              const fmpz *poly, long len, long val, 
                              const char *var, 
                              const padic_ctx_t ctx)
{
    long i;
    padic_t x;

    _padic_init(x);

    if (len == 0)
    {
        fputc('0', file);
    }
    else if (len == 1)
    {
        fmpz_set(padic_unit(x), poly);
        padic_val(x) = val;

        padic_fprint(file, x, ctx);
    }
    else if (len == 2)
    {
        fmpz_set(padic_unit(x), poly + 1);
        padic_val(x) = val;
        _padic_canonicalise(x, ctx);

        if (_padic_is_one(x))
        {
            fprintf(file, "%s", var);
        }
        else if (*(padic_unit(x)) == -1L && padic_val(x) == 0)
        {
            fprintf(file, "-%s", var);
        }
        else
        {
            fputc('(', file);
            padic_fprint(file, x, ctx);
            fputc(')', file);
            fprintf(file, "*%s", var);
        }

        fmpz_abs(padic_unit(x), poly);
        padic_val(x) = val;
        _padic_canonicalise(x, ctx);

        if (fmpz_sgn(poly) > 0)
        {
            fputc('+', file);
        }
        else if (fmpz_sgn(poly) < 0)
        {
            fputc('-', file);
        }
        fputc('(', file);
        padic_fprint(file, x, ctx);
        fputc(')', file);
    }
    else  /* len >= 3 */
    {
        i = len - 1;  /* i >= 2 */
        {
            fmpz_set(padic_unit(x), poly + i);
            padic_val(x) = val;
            _padic_canonicalise(x, ctx);

            if (_padic_is_one(x))
               fprintf(file, "%s^%ld", var, i);
            else if (*(padic_unit(x)) == -1L && padic_val(x) == 0)
               fprintf(file, "-%s^%ld", var, i);
            else
            {
                fputc('(', file);
                padic_fprint(file, x, ctx);
                fputc(')', file);
               fprintf(file, "*%s^%ld", var, i);
            }
            --i;
        }

        for (; i > 1; --i)
        {
            if (*(poly + i) == 0)
                continue;

            fmpz_abs(padic_unit(x), poly + i);
            padic_val(x) = val;
            _padic_canonicalise(x, ctx);

            if (fmpz_sgn(poly + i) > 0)
                fputc('+', file);
            else
                fputc('-', file);

            if (_padic_is_one(x))
                fprintf(file, "%s^%ld", var, i);
            else
            {
                fputc('(', file);
                padic_fprint(file, x, ctx);
                fputc(')', file);
                fprintf(file, "*%s^%ld", var, i);
            }
        }

        if (*(poly + 1))
        {
            fmpz_abs(padic_unit(x), poly + 1);
            padic_val(x) = val;
            _padic_canonicalise(x, ctx);

            fputc(fmpz_sgn(poly + 1) > 0 ? '+' : '-', file);

            if (_padic_is_one(x))
                fputs(var, file);
            else
            {
                fputc('(', file);
                padic_fprint(file, x, ctx);
                fputc(')', file);
                fputc('*', file);
                fputs(var, file);
            }
        }
        if (*(poly))
        {
            fmpz_abs(padic_unit(x), poly);
            padic_val(x) = val;
            _padic_canonicalise(x, ctx);

            fputc(fmpz_sgn(poly) > 0 ? '+' : '-', file);
            fputc('(', file);
            padic_fprint(file, x, ctx);
            fputc(')', file);
        }
    }

    _padic_clear(x);

    return 1;
}

int padic_poly_fprint_pretty(FILE *file, 
                             const padic_poly_t poly, const char *var, 
                             const padic_ctx_t ctx)
{
    return _padic_poly_fprint_pretty(file, 
        poly->coeffs, poly->length, poly->val, var, ctx);
}

