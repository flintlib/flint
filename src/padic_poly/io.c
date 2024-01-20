/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz_vec.h"
#include "padic.h"
#include "padic_poly.h"

int _padic_poly_fprint(FILE *file, const fmpz *poly, slong val, slong len,
                       const padic_ctx_t ctx)
{
    slong i, v;
    fmpz_t u;

    if (len == 0)
    {
        flint_fprintf(file, "0");
        return 1;
    }

    fmpz_init(u);

    flint_fprintf(file, "%wd ", len);

    for (i = 0; i < len; i++)
    {
        flint_fprintf(file, " ");

        if (fmpz_is_zero(poly + i))
        {
            flint_fprintf(file, "0");
        }
        else
        {
            v = val + fmpz_remove(u, poly + i, ctx->p);

            _padic_fprint(file, u, v, ctx);
        }
    }

    fmpz_clear(u);

    return 1;
}

int _padic_poly_fprint_pretty(FILE *file,
                              const fmpz *poly, slong len, slong val,
                              const char *var,
                              const padic_ctx_t ctx)
{
    slong i;
    padic_t x;

    padic_init(x);  /* Precision is not used anywhere! */

    if (len == 0)
    {
        fputc('0', file);
    }
    else if (len == 1)
    {
        _padic_fprint(file, poly + 0, val, ctx);
    }
    else if (len == 2)
    {
        fmpz_set(padic_unit(x), poly + 1);
        padic_val(x) = val;
        _padic_canonicalise(x, ctx);

        if (padic_is_one(x))
        {
            flint_fprintf(file, "%s", var);
        }
        else if (*(padic_unit(x)) == WORD(-1) && padic_val(x) == 0)
        {
            flint_fprintf(file, "-%s", var);
        }
        else
        {
            fputc('(', file);
            padic_fprint(file, x, ctx);
            fputc(')', file);
            flint_fprintf(file, "*%s", var);
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

            if (padic_is_one(x))
               flint_fprintf(file, "%s^%wd", var, i);
            else if (*(padic_unit(x)) == WORD(-1) && padic_val(x) == 0)
               flint_fprintf(file, "-%s^%wd", var, i);
            else
            {
                fputc('(', file);
                padic_fprint(file, x, ctx);
                fputc(')', file);
                flint_fprintf(file, "*%s^%wd", var, i);
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

            if (padic_is_one(x))
                flint_fprintf(file, "%s^%wd", var, i);
            else
            {
                fputc('(', file);
                padic_fprint(file, x, ctx);
                fputc(')', file);
                flint_fprintf(file, "*%s^%wd", var, i);
            }
        }

        if (*(poly + 1))
        {
            fmpz_abs(padic_unit(x), poly + 1);
            padic_val(x) = val;
            _padic_canonicalise(x, ctx);

            fputc(fmpz_sgn(poly + 1) > 0 ? '+' : '-', file);

            if (padic_is_one(x))
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

    padic_clear(x);
    return 1;
}

int padic_poly_fprint(FILE *file, const padic_poly_t poly, const padic_ctx_t ctx) { return _padic_poly_fprint(file, poly->coeffs, poly->val, poly->length, ctx); }
int padic_poly_fprint_pretty(FILE * file, const padic_poly_t poly, const char * var, const padic_ctx_t ctx) { return _padic_poly_fprint_pretty(file, poly->coeffs, poly->length, poly->val, var, ctx); }
int _padic_poly_print(const fmpz *poly, slong val, slong len, const padic_ctx_t ctx) { return _padic_poly_fprint(stdout, poly, val, len, ctx); }
int padic_poly_print(const padic_poly_t poly, const padic_ctx_t ctx) { return padic_poly_fprint(stdout, poly, ctx); }
int _padic_poly_print_pretty(const fmpz *poly, slong val, slong len, const char *var, const padic_ctx_t ctx) { return _padic_poly_fprint_pretty(stdout, poly, val, len, var, ctx); }
int padic_poly_print_pretty(const padic_poly_t poly, const char *var, const padic_ctx_t ctx) { return padic_poly_fprint_pretty(stdout, poly, var, ctx); }

int padic_poly_debug(const padic_poly_t poly)
{
    flint_printf("(alloc = %wd, length = %wd, val = %wd, N = %wd, vec = ",
        poly->alloc, poly->length, poly->val, poly->N);
    if (poly->coeffs)
    {
        flint_printf("{");
        _fmpz_vec_print(poly->coeffs, poly->alloc);
        flint_printf("}");
    }
    else
    {
        flint_printf("NULL");
    }
    flint_printf(")");

    return 1;
}
