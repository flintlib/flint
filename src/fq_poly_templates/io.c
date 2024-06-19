/*
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"
#include "fmpz.h"

/* printing *******************************************************************/

/*
    TODO:  When op has only one non-zero term, this
    function would omit the parentheses.
 */

static void
__TEMPLATE(T, print)(FILE * file, const TEMPLATE(T, struct) * op,
                      const TEMPLATE(T, ctx_t) ctx)
{
    fputc('(', file);
    TEMPLATE(T, fprint_pretty) (file, op, ctx);
    fputc(')', file);
}

int
_TEMPLATE(T, poly_fprint)(FILE * file, const TEMPLATE(T, struct) * poly,
                           slong len, const TEMPLATE(T, ctx_t) ctx)
{
    int r;
    slong i;

    r = flint_fprintf(file, "%wd ", len);
    if (r <= 0)
        return r;

    if (len == 0)
        return r;

    for (i = 0; (r > 0) && (i < len); i++)
    {
        r = flint_fprintf(file, " ");
        if (r <= 0)
            return r;
        r = TEMPLATE(T, fprint) (file, poly + i, ctx);
        if (r <= 0)
            return r;
    }

    return r;
}

int
_TEMPLATE(T, poly_fprint_pretty)(FILE * file,
                                  const TEMPLATE(T, struct) * poly, slong len,
                                  const char *x, const TEMPLATE(T, ctx_t) ctx)
{
    if (len == 0)
    {
        fputc('0', file);
    }
    else if (len == 1)
    {
        TEMPLATE(T, fprint_pretty) (file, poly + 0, ctx);
    }
    else if (len == 2)
    {
        if (TEMPLATE(T, is_one) (poly + 1, ctx))
            flint_fprintf(file, "%s", x);
        else
        {
            __TEMPLATE(T, print) (file, poly + 1, ctx);
            flint_fprintf(file, "*%s", x);
        }
        if (!TEMPLATE(T, is_zero) (poly + 0, ctx))
        {
            fputc('+', file);
            __TEMPLATE(T, print) (file, poly + 0, ctx);
        }
    }
    else
    {
        slong i = len - 1;

        {
            if (TEMPLATE(T, is_one) (poly + i, ctx))
                flint_fprintf(file, "%s^%wd", x, i);
            else
            {
                __TEMPLATE(T, print) (file, poly + i, ctx);
                flint_fprintf(file, "*%s^%wd", x, i);
            }
            --i;
        }

        for (; i > 1; --i)
        {
            if (TEMPLATE(T, is_zero) (poly + i, ctx))
                continue;

            if (TEMPLATE(T, is_one) (poly + i, ctx))
                flint_fprintf(file, "+%s^%wd", x, i);
            else
            {
                fputc('+', file);
                __TEMPLATE(T, print) (file, poly + i, ctx);
                flint_fprintf(file, "*%s^%wd", x, i);
            }
        }

        if (!TEMPLATE(T, is_zero) (poly + 1, ctx))
        {
            if (TEMPLATE(T, is_one) (poly + 1, ctx))
            {
                fputc('+', file);
                fputs(x, file);
            }
            else
            {
                fputc('+', file);
                __TEMPLATE(T, print) (file, poly + 1, ctx);
                fputc('*', file);
                fputs(x, file);
            }
        }
        if (!TEMPLATE(T, is_zero) (poly + 0, ctx))
        {
            fputc('+', file);
            __TEMPLATE(T, print) (file, poly + 0, ctx);
        }
    }

    return 1;
}

int TEMPLATE(T, poly_fprint)(FILE * file, const TEMPLATE(T, poly_t) poly, const TEMPLATE(T, ctx_t) ctx) { return _TEMPLATE(T, poly_fprint)(file, poly->coeffs, poly->length, ctx); }
int TEMPLATE(T, poly_fprint_pretty)(FILE * file, const TEMPLATE(T, poly_t) poly, const char *x, const TEMPLATE(T, ctx_t) ctx) { return _TEMPLATE(T, poly_fprint_pretty)(file, poly->coeffs, poly->length, x, ctx); }
int _TEMPLATE(T, poly_print)(const TEMPLATE(T, struct) *poly, slong len, const TEMPLATE(T, ctx_t) ctx) { return _TEMPLATE(T, poly_fprint)(stdout, poly, len, ctx); }
int TEMPLATE(T, poly_print)(const TEMPLATE(T, poly_t) poly, const TEMPLATE(T, ctx_t) ctx) { return TEMPLATE(T, poly_fprint)(stdout, poly, ctx); }
int _TEMPLATE(T, poly_print_pretty)(const TEMPLATE(T, struct) *poly, slong len, const char *x, const TEMPLATE(T, ctx_t) ctx) { return _TEMPLATE(T, poly_fprint_pretty)(stdout, poly, len, x, ctx); }
int TEMPLATE(T, poly_print_pretty)(const TEMPLATE(T, poly_t) poly, const char *x, const TEMPLATE(T, ctx_t) ctx) { return TEMPLATE(T, poly_fprint_pretty)(stdout, poly, x, ctx); }

#endif
