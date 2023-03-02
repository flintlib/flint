/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <string.h>
#include <math.h>

char *_TEMPLATE(T, poly_get_str_pretty) (const TEMPLATE(T, struct) * poly,
                                         slong len, const char *x,
                                         const TEMPLATE(T, ctx_t) ctx)
{
    char *str, **coeffstrs;
    size_t off;
    slong i, bound, nz;

    if (len == 0)
    {
        str = flint_malloc(2);
        str[0] = '0';
        str[1] = '\0';
        return str;
    }

    if (len == 1)
    {
        return TEMPLATE(T, get_str_pretty) (poly, ctx);
    }

    coeffstrs = (char **)flint_malloc(len * sizeof(char *));

    nz = 0;
    bound = 1;
    for (i = 0; i < len; i++)
        if (!TEMPLATE(T, is_zero) (poly + i, ctx))
        {
            coeffstrs[i] = TEMPLATE(T, get_str_pretty) (poly + i, ctx);
            bound += strlen(coeffstrs[i]);
            nz++;
        }
    bound += nz * (5 + strlen(x) + (slong) (ceil(log10(len))));

    str = flint_malloc(bound);
    off = 0;
    i = len - 1;

    if (TEMPLATE(T, is_one) (poly + i, ctx))
    {
    }
    else
        off += flint_sprintf(str + off, "(%s)*", coeffstrs[i]);

    if (i > 1)
        off += flint_sprintf(str + off, "%s^%wd", x, i);
    else
        off += flint_sprintf(str + off, "%s", x);

    for (--i; i > 0; --i)
    {
        if (TEMPLATE(T, is_zero) (poly + i, ctx))
            continue;
        if (!TEMPLATE(T, is_one) (poly + i, ctx))
        {
            off += flint_sprintf(str + off, "+(%s)*", coeffstrs[i]);
        }
        else
        {
            off += flint_sprintf(str + off, "+");
        }
        if (i > 1)
            off += flint_sprintf(str + off, "%s^%wd", x, i);
        else
            off += flint_sprintf(str + off, "%s", x);
    }

    if (!TEMPLATE(T, is_zero) (poly + i, ctx))
    {
        off += flint_sprintf(str + off, "+(%s)", coeffstrs[i]);
    }

    for (i = 0; i < len; i++)
        if (!TEMPLATE(T, is_zero) (poly + i, ctx))
        {
            flint_free(coeffstrs[i]);
        }

    flint_free(coeffstrs);

    return str;
}

char *TEMPLATE(T, poly_get_str_pretty) (const TEMPLATE(T, poly_t) poly,
                                        const char *x,
                                        const TEMPLATE(T, ctx_t) ctx)
{
    return _TEMPLATE(T, poly_get_str_pretty) (poly->coeffs, poly->length, x,
                                              ctx);
}


#endif
