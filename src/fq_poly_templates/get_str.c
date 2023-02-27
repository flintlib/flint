/*
    Copyright (C) 2010 Sebastian Pancratz
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

char *_TEMPLATE(T, poly_get_str) (const TEMPLATE(T, struct) * poly, slong len,
                                  const TEMPLATE(T, ctx_t) ctx)
{
    char *str, **coeffstrs;
    size_t off;
    slong i, bound, nz;

    if (len == 0)
    {
        str = (char *)flint_malloc(2 * sizeof(char));
        str[0] = '0';
        str[1] = '\0';
        return str;
    }

    coeffstrs = (char **)flint_malloc(len * sizeof(char *));

    nz = 0;
    bound = (slong) (ceil(log10((double)(len + 1)))) + 2;
    for (i = 0; i < len; i++)
    {
        if (!TEMPLATE(T, is_zero) (poly + i, ctx))
        {
            coeffstrs[i] = TEMPLATE(T, get_str) (poly + i, ctx);
            bound += 1 + strlen(coeffstrs[i]);
            nz++;
        }
        else
        {
            bound += 2;
        }
    }

    str = (char *)flint_malloc(bound * sizeof(char));
    off = 0;

    off += flint_sprintf(str + off, "%wd ", len);
    for (i = 0; i < len; i++)
    {
        if (TEMPLATE(T, is_zero) (poly + i, ctx))
        {
            off += flint_sprintf(str + off, " 0");
        }
        else
        {
            off += flint_sprintf(str + off, " %s", coeffstrs[i]);
            flint_free(coeffstrs[i]);
        }
    }

    flint_free(coeffstrs);

    return str;
}

char *TEMPLATE(T, poly_get_str) (const TEMPLATE(T, poly_t) poly,
                                 const TEMPLATE(T, ctx_t) ctx)
{
    return _TEMPLATE(T, poly_get_str) (poly->coeffs, poly->length, ctx);
}


#endif
