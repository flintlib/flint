/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <stdio.h>
#include <gmp.h>
#include "fmpz.h"
int
_TEMPLATE(T, poly_fprint) (FILE * file, const TEMPLATE(T, struct) * poly,
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
TEMPLATE(T, poly_fprint) (FILE * file, const TEMPLATE(T, poly_t) poly,
                          const TEMPLATE(T, ctx_t) ctx)
{
    return _TEMPLATE(T, poly_fprint) (file, poly->coeffs, poly->length, ctx);
}


#endif
