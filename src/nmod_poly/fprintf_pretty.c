/*
    Copyright (C) 2014 Jean-Pierre Flori

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"


int nmod_poly_fprint_pretty(FILE * f, const nmod_poly_t a, const char * x)
{
    size_t r;
    slong i;

    if (a->length == 0)
    {
        r = fputc('0', f);
        r = ((int) r != EOF) ? 1 : EOF;
        return r;
    }
    else if (a->length == 1)
    {
        r = flint_fprintf(f, "%wu", a->coeffs[0]);
        return r;
    }

    i = a->length - 1;
    r = 1;
    if (i == 1)
    {
        switch (a->coeffs[1])
        {
            case UWORD(0):
                break;
            case UWORD(1):
                r = flint_fprintf(f, "%s", x);
                break;
            default:
                r = flint_fprintf(f, "%wu*%s", a->coeffs[1], x);
        }
        --i;
    }
    else
    {
        switch (a->coeffs[i])
        {
            case UWORD(0):
                break;
            case UWORD(1):
                r = flint_fprintf(f, "%s^%wd", x, i);
                break;
            default:
                r = flint_fprintf(f, "%wu*%s^%wd", a->coeffs[i], x, i);
        }
        --i;
    }
    for (; (r > 0) && (i > 1); --i)
    {
        switch (a->coeffs[i])
        {
            case UWORD(0):
                break;
            case UWORD(1):
                r = flint_fprintf(f, "+%s^%wd", x, i);
                break;
            default:
                r = flint_fprintf(f, "+%wu*%s^%wd", a->coeffs[i], x, i);
        }
    }
    if (r > 0 && i == 1)
    {   
        switch (a->coeffs[1])
        {   
            case UWORD(0):
                break;
            case UWORD(1):
                r = flint_fprintf(f, "+%s", x);
                break;
            default:
                r = flint_fprintf(f, "+%wu*%s", a->coeffs[1], x);
        }
    }
    if (r > 0)
    {
        if (a->coeffs[0] != UWORD(0))
            r = flint_fprintf(f, "+%wu", a->coeffs[0]);
    }

    return (int) r;
}

