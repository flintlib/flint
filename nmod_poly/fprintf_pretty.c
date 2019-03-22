/*
    Copyright (C) 2014 Jean-Pierre Flori

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"


static int _nmod_poly_fprint_pretty(FILE * f,
                                mp_limb_t * Acoeff, slong Alen, const char * x)
{
    size_t r;
    slong i;

    if (Alen == 0)
    {
        r = fputc('0', f);
        r = ((int) r != EOF) ? 1 : EOF;
        return r;
    }
    else if (Alen == 1)
    {
        r = flint_fprintf(f, "%wu", Acoeff[0]);
        return r;
    }

    i = Alen - 1;
    r = 1;
    if (i == 1)
    {
        switch (Acoeff[1])
        {
            case UWORD(0):
                break;
            case UWORD(1):
                r = flint_fprintf(f, "%s", x);
                break;
            default:
                r = flint_fprintf(f, "%wu*%s", Acoeff[1], x);
        }
        --i;
    }
    else
    {
        switch (Acoeff[i])
        {
            case UWORD(0):
                break;
            case UWORD(1):
                r = flint_fprintf(f, "%s^%wd", x, i);
                break;
            default:
                r = flint_fprintf(f, "%wu*%s^%wd", Acoeff[i], x, i);
        }
        --i;
    }
    for (; (r > 0) && (i > 1); --i)
    {
        switch (Acoeff[i])
        {
            case UWORD(0):
                break;
            case UWORD(1):
                r = flint_fprintf(f, "+%s^%wd", x, i);
                break;
            default:
                r = flint_fprintf(f, "+%wu*%s^%wd", Acoeff[i], x, i);
        }
    }
    if (r > 0 && i == 1)
    {   
        switch (Acoeff[1])
        {   
            case UWORD(0):
                break;
            case UWORD(1):
                r = flint_fprintf(f, "+%s", x);
                break;
            default:
                r = flint_fprintf(f, "+%wu*%s", Acoeff[1], x);
        }
    }
    if (r > 0)
    {
        if (Acoeff[0] != UWORD(0))
            r = flint_fprintf(f, "+%wu", Acoeff[0]);
    }

    return (int) r;
}


int nmod_polydr_fprint_pretty(FILE * f, const nmod_polydr_t A, const char * x,
                                                          const nmod_ctx_t ctx)
{
    return _nmod_poly_fprint_pretty(f, A->coeffs, A->length, x);
}

int nmod_poly_fprint_pretty(FILE * f, const nmod_poly_t A, const char * x)
{
    return _nmod_poly_fprint_pretty(f, A->coeffs, A->length, x);
}
