/*
    Copyright (C) 2007, David Howden.
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "flint.h"
#include "nmod_poly.h"

static char * _nmod_poly_get_str_pretty(const mp_limb_t * Acoeff, slong Alen,
                                                                const char * x)
{
    slong i;
    char * buf, * ptr;

    slong size = 0;

    if (Alen == 0)
    {
        buf = (char *) flint_malloc(2);
        buf[0] = '0';
        buf[1] = '\0';
        return buf;
    }
    else if (Alen == 1)
    {
        size = (ulong) ceil(0.30103*FLINT_BIT_COUNT(Acoeff[0])) + 1;
        buf = (char *) flint_malloc(size);
        flint_sprintf(buf, "%wu", Acoeff[0]);
        return buf;
    }

    for (i = 0; i < Alen; i++)
    {
        if (Acoeff[i]) /* log(2)/log(10) < 0.30103, +3 for +*^ or null*/
            size += (ulong) ceil(0.30103*FLINT_BIT_COUNT(Acoeff[i])) + 
                    (ulong) ceil(0.30103*FLINT_BIT_COUNT(i)) + strlen(x) + 3; 
    }

    buf = (char *) flint_malloc(size);  
    ptr = buf;
    --i;
    if (i == 1)
    {
        switch (Acoeff[1])
        {
            case UWORD(1):
                ptr += flint_sprintf(ptr, "%s", x);
                break;
            default:
                ptr += flint_sprintf(ptr, "%wu*%s", Acoeff[1], x);
        }
        --i;
    }
    else
    {
        switch (Acoeff[i])
        {
            case UWORD(1):
                ptr += flint_sprintf(ptr, "%s^%wd", x, i);
                break;
            default:
                ptr += flint_sprintf(ptr, "%wu*%s^%wd", Acoeff[i], x, i);
        }
        --i;
    }
    for (; i > 1; --i)
    {
        switch (Acoeff[i])
        {
            case UWORD(0):
                break;
            case UWORD(1):
                ptr += flint_sprintf(ptr, "+%s^%wd", x, i);
                break;
            default:
                ptr += flint_sprintf(ptr, "+%wu*%s^%wd", Acoeff[i], x, i);
        }
        
    }
    if (i == 1)
    {   
        switch (Acoeff[1])
        {   
            case UWORD(0):
                break;
            case UWORD(1):
                ptr += flint_sprintf(ptr, "+%s", x);
                break;
            default:
                ptr += flint_sprintf(ptr, "+%wu*%s", Acoeff[1], x);
        }
    }
    {
        if (Acoeff[0] != UWORD(0))
            ptr += flint_sprintf(ptr, "+%wu", Acoeff[0]);
    }

    return buf;
}

char * nmod_polydr_get_str_pretty(const nmod_polydr_t A, const char * x,
                                                          const nmod_ctx_t ctx)
{
    return _nmod_poly_get_str_pretty(A->coeffs, A->length, x);
}

char * nmod_poly_get_str_pretty(const nmod_poly_t A, const char * x)
{
    return _nmod_poly_get_str_pretty(A->coeffs, A->length, x);
}
