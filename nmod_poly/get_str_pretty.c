/*
    Copyright (C) 2007, David Howden.
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "flint.h"
#include "flint-impl.h"

char * nmod_poly_get_str_pretty(const nmod_poly_t poly, const char * x)
{
    slong i;
    char * buf, * ptr;

    slong size = 0;

    if (poly->length == 0)
    {
        buf = (char *) flint_malloc(2);
        buf[0] = '0';
        buf[1] = '\0';
        return buf;
    }
    else if (poly->length == 1)
    {
        size = (ulong) ceil(0.30103*FLINT_BIT_COUNT(poly->coeffs[0])) + 1;
        buf = (char *) flint_malloc(size);
        sprintf(buf, WORD_FMT "u", poly->coeffs[0]);
        return buf;
    }

    for (i = 0; i < poly->length; i++)
    {
        if (poly->coeffs[i]) /* log(2)/log(10) < 0.30103, +3 for +*^ or null*/
            size += (ulong) ceil(0.30103*FLINT_BIT_COUNT(poly->coeffs[i])) + 
                    (ulong) ceil(0.30103*FLINT_BIT_COUNT(i)) + strlen(x) + 3; 
    }

    buf = (char *) flint_malloc(size);  
    ptr = buf;
    --i;
    if (i == 1)
    {
        switch (poly->coeffs[1])
        {
            case UWORD(1):
                ptr += sprintf(ptr, "%s", x);
                break;
            default:
                ptr += sprintf(ptr, WORD_FMT "u*%s", poly->coeffs[1], x);
        }
        --i;
    }
    else
    {
        switch (poly->coeffs[i])
        {
            case UWORD(1):
                ptr += sprintf(ptr, "%s^" WORD_FMT "d", x, i);
                break;
            default:
                ptr += sprintf(ptr, WORD_FMT "u*%s^" WORD_FMT "d", poly->coeffs[i], x, i);
        }
        --i;
    }
    for (; i > 1; --i)
    {
        switch (poly->coeffs[i])
        {
            case UWORD(0):
                break;
            case UWORD(1):
                ptr += sprintf(ptr, "+%s^" WORD_FMT "d", x, i);
                break;
            default:
                ptr += sprintf(ptr, "+" WORD_FMT "u*%s^" WORD_FMT "d", poly->coeffs[i], x, i);
        }
        
    }
    if (i == 1)
    {   
        switch (poly->coeffs[1])
        {   
            case UWORD(0):
                break;
            case UWORD(1):
                ptr += sprintf(ptr, "+%s", x);
                break;
            default:
                ptr += sprintf(ptr, "+" WORD_FMT "u*%s", poly->coeffs[1], x);
        }
    }
    {
        if (poly->coeffs[0] != UWORD(0))
            ptr += sprintf(ptr, "+" WORD_FMT "u", poly->coeffs[0]);
    }

    return buf;
}

