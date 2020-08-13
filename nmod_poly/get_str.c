/*
    Copyright (C) 2007, David Howden.
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "flint.h"
#include "nmod_poly.h"

char * nmod_poly_get_str(const nmod_poly_t poly)
{
    slong i;
    char * buf, * ptr;

    /* estimate for the length, n and three spaces */
#if FLINT64
    slong size = 21*2 + 1;
#else
    slong size = 11*2 + 1;
#endif

    for (i = 0; i < poly->length; i++)
    {
        if (poly->coeffs[i]) /* log(2)/log(10) < 0.30103, +1 for space/null */
            size += (ulong) ceil(0.30103*FLINT_BIT_COUNT(poly->coeffs[i])) + 1;
        else size += 2;
    }

    buf = (char *) flint_malloc(size);  
    ptr = buf + flint_sprintf(buf, "%wd %wu", poly->length, poly->mod.n);
   
    if (poly->length)
        ptr += flint_sprintf(ptr, " ");

    for (i = 0; i < poly->length; i++)
        ptr += flint_sprintf(ptr, " %wu", poly->coeffs[i]);
   
    return buf;
}

