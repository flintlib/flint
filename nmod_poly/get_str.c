/*
    Copyright (C) 2007, David Howden.
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "flint.h"
#include "nmod_poly.h"

static char * _nmod_poly_get_str(const mp_limb_t * Acoeff, slong Alen, nmod_t mod)
{
    slong i;
    char * buf, * ptr;

    /* estimate for the length, n and three spaces */
#if FLINT64
    slong size = 21*2 + 1;
#else
    slong size = 11*2 + 1;
#endif

    for (i = 0; i < Alen; i++)
    {
        if (Acoeff[i]) /* log(2)/log(10) < 0.30103, +1 for space/null */
            size += (ulong) ceil(0.30103*FLINT_BIT_COUNT(Acoeff[i])) + 1;
        else size += 2;
    }

    buf = (char *) flint_malloc(size);  
    ptr = buf + flint_sprintf(buf, "%wd %wu", Alen, mod.n);
   
    if (Alen)
        ptr += flint_sprintf(ptr, " ");

    for (i = 0; i < Alen; i++)
        ptr += flint_sprintf(ptr, " %wu", Acoeff[i]);
   
    return buf;
}

char * nmod_polydr_get_str(const nmod_polydr_t poly, const nmod_ctx_t ctx)
{
    return _nmod_poly_get_str(poly->coeffs, poly->length, ctx->mod);
}

char * nmod_poly_get_str(const nmod_poly_t poly)
{
    return _nmod_poly_get_str(poly->coeffs, poly->length, poly->mod);
}
