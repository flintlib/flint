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
#include <stdio.h>
#include "flint.h"
#include "nmod_poly.h"

int nmod_polydr_fread(FILE * f, nmod_polydr_t poly, const nmod_ctx_t ctx)
{
    slong i, length;
    mp_limb_t n;

    if (flint_fscanf(f, "%wd %wu", &length, &n) != 2)
        return 0;

    if (n != ctx->mod.n)
        return 0;
    
    nmod_polydr_clear(poly, ctx);
    nmod_polydr_init(poly, ctx); 
    nmod_polydr_fit_length(poly, length, ctx);
    poly->length = length;
    
    for (i = 0; i < length; i++)
    {
        if (!flint_fscanf(f, "%wu", &poly->coeffs[i]))
        {
            poly->length = i;
            return 0;
        }
    }
   
    _nmod_polydr_normalise(poly);
   
    return 1;
}

int nmod_poly_fread(FILE * f, nmod_poly_t poly)
{
    slong i, length;
    mp_limb_t n;

    if (flint_fscanf(f, "%wd %wu", &length, &n) != 2)
        return 0;
    
    nmod_poly_clear(poly);
    nmod_poly_init(poly,n); 
    nmod_poly_fit_length(poly, length);
    poly->length = length;
    
    for (i = 0; i < length; i++)
    {
        if (!flint_fscanf(f, "%wu", &poly->coeffs[i]))
        {
            poly->length = i;
            return 0;
        }
    }
   
    _nmod_poly_normalise(poly);
   
    return 1;
}
