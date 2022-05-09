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
#include "flint-impl.h"
#include "nmod_poly_mini.h"

int nmod_poly_set_str(nmod_poly_t poly, const char * s)
{
    const char * whitespace = " \t\n\r";
    slong i, length;
    ulong n;

    if (sscanf(s, WORD_FMT "d " WORD_FMT "u", &length, &n) != 2)
        return 0;
      
    /* jump past length (n will be skipped in first loop iter)  */
    s += strcspn(s, whitespace);
    s += strspn(s, whitespace);
    
    nmod_poly_fit_length(poly, length);
    poly->length = length;
    
    for (i = 0; i < length; i++)
    {
        s += strcspn(s, whitespace); /* jump to next whitespace */
        s += strspn(s, whitespace); /* skip whitespace */
      
        if (!sscanf(s, WORD_FMT "u", &poly->coeffs[i]))
        {
            poly->length = i;
            return 0;
        }
    }
   
    _nmod_poly_normalise(poly);
   
    return 1;
}
