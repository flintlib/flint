/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2007, David Howden.
    Copyright (C) 2010 William Hart

******************************************************************************/

#include <gmp.h>
#include <stdio.h>
#include <string.h>
#include "flint.h"
#include "nmod_poly.h"

int nmod_poly_set_str(nmod_poly_t poly, const char * s)
{
    const char * whitespace = " \t\n\r";
    len_t i, length;
    mp_limb_t n;

    if (sscanf(s, "%ld %lu", &length, &n) != 2)
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
      
        if (!sscanf(s, "%lu", &poly->coeffs[i]))
        {
            poly->length = i;
            return 0;
        }
    }
   
    _nmod_poly_normalise(poly);
   
    return 1;
}
