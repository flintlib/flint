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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "flint.h"
#include "nmod_poly.h"

char * nmod_poly_get_str(const nmod_poly_t poly)
{
    long i;
    char * buf, * ptr;

    /* estimate for the length, n and three spaces */
#if FLINT64
    long size = 21*2 + 1;
#else
    long size = 11*2 + 1;
#endif

    for (i = 0; i < poly->length; i++)
    {
        if (poly->coeffs[i]) /* log(2)/log(10) < 0.30103, +1 for space/null */
            size += (ulong) ceil(0.30103*FLINT_BIT_COUNT(poly->coeffs[i])) + 1;
        else size += 2;
    }

    buf = (char *) flint_malloc(size);  
    ptr = buf + sprintf(buf, "%ld %lu", poly->length, poly->mod.n);
   
    if (poly->length)
        ptr += sprintf(ptr, " ");

    for (i = 0; i < poly->length; i++)
        ptr += sprintf(ptr, " %lu", poly->coeffs[i]);
   
    return buf;
}

