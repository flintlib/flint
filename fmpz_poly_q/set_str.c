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

    Copyright (C) 2010, 2011 Sebastian Pancratz
   
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "fmpz_poly_q.h"

/**
 * \ingroup  StringConversions
 * 
 * Sets the rational function \c rop to the value specified by the 
 * null-terminated string \c s.
 *
 * This method has now already been somewhat improved and is not very tolerant 
 * in the handling of malformed input.  It expects either legitimate input for 
 * an \c fmpz_poly_t element, or two such inputs separated by a <tt>/</tt> 
 * only, in which case it is also assumed that the second polynomial is 
 * non-zero.
 * 
 * The rational function is brought into canonical form by calling 
 * #fmpz_poly_q_canonicalize() in this function.
 * 
 * Returns \c 0 if the string represents a valid rational function and 
 * \c non-zero otherwise.
 */
int fmpz_poly_q_set_str(fmpz_poly_q_t rop, const char *s)
{
    int ans, i, m;
    size_t len;
    char * numstr;
    
    len = strlen(s);
    
    for (m = 0; m < len; m++)
    {
        if (s[m] == '/')
            break;
    }
    
    if (m == len)
    {
        ans = fmpz_poly_set_str(rop->num, s);
        fmpz_poly_set_si(rop->den, 1);
        return ans;
    }
    else
    {
        numstr = flint_malloc(m + 1);
        if (!numstr)
        {
            flint_printf("Exception (fmpz_poly_q_set_str). Memory allocation failed.\n");
            flint_abort();
        }
        
        for (i = 0; i < m; i++)
            numstr[i] = s[i];
        numstr[i] = '\0';
        
        ans  = fmpz_poly_set_str(rop->num, numstr);
        ans |= fmpz_poly_set_str(rop->den, s + (m + 1));
        if (ans == 0)
            fmpz_poly_q_canonicalise(rop);
        else
            fmpz_poly_q_zero(rop);
        flint_free(numstr);
        return ans;
    }
}
