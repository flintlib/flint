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

    Authored 2015 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include <math.h>
#include <gmp.h>
#include "fmpz_sparse.h"

/* note: caller's responsibility to deallocate */
char * fmpz_sparse_get_str(const fmpz_sparse_t poly)
{
    slong i, bound;
    char *str, *strbase;
    const fmpz *cptr, *eptr;

    if (poly->length == 0)
    {
        str = (char *) flint_malloc(2 * sizeof(char));
        str[0] = '0';
        str[1] = '\0';
        return str;
    }

    bound = (slong) (ceil(log10((double) (poly->length + 1)))) + 3;
    for (i=0, cptr=poly->coeffs, eptr=poly->expons; i<poly->length; 
         ++i, ++cptr, ++eptr)
    {
        bound += fmpz_sizeinbase(cptr, 10) + 2;
        bound += fmpz_sizeinbase(eptr, 10) + 2;
    }

    strbase = (char*) flint_malloc(bound * sizeof(char));
    str = strbase;

    str += flint_sprintf(str, "%wd ", poly->length);
    for (i=0, cptr=poly->coeffs, eptr=poly->expons; i<poly->length; 
         ++i, ++cptr, ++eptr)
    {
        if (!COEFF_IS_MPZ(*cptr))
            str += flint_sprintf(str, " %wd", *cptr);
        else
            str += gmp_sprintf(str, " %Zd", COEFF_TO_PTR(*cptr));
        if (!COEFF_IS_MPZ(*eptr))
            str += flint_sprintf(str, " %wd", *eptr);
        else
            str += gmp_sprintf(str, " %Zd", COEFF_TO_PTR(*eptr));
    }

    FLINT_ASSERT ((strbase + bound) <= str);
    return strbase;
}
