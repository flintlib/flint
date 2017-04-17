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

#include <ctype.h>
#include "fmpz_sparse.h"

/* returns 0 on success */
int fmpz_sparse_set_str(fmpz_sparse_t poly, const char * str)
{
    int ans = 0;
    char * w;
    slong i, len;

    if (!isdigit((unsigned char) str[0])) return -1;
    len = atol(str);
    if (len < 0) return -1;

    if (len == 0)
    {
        fmpz_sparse_zero(poly);
        return 0;
    }

    _fmpz_sparse_reserve(poly, len);
    for (i=poly->length; i<len; ++i)
    {
        fmpz_init(poly->coeffs+i);
        fmpz_init(poly->expons+i);
    }
    for (i=len+1; i<poly->length; ++i) 
    {
        fmpz_clear(poly->coeffs+i);
        fmpz_clear(poly->expons+i);
    }

    poly->length = len;

    while (*str++ != ' ') ;

    /* Find length of longest integer coeff or expon */
    {
        const char * s = str;
        slong max;
        for (max=0; *s != '\0';)
        {
            slong cur;
            for (s++, cur=1; *s != ' ' && *s != '\0'; s++, cur++) ;
            if (max < cur) max = cur;
        }

        w = flint_malloc(max+1);
    }

    for (i=0; i<poly->length; ++i)
    {
        char * v;

        for (str++, v=w; *str != ' ' && *str != '\0';)
            *v++ = *str++;
        *v = '\0';

        ans = fmpz_set_str(poly->coeffs+i, w, 10);
        if (ans) break;

        for (str++, v=w; *str != ' ' && *str != '\0';)
            *v++ = *str++;
        *v = '\0';

        ans = fmpz_set_str(poly->expons+i, w, 10);
        if (ans) break;
    }

    flint_free(w);

    if (ans)
    {
        fmpz_sparse_zero(poly);
        return -1;
    }
    else return 0;
}
