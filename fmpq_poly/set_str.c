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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

int
_fmpq_poly_set_str(fmpz * poly, fmpz_t den, const char * str)
{
    char * w;
    long i, len;
    mpq_t * a;

    len = atol(str);
    if (len < 0)
        return -1;
    if (len == 0)
    {
        fmpz_one(den);
        return 0;
    }

    a = (mpq_t *) flint_malloc(len * sizeof(mpq_t));

    while (*str++ != ' ')
        ;

    /* Find maximal gap between spaces and allocate w */
    {
        const char * s = str;
        long max;
        for (max = 0; *s != '\0';)
        {
            long cur;
            for (s++, cur = 1; *s != ' ' && *s != '\0'; s++, cur++) ;
            if (max < cur)
                max = cur;
        }

        w = (char *) flint_malloc((max + 1) * sizeof(char));
    }

    for (i = 0; i < len; i++)
    {
        char * v;
        int ans;

        for (str++, v = w; *str != ' ' && *str != '\0';)
            *v++ = *str++;
        *v = '\0';
        mpq_init(a[i]);
        ans = mpq_set_str(a[i], w, 10);

        /* If the format is not correct, clear up and return -1 */
        if (ans)
        {
            int j;
            for (j = 0; j <= i; j++)
                mpq_clear(a[j]);
            flint_free(a);
            flint_free(w);
            return -1;
        }
    }

    _fmpq_poly_set_array_mpq(poly, den, (const mpq_t *) a, len);

    for (i = 0; i < len; i++)
        mpq_clear(a[i]);
    flint_free(a);
    flint_free(w);

    return 0;
}

int
fmpq_poly_set_str(fmpq_poly_t poly, const char * str)
{
    int ans;
    long len;

    len = atol(str);
    if (len < 0)
        return -1;
    if (len == 0)
    {
        fmpq_poly_zero(poly);
        return 0;
    }

    fmpq_poly_fit_length(poly, len);
    
    ans = _fmpq_poly_set_str(poly->coeffs, poly->den, str);
    
    if (ans == 0)
    {
        _fmpq_poly_set_length(poly, len);
        _fmpq_poly_normalise(poly);
    }
    else
    {
        _fmpz_vec_zero(poly->coeffs, len);
        fmpz_one(poly->den);
        _fmpq_poly_set_length(poly, 0);
    }
    
    return ans;
}
