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
#include "fmpz_poly.h"

int
_fmpz_poly_set_str(fmpz * poly, const char *str)
{
    char * w;
    len_t i, len;

    if (!isdigit(str[0]))
        return -1;
    len = atol(str);
    if (len < 0)
        return -1;
    if (len == 0)
        return 0;

    while (*str++ != ' ')
        ;

    /* Find maximal gap between spaces and allocate w */
    {
        const char * s = str;
        len_t max;
        for (max = 0; *s != '\0';)
        {
            len_t cur;
            for (s++, cur = 1; *s != ' ' && *s != '\0'; s++, cur++) ;
            if (max < cur)
                max = cur;
        }

        w = flint_malloc(max + 1);
    }

    for (i = 0; i < len; i++)
    {
        char * v;
        int ans;

        for (str++, v = w; *str != ' ' && *str != '\0';)
            *v++ = *str++;
        *v = '\0';
        ans = fmpz_set_str(poly++, w, 10);

        if (ans)
        {
            flint_free(w);
            return -1;
        }
    }

    flint_free(w);
    return 0;
}

int
fmpz_poly_set_str(fmpz_poly_t poly, const char * str)
{
    int ans;
    len_t len;

    if (!isdigit(str[0]))
        return -1;
    len = atol(str);
    if (len < 0)
        return -1;
    if (len == 0)
    {
        fmpz_poly_zero(poly);
        return 0;
    }

    fmpz_poly_fit_length(poly, len);

    ans = _fmpz_poly_set_str(poly->coeffs, str);

    if (ans == 0)
    {
        _fmpz_poly_set_length(poly, len);
        _fmpz_poly_normalise(poly);
    }
    else
    {
        _fmpz_vec_zero(poly->coeffs, len);
        _fmpz_poly_set_length(poly, 0);
    }

    return ans;
}
