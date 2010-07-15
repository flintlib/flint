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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

int fmpz_poly_verify_format(const char * str)
{
    unsigned long i, len;
    const char * s = str;
    
    if (s == NULL || !isdigit(s[0]))
        return 0;
    
    if (s[0] == '0')
        return s[1] == '\0';
    
    for (s++; *s != '\0' && *s != ' '; s++)
        if (!(isdigit(*s)))
            return 0;
    
    if (*s++ == '\0')
        return 0;
    
    /* TODO:  Technically, we should check that the number fits in a long. */
    len = (ulong) atol(str);
    
    /* Check that s is len times " [-]####" followed by '\0'. */
    for (i = 0; i < len; i++)
    {
        if (*s++ != ' ') 
            return 0;
        if (*s == '-')
        {
            char c = *++s;
            if (!isdigit(c) || c == '0')
                return 0;
        }
        else if (isdigit(*s))
        {
            if (*s == '0')
            {
                char c = *(s + 1);
                if (c != ' ' && c != '\0')
                    return 0;
            }
        }
        else return 0;
        while (isdigit(*++s)) ;
    }
    if (*s != '\0')
        return 0;
    return 1;
}

void _fmpz_poly_from_string(fmpz * poly, const char * str)
{
    const char * s = str;
    char * v, * w;
    ulong cur, i, len, max;
    
    len = (ulong) atol(str);
    if (len == 0UL)
        return;
    
    while (s++, *str++ != ' ') ;
    
    /* Find maximal gap between spaces */
    for (max = 0; *s != '\0'; )
    {
        for (s++, cur = 1; *s != ' ' && *s != '\0'; s++, cur++) ;
        if (max < cur)
            max = cur;
    }
    
    w = (char *) malloc((max + 1) * sizeof(char));
    
    for (i = 0; i < len; i++)
    {
        for (str++, v = w; *str != ' ' && *str != '\0'; )
            *v++ = *str++;
        *v = '\0';
        fmpz_set_str(poly++, w, 10);
    }
    
    free(w);
}

int fmpz_poly_from_string(fmpz_poly_t poly, const char * str)
{
    int check = fmpz_poly_verify_format(str);
    ulong len;
    
    if (!check)
        return 0;
    
    len = (ulong) atol(str);
    fmpz_poly_fit_length(poly, len);
    
    _fmpz_poly_from_string(poly->coeffs, str);
    
    _fmpz_poly_set_length(poly, len);
    _fmpz_poly_normalise(poly);
    return 1;
}

