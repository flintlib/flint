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
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_pow_addchains(fmpz * res, const fmpz * poly, len_t len, 
                                                          const int * a, int n)
{
    int *b;
    len_t lenm1 = len - 1, lenv;
    fmpz *v;

    /*
       Compute partial sums
     */
    {
        int i;
        b = (int *) flint_malloc(n * sizeof(int));
        b[0] = 0;
        for (i = 1; i < n; i++)
            b[i] = b[i-1] + a[i];
    }
    
    /*
       Allocate memory for the polynomials f^{a[1]}, ..., f^{a[n-1]}
     */
    lenv = lenm1 * b[n-1] + n - 1;
    v = _fmpz_vec_init(lenv);
    
    /* 
       Compute f^{a[1]}, ..., f^{a[n-1]}
     */
    {
        int d, i, j;
        
        _fmpz_poly_sqr(v, poly, len);
        
        for (i = 1; i < n-1; i++)
        {
            d = a[i+1] - a[i];
            if (d == 1)
            {
                _fmpz_poly_mul(v + lenm1 * b[i] + (i),
                               v + lenm1 * b[i-1], lenm1 * a[i] + 1,
                               poly, len);
            }
            else
            {
                for (j = i; a[j] != d; j--) ;
                _fmpz_poly_mul(v + lenm1 * b[i] + (i),
                               v + lenm1 * b[i-1], lenm1 * a[i] + 1,
                               v + lenm1 * b[j-1] + (j-1), lenm1 * a[j] + 1);
            }
        }

        /*
           Deal with the final product stored in res, i == n-1
         */
        {
            d = a[i+1] - a[i];
            if (d == 1)
            {
                _fmpz_poly_mul(res,
                               v + lenm1 * b[i-1], lenm1 * a[i] + 1,
                               poly, len);
            }
            else
            {
                for (j = i; a[j] != d; j--) ;
                _fmpz_poly_mul(res,
                               v + lenm1 * b[i-1], lenm1 * a[i] + 1,
                               v + lenm1 * b[j-1] + (j-1), lenm1 * a[j] + 1);
            }
        }
        
    }

    flint_free(b);
    _fmpz_vec_clear(v, lenv);
}

void fmpz_poly_pow_addchains(fmpz_poly_t res, const fmpz_poly_t poly, ulong e)
{
    const len_t len = poly->length;
    
    if ((len < 2) | (e < 3UL))
    {
        if (e == 0UL)
            fmpz_poly_set_ui(res, 1);
        else if (len == 0)
            fmpz_poly_zero(res);
        else if (len == 1)
        {
            fmpz_poly_fit_length(res, 1);
            fmpz_pow_ui(res->coeffs, poly->coeffs, e);
            _fmpz_poly_set_length(res, 1);
        }
        else if (e == 1UL)
            fmpz_poly_set(res, poly);
        else  /* e == 2UL */
            fmpz_poly_sqr(res, poly);
        return;
    }
    
    if (e <= 148UL)
    {
        /*
           An array storing a tree with shortest addition chains (star chains, 
           in fact) for all integers up to and including 148.  
           
           Let A denote the array.  The entry A[0] is present to provide 
           1-based indexing.  The integer 1 is the root of the tree and the 
           entry A[1] is irrelevant.  For integers i >= 2, A[i] is the parent 
           of i.
           
           We can iterate through an addition chain for n, where 0 < n < 148, 
           in the array shortest_addchains_148 as follows:
           
               Visit n
               while ((n = shortest_addchains_148[n]))
               {
                   Visit n
               }
         */
        static const int shortest_addchains_148[149] = {
             0,  0,  1,  2,  2,  3,  3,  5,  4,  8,
             5, 10,  6,  9,  7, 12,  8,  9, 16, 18,
            10, 15, 11, 20, 12, 17, 13, 24, 14, 25, 
            15, 28, 16, 32, 17, 26, 18, 36, 19, 27,
            20, 40, 21, 34, 22, 30, 23, 46, 24, 33,
            25, 48, 26, 37, 27, 54, 28, 49, 29, 56, 
            30, 52, 31, 51, 32, 64, 33, 66, 34, 68,
            35, 70, 36, 72, 66, 60, 38, 43, 39, 78,
            40, 65, 41, 80, 42, 80, 43, 86, 44, 88,
            45, 90, 46, 92, 47, 92, 48, 96, 49, 96, 
            50,100, 51,102, 52,102, 53, 74, 54,108,
            55,108, 56,104, 57,112, 58,104, 59,112,
            60,120, 61,120, 62,100, 63,126, 64,128,
            65,130,128,132, 67, 90, 68,136, 69,138,
            70,140, 71,117, 72,144, 73, 99, 74 
        };
        
        int a[11], i = 11, n = (int) e;
        len_t rlen = (len_t) e * (len - 1) + 1;

        /*
           Copy the addition chain into 1 = a[0] < a[1] < ... < a[n]
         */
        a[--i] = n;
        while ((n = shortest_addchains_148[n]))
            a[--i] = n;
        n = 10 - i;

        if (res != poly)
        {
            fmpz_poly_fit_length(res, rlen);
            _fmpz_poly_pow_addchains(res->coeffs, poly->coeffs, len, a + i, n);
            _fmpz_poly_set_length(res, rlen);
        }
        else
        {
            fmpz_poly_t t;
            fmpz_poly_init2(t, rlen);
            _fmpz_poly_pow_addchains(t->coeffs, poly->coeffs, len, a + i, n);
            _fmpz_poly_set_length(t, rlen);
            fmpz_poly_swap(res, t);
            fmpz_poly_clear(t);
        }
    }
    else
    {
        printf("Exception (fmpz_poly_addchains). Powering via chains not implemented for e > 148.\n");
        abort();
    }
}
