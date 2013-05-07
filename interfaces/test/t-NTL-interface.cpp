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

    Copyright (C) 2011 William Hart

******************************************************************************/

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"
#include "NTL-interface.h"

NTL_CLIENT

int test_ZZ_to_fmpz()
{
    int i, result;
    flint_rand_t state;
    mp_bitcnt_t bits, randbits;
    fmpz_t int1, int2;
   
    ZZ z;

    printf("ZZ_to_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    /* Check conversion */
    for (i = 0; i < 10000; i++)
    {
        bits = n_randint(state, 1000) + 1;
        randbits = n_randint(state, bits);
      
        fmpz_init(int1);
        fmpz_init(int2);
        
        fmpz_randbits(int1, state, randbits);

        fmpz_get_ZZ(z, int1);
        fmpz_set_ZZ(int2, z);
      
        result = fmpz_equal(int1, int2);
        if (!result)
        {
           printf("FAIL:\n");
           printf("int1 = %ld  ", *int1); fmpz_print(int1); printf("\n");
           printf("int2 = %ld  ", *int2); fmpz_print(int2); printf("\n");
           return 0;
        }

        fmpz_clear(int1);
        fmpz_clear(int2);
    }    

    return 1;
}

int test_ZZX_to_fmpz_poly()
{
    fmpz_poly_t f_poly1, f_poly2;
    long length;
    mp_bitcnt_t bits;
    flint_rand_t state;
    int i, result;
   
    printf("ZZX_to_fmpz_poly....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and c */
    for (i = 0; i < 10000; i++)
    {
        bits = n_randint(state, 1000) + 1;
        length = n_randint(state, 1000);

        fmpz_poly_init(f_poly1);
        fmpz_poly_init(f_poly2);
        
        ZZX ZZX_poly;
          
        fmpz_poly_randtest(f_poly1, state, length, bits);
          
        fmpz_poly_get_ZZX(ZZX_poly, f_poly1);
        fmpz_poly_set_ZZX(f_poly2, ZZX_poly);
          
        result = fmpz_poly_equal(f_poly1, f_poly2);  
        if (!result)
        {
           printf("FAIL:\n");
           printf("f_poly1 = "); fmpz_poly_print(f_poly1); printf("\n");
           printf("f_poly2 = "); fmpz_poly_print(f_poly2); printf("\n");
           return 0;
        }
          
        fmpz_poly_clear(f_poly1);
        fmpz_poly_clear(f_poly2);
    }
      
    return 1;
}

int
main(void)
{
    int r1, r2;

    if ((r1 = test_ZZ_to_fmpz())) printf("PASS\n");
    if ((r2 = test_ZZX_to_fmpz_poly())) printf("PASS\n");

    if (!r1 || !r2) abort();

    _fmpz_cleanup();
    
    return 0;
}
