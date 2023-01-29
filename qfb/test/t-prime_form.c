/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "qfb.h"

int main(void)
{
    int result;
    flint_rand_t state;
    slong i;
    
    printf("prime_form....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 1; i < 50000; i++) 
    {
        fmpz_t D, D2, p;
        qfb_t r, s;
        ulong n, q, Dmodp;
        
        fmpz_init(D);
        fmpz_init(D2);
        fmpz_init(p);
        qfb_init(r);
        qfb_init(s);
            
        do
        {
           fmpz_randtest_unsigned(r->a, state, 100);
           if (fmpz_is_zero(r->a))
              fmpz_set_ui(r->a, 1);
 
           fmpz_randtest(r->b, state, 100);
           fmpz_randtest(r->c, state, 100);

           qfb_discriminant(D, r);
        } while (fmpz_sgn(D) >= 0);

        do
        {
           n = n_randprime(state, n_randint(state, FLINT_BITS - 1) + 2, 0);
           fmpz_set_ui(p, n);
           Dmodp = fmpz_fdiv_ui(D, n);
        } while ((mp_limb_signed_t) Dmodp < 0 /* Jacobi can't handle this */
         || (n == 2 && ((q = fmpz_fdiv_ui(D, 8)) == 2 || q == 3 || q == 5))
         || (n != 2 && Dmodp != 0 && n_jacobi(Dmodp, n) < 0));

        qfb_prime_form(s, D, p);
        qfb_discriminant(D2, s);

        result = (fmpz_equal(s->a, p) && fmpz_equal(D2, D) && fmpz_sgn(s->b) >= 0);
        if (!result)
        {
           printf("FAIL:\n");
           qfb_print(s); printf("\n");
           printf("p = \n"); fmpz_print(p); printf("\n");
           printf("D = \n"); fmpz_print(D); printf("\n");
           printf("D2 = \n"); fmpz_print(D2); printf("\n");
           abort();
        }
           
        fmpz_clear(D);
        fmpz_clear(D2);
        fmpz_clear(p);
        qfb_clear(r);
        qfb_clear(s);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
