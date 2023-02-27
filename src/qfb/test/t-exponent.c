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
    qfb * forms;
    slong i, k, i1, d, num;

    printf("exponent....");
    fflush(stdout);

    flint_randinit(state);

    /* Check exponent divides class number */
    for (i = 1; i < 1000; i++) 
    {
        fmpz_t D, exp;
        qfb_t pow;
        slong e;
        
        d = n_randint(state, 100000) + 1;
        num = qfb_reduced_forms(&forms, -d);
        
        if (num)
        {
           fmpz_init(D);
           fmpz_init(exp);
           qfb_init(pow);

           fmpz_set_si(D, -d);
           
           result = qfb_exponent(exp, D, 1000000, 100000, 20);
           if (!result)
           {
              printf("FAIL:\n");
              printf("Exponent not computed\n");
              printf("Discriminant: "); fmpz_print(D); printf("\n");
              abort();
           }
           
           e = fmpz_get_si(exp);
           result = ((num % e) == 0);
           if (!result)
           {
              printf("FAIL:\n");
              printf("Exponent does not divide class number\n");
              printf("Discriminant: "); fmpz_print(D); printf("\n");
              printf("%ld does not divide %ld\n", e, num);
              abort();
           }
           
           for (k = 0; k < 5; k++)
           {
              i1 = n_randint(state, num);
              qfb_pow_ui(pow, forms + i1, D, e);
              result = (qfb_is_principal_form(pow, D));
              if (!result)
              {
                 printf("FAIL:\n");
                 printf("Exponent does not annihilate form\n");
                 printf("Discriminant: "); fmpz_print(D); printf("\n");
                 printf("Form: "); qfb_print(forms + i1); printf("\n");
                 printf("Exponent: %ld\n", e);
                 abort();
              }
           }

           qfb_clear(pow);
           fmpz_clear(D);
           fmpz_clear(exp);
        }

        qfb_array_clear(&forms, num);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
