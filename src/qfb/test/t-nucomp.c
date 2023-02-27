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
    slong i, k, i1, i2, i3, num;

    printf("nucomp....");
    fflush(stdout);

    flint_randinit(state);

    /* Check discriminants of composed forms are correct */
    for (i = 1; i < 10000; i++) 
    {
        qfb_t r;
        fmpz_t root, D;

        num = qfb_reduced_forms(&forms, -i);
        
        if (num)
        {
           fmpz_init(root);
           fmpz_init(D);
           qfb_init(r);
              
           fmpz_set_ui(root, i);
           fmpz_root(root, root, 4);

           for (k = 0; k < 20; k++)
           {
              i1 = n_randint(state, num);
              i2 = n_randint(state, num);

              fmpz_set_si(D, -i);
              qfb_nucomp(r, forms + i1, forms + i2, D, root);

              qfb_discriminant(D, r);
              result = (fmpz_cmp_si(D, -i) == 0);
              if (!result)
              {
                 printf("FAIL:\n");
                 printf("Incorrect discriminant\n");
                 fmpz_print(D); printf(" should be %ld\n", -i);
                 abort();
              }
           }
           
           fmpz_clear(D);
           fmpz_clear(root);
           qfb_clear(r);
        }

        qfb_array_clear(&forms, num);
    }

    /* Associativity test (Q1 o Q2) o Q3 == Q1 o (Q2 o Q3) */
    for (i = 1; i < 10000; i++) 
    {
        qfb_t r, s;
        fmpz_t root, t1, t2, D;

        num = qfb_reduced_forms(&forms, -i);
        
        if (num)
        {
           fmpz_init(root);
           fmpz_init(t1);
           fmpz_init(t2);
           fmpz_init(D);
           qfb_init(r);
           qfb_init(s);
              
           fmpz_set_si(D, -i);
           fmpz_set_ui(root, i);
           fmpz_root(root, root, 4);

           for (k = 0; k < 20; k++)
           {
              i1 = n_randint(state, num);
              i2 = n_randint(state, num);
              i3 = n_randint(state, num);

              qfb_nucomp(r, forms + i1, forms + i2, D, root);
              qfb_reduce(r, r, D);

              qfb_nucomp(r, r, forms + i3, D, root);
              qfb_reduce(r, r, D);

              qfb_nucomp(s, forms + i2, forms + i3, D, root);
              qfb_reduce(s, s, D);

              qfb_nucomp(s, forms + i1, s, D, root);
              qfb_reduce(s, s, D);

              result = (qfb_equal(r, s));
              if (!result)
              {
                 printf("FAIL:\n");
                 printf("Associativity fails\n");
                 printf("r = "); qfb_print(r); printf("\n");
                 printf("s = "); qfb_print(s); printf("\n");
                 abort();
              }
           }
           
           fmpz_clear(t1);
           fmpz_clear(t2);
           fmpz_clear(D);
           fmpz_clear(root);
           qfb_clear(r);
           qfb_clear(s);
        }

        qfb_array_clear(&forms, num);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
