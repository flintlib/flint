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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qfb.h"

int main(void)
{
    int result;
    flint_rand_t state;
    qfb * forms;
    long i, k, i1, i2, i3, num;

    printf("nucomp....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 1; i < 10000; i++) 
    {
        qfb_t r;
        fmpz_t root, t1, t2;

        num = qfb_reduced_forms(&forms, -i);
        
        if (num)
        {
           fmpz_init(root);
           fmpz_init(t1);
           fmpz_init(t2);
           qfb_init(r);
              
           fmpz_set_ui(root, i);
           fmpz_root(root, root, 4);

           for (k = 0; k < 20; k++)
           {
              i1 = n_randint(state, num);
              i2 = n_randint(state, num);

              qfb_nucomp(r, forms + i1, forms + i2, root);

              fmpz_mul(t1, r->b, r->b);
              fmpz_mul(t2, r->a, r->c);
              fmpz_mul_2exp(t2, t2, 2);
              fmpz_sub(t2, t1, t2);
              result = (fmpz_cmp_si(t2, -i) == 0);
              if (!result)
              {
                 printf("FAIL:\n");
                 printf("Incorrect discriminant\n");
                 fmpz_print(t2); printf(" should be %ld\n", -i);
                 abort();
              }
           }
           
           fmpz_clear(t1);
           fmpz_clear(t2);
           fmpz_clear(root);
           qfb_clear(r);
        }

        qfb_array_clear(&forms, num);
    }

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

              qfb_nucomp(r, forms + i1, forms + i2, root);
              qfb_reduce(r, r, D);

              qfb_nucomp(r, r, forms + i3, root);
              qfb_reduce(r, r, D);

              qfb_nucomp(s, forms + i2, forms + i3, root);
              qfb_reduce(s, s, D);

              qfb_nucomp(s, forms + i1, s, root);
              qfb_reduce(s, s, D);

              result = (qfb_equal(r, s));
              if (!result)
              {
                 printf("FAIL:\n");
                 printf("Transitivity fails\n");
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
