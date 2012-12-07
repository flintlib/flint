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

    Copyright (C) 2012 William Hart

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
    long i, k, i1, num;

    printf("nudupl....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 1; i < 10000; i++) 
    {
        qfb_t r, s;
        fmpz_t root, D;

        num = qfb_reduced_forms(&forms, -i);
        
        if (num)
        {
           fmpz_init(root);
           fmpz_init(D);
           qfb_init(r);
           qfb_init(s);
              
           fmpz_set_ui(root, i);
           fmpz_root(root, root, 4);

           for (k = 0; k < 20; k++)
           {
              i1 = n_randint(state, num);
              
              qfb_discriminant(D, forms + i1);

              qfb_nucomp(r, forms + i1, forms + i1, root);
              qfb_reduce(r, r, D);
              
              qfb_nudupl(s, forms + i1, root);
              qfb_reduce(s, s, D);
              
              result = (qfb_equal(r, s));
              if (!result)
              {
                 printf("FAIL:\n");
                 printf("NUDUPL does not agree with NUCOMP\n");
                 printf("f = "); qfb_print(forms + i1); printf("\n");
                 qfb_print(r); printf(" != "); qfb_print(s);
                 abort();
              }
           }
           
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
