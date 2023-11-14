/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "qfb.h"

TEST_FUNCTION_START(qfb_exponent_element, state)
{
    int result;
    qfb * forms;
    slong i, k, i1, d, num;

    /* Check correct exponent is returned */
    for (i = 1; i < 500; i++)
    {
        qfb_t r, pow;
        fmpz_t root, D, exp1, exp2;

        d = n_randint(state, 100000);
        num = qfb_reduced_forms(&forms, -d);

        if (num)
        {
           fmpz_init(D);
           fmpz_init(exp1);
           fmpz_init(exp2);
           fmpz_init(root);
           qfb_init(r);
           qfb_init(pow);

           fmpz_set_ui(root, d);
           fmpz_root(root, root, 4);

           for (k = 0; k < 3; k++)
           {
              i1 = n_randint(state, num);
              fmpz_set_si(D, -d);

              fmpz_set_ui(exp1, 1);
              qfb_set(pow, forms + i1);

              while (!qfb_is_principal_form(pow, D))
              {
                 qfb_nucomp(pow, pow, forms + i1, D, root);
                 qfb_reduce(pow, pow, D);
                 fmpz_add_ui(exp1, exp1, 1);
              }

              result = qfb_exponent_element(exp2, forms + i1, D, 1000000, 100000);
              if (!result)
              {
                 printf("FAIL:\n");
                 printf("Exponent not found\n");
                 printf("Discriminant: "); fmpz_print(D); printf("\n");
                 printf("Form: "); qfb_print(forms + i1); printf("\n");
                 flint_abort();
              }

              result = (fmpz_cmp(exp1, exp2) == 0);
              if (!result)
              {
                 printf("FAIL:\n");
                 printf("Incorrect exponent\n");
                 printf("Discriminant: "); fmpz_print(D); printf("\n");
                 printf("Form: "); qfb_print(forms + i1); printf("\n");
                 printf("Exponent ");
                 fmpz_print(exp2); printf(" should be ");
                 fmpz_print(exp1); printf("\n");
                 flint_abort();
              }
           }

           fmpz_clear(root);
           fmpz_clear(D);
           fmpz_clear(exp1);
           fmpz_clear(exp2);
           qfb_clear(r);
        }

        qfb_array_clear(&forms, num);
    }

    TEST_FUNCTION_END(state);
}
