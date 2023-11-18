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

TEST_FUNCTION_START(qfb_nucomp, state)
{
    int result;
    qfb * forms;
    slong i, k, i1, i2, i3, num;

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
                 fmpz_print(D); flint_printf(" should be %wd\n", -i);
                 flint_abort();
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
                 flint_abort();
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

    TEST_FUNCTION_END(state);
}
