/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "qfb.h"

TEST_FUNCTION_START(qfb_exponent_grh, state)
{
    int result;
    qfb * forms;
    slong i, k, i1, d, num;

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

           result = qfb_exponent_grh(exp, D, 1000000, 100000);
           if (!result)
           {
              printf("FAIL:\n");
              printf("Exponent not computed\n");
              printf("Discriminant: "); fmpz_print(D); printf("\n");
              flint_abort();
           }

           e = fmpz_get_si(exp);
           result = ((num % e) == 0);
           if (!result)
           {
              printf("FAIL:\n");
              printf("Exponent does not divide class number\n");
              printf("Discriminant: "); fmpz_print(D); printf("\n");
              flint_printf("%wd does not divide %wd\n", e, num);
              flint_abort();
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
                 flint_printf("Exponent: %wd\n", e);
                 flint_abort();
              }
           }

           qfb_clear(pow);
           fmpz_clear(D);
           fmpz_clear(exp);
        }

        qfb_array_clear(&forms, num);
    }

    /* Small B1 (fewer than 1024 primes up to B1): if an exponent is
       returned, it must annihilate the forms */
    for (i = 0; i < 20; i++)
    {
        fmpz_t D, exp;
        qfb_t pow;
        ulong B1 = 2 + n_randint(state, 8000);

        d = (i == 0) ? 1000003 : n_randint(state, 100000) + 1;
        if (i == 0)
            B1 = 1000;
        num = qfb_reduced_forms(&forms, -d);

        if (num)
        {
           fmpz_init(D);
           fmpz_init(exp);
           qfb_init(pow);

           fmpz_set_si(D, -d);

           if (qfb_exponent_grh(exp, D, B1, 100))
           {
              for (k = 0; k < FLINT_MIN(num, 5); k++)
              {
                 qfb_pow(pow, forms + k, D, exp);
                 qfb_reduce(pow, pow, D);
                 if (!qfb_is_principal_form(pow, D))
                 {
                    printf("FAIL (small B1):\n");
                    printf("Discriminant: "); fmpz_print(D); printf("\n");
                    flint_printf("B1 = %wu\n", B1);
                    flint_abort();
                 }
              }
           }

           qfb_clear(pow);
           fmpz_clear(D);
           fmpz_clear(exp);
        }

        qfb_array_clear(&forms, num);
    }

    TEST_FUNCTION_END(state);
}
