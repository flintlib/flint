/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

int
main(void)
{
    slong i, j, k;
    FLINT_TEST_INIT(state);

    flint_printf("crt....");
    fflush(stdout);

    /* test internal interface */
    {
        fmpz_multi_crt_t P;
        fmpz ** moduli, ** inputs, * outputs;
        slong moduli_count = 1000;

        moduli = (fmpz **) flint_malloc(moduli_count*sizeof(fmpz *));
        inputs = (fmpz **) flint_malloc(moduli_count*sizeof(fmpz *));
        outputs = (fmpz *) flint_malloc(moduli_count*sizeof(fmpz));

        for (k = 0; k < moduli_count; k++)
        {
            moduli[k] = (fmpz *) flint_malloc(sizeof(fmpz));
            inputs[k] = (fmpz *) flint_malloc(sizeof(fmpz));

            fmpz_init(moduli[k]);
            fmpz_init(inputs[k]);
            fmpz_init(outputs + k);

            fmpz_set_ui(moduli[k], n_nth_prime(k + 2));
            fmpz_sub_ui(inputs[k], moduli[k], 1);
        }

        fmpz_multi_crt_init(P);
        if (!fmpz_multi_crt_precompute_p(P,
                                  (const fmpz * const *) moduli, moduli_count))
        {
            printf("FAIL\n");
            flint_printf("Check simple example\n");
            fflush(stdout);
            flint_abort();
        }

        FLINT_ASSERT(_fmpz_multi_crt_local_size(P) <= moduli_count);

        for (k = 0; k < 1; k++)
        {
            _fmpz_multi_crt_run_p(outputs, P, (const fmpz * const *) inputs);
        }

        for (k = 0; k < moduli_count; k++)
        {
            fmpz_clear(moduli[k]);
            fmpz_clear(inputs[k]);
            fmpz_clear(outputs + k);
            flint_free(moduli[k]);
            flint_free(inputs[k]);
        }

        flint_free(moduli);
        flint_free(inputs);
        flint_free(outputs);

        fmpz_multi_crt_clear(P);
    }

    /* test pointer interface */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_multi_crt_t P;
        fmpz_t t, p;
        slong moduli_count;
        flint_bitcnt_t moduli_bits;
        fmpz ** moduli, ** inputs;
        fmpz_t output;

        fmpz_init(t);
        fmpz_init(p);
        fmpz_init(output);

        fmpz_multi_crt_init(P);

        for (j = 0; j < 4; j++)
        {
            moduli_bits = n_randint(state, 100) + 1;
            moduli_count = n_randint(state, 50) + 1;

            moduli = (fmpz **) flint_malloc(moduli_count*sizeof(fmpz *));
            inputs = (fmpz **) flint_malloc(moduli_count*sizeof(fmpz *));
            for (k = 0; k < moduli_count; k++)
            {
                moduli[k] = (fmpz *) flint_malloc(sizeof(fmpz));
                inputs[k] = (fmpz *) flint_malloc(sizeof(fmpz));
                fmpz_init(moduli[k]);
                fmpz_init(inputs[k]);
                fmpz_randtest_unsigned(moduli[k], state, moduli_bits);
                fmpz_randtest(inputs[k], state, moduli_bits);
            }

            if (fmpz_multi_crt_precompute_p(P,
                                  (const fmpz * const *) moduli, moduli_count))
            {
                fmpz_multi_crt_precomp_p(output, P,
                                                (const fmpz * const *) inputs);

                fmpz_one(p);
                for (k = 0; k < moduli_count; k++)
                {
                    fmpz_mul(p, p, moduli[k]);

                    fmpz_sub(t, output, inputs[k]);
                    if (!fmpz_divisible(t, moduli[k]))
                    {
                        printf("FAIL\n");
                        flint_printf("Check remainder "
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                        fflush(stdout);
                        flint_abort();
                    }
                }
                fmpz_add(output, output, output);
                if (fmpz_cmpabs(output, p) > 0)
                {
                    printf("FAIL\n");
                    flint_printf("Check output size "
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                    fflush(stdout);
                    flint_abort();
                }
            }
            else
            {
                /* check if it was ok to fail on these moduli */
                int ok = 0;

                fmpz_one(p);
                for (k = 0; !ok && k < moduli_count; k++)
                {
                    fmpz_mul(p, p, moduli[k]);
                    if (fmpz_is_zero(moduli[k]) || fmpz_is_pm1(moduli[k]))
                    {
                        ok = 1;
                    }
                }

                for (k = 0; !ok && k < moduli_count; k++)
                {
                    fmpz_divexact(t, p, moduli[k]);
                    fmpz_gcd(t, t, moduli[k]);
                    if (fmpz_cmp_ui(t, 1) > 0)
                    {
                        ok = 1;
                    }
                }

                if (!ok)
                {
                    printf("FAIL\n");
                    flint_printf("Check crt failure i = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }

            for (k = 0; k < moduli_count; k++)
            {
                fmpz_clear(moduli[k]);
                fmpz_clear(inputs[k]);
                flint_free(moduli[k]);
                flint_free(inputs[k]);
            }
            flint_free(moduli);
            flint_free(inputs);
        }

        fmpz_clear(t);
        fmpz_clear(p);
        fmpz_clear(output);

        fmpz_multi_crt_clear(P);
    }

    /* test flat interface */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_multi_crt_t P;
        fmpz_t t, p;
        slong moduli_count;
        flint_bitcnt_t moduli_bits;
        fmpz * moduli, * inputs;
        fmpz_t output;

        fmpz_init(t);
        fmpz_init(p);
        fmpz_init(output);

        fmpz_multi_crt_init(P);

        for (j = 0; j < 4; j++)
        {
            moduli_bits = n_randint(state, 100) + 1;
            moduli_count = n_randint(state, 50) + 1;

            moduli = (fmpz *) flint_malloc(moduli_count*sizeof(fmpz));
            inputs = (fmpz *) flint_malloc(moduli_count*sizeof(fmpz));
            for (k = 0; k < moduli_count; k++)
            {
                fmpz_init(moduli + k);
                fmpz_init(inputs + k);
                fmpz_randtest_unsigned(moduli + k, state, moduli_bits);
                fmpz_randtest(inputs + k, state, moduli_bits);
            }

            if (fmpz_multi_crt_precompute(P, moduli, moduli_count))
            {
                fmpz_multi_crt_precomp(output, P, inputs);

                fmpz_one(p);
                for (k = 0; k < moduli_count; k++)
                {
                    fmpz_mul(p, p, moduli + k);

                    fmpz_sub(t, output, inputs + k);
                    if (!fmpz_divisible(t, moduli + k))
                    {
                        printf("FAIL\n");
                        flint_printf("Check remainder flat "
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                        fflush(stdout);
                        flint_abort();
                    }
                }
                fmpz_add(output, output, output);
                if (fmpz_cmpabs(output, p) > 0)
                {
                    printf("FAIL\n");
                    flint_printf("Check output size flat "
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                    fflush(stdout);
                    flint_abort();
                }
            }
            else
            {
                /* check if it was ok to fail on these moduli */
                int ok = 0;

                fmpz_one(p);
                for (k = 0; !ok && k < moduli_count; k++)
                {
                    fmpz_mul(p, p, moduli + k);
                    if (fmpz_is_zero(moduli + k) || fmpz_is_pm1(moduli + k))
                    {
                        ok = 1;
                    }
                }

                for (k = 0; !ok && k < moduli_count; k++)
                {
                    fmpz_divexact(t, p, moduli + k);
                    fmpz_gcd(t, t, moduli + k);
                    if (fmpz_cmp_ui(t, 1) > 0)
                    {
                        ok = 1;
                    }
                }

                if (!ok)
                {
                    printf("FAIL\n");
                    flint_printf("Check flat crt failure i = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }

            for (k = 0; k < moduli_count; k++)
            {
                fmpz_clear(moduli + k);
                fmpz_clear(inputs + k);
            }
            flint_free(moduli);
            flint_free(inputs);
        }

        fmpz_clear(t);
        fmpz_clear(p);
        fmpz_clear(output);

        fmpz_multi_crt_clear(P);
    }

    /* test lazy interface */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_t t, p;
        slong moduli_count;
        flint_bitcnt_t moduli_bits;
        fmpz * moduli, * inputs;
        fmpz_t output;

        fmpz_init(t);
        fmpz_init(p);
        fmpz_init(output);

        for (j = 0; j < 4; j++)
        {
            moduli_bits = n_randint(state, 100) + 1;
            moduli_count = n_randint(state, 50) + 1;

            moduli = (fmpz *) flint_malloc(moduli_count*sizeof(fmpz));
            inputs = (fmpz *) flint_malloc(moduli_count*sizeof(fmpz));
            for (k = 0; k < moduli_count; k++)
            {
                fmpz_init(moduli + k);
                fmpz_init(inputs + k);
                fmpz_randtest_unsigned(moduli + k, state, moduli_bits);
                fmpz_randtest(inputs + k, state, moduli_bits);
            }

            if (fmpz_multi_crt(output, moduli, inputs, moduli_count))
            {
                fmpz_one(p);
                for (k = 0; k < moduli_count; k++)
                {
                    fmpz_mul(p, p, moduli + k);

                    fmpz_sub(t, output, inputs + k);
                    if (!fmpz_divisible(t, moduli + k))
                    {
                        printf("FAIL\n");
                        flint_printf("Check remainder lazy "
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                        fflush(stdout);
                        flint_abort();
                    }
                }
                fmpz_add(output, output, output);
                if (fmpz_cmpabs(output, p) > 0)
                {
                    printf("FAIL\n");
                    flint_printf("Check output size lazy "
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                    fflush(stdout);
                    flint_abort();
                }
            }
            else
            {
                /* check if it was ok to fail on these moduli */
                int ok = 0;

                fmpz_one(p);
                for (k = 0; !ok && k < moduli_count; k++)
                {
                    fmpz_mul(p, p, moduli + k);
                    if (fmpz_is_zero(moduli + k) || fmpz_is_pm1(moduli + k))
                    {
                        ok = 1;
                    }
                }

                for (k = 0; !ok && k < moduli_count; k++)
                {
                    fmpz_divexact(t, p, moduli + k);
                    fmpz_gcd(t, t, moduli + k);
                    if (fmpz_cmp_ui(t, 1) > 0)
                    {
                        ok = 1;
                    }
                }

                if (!ok)
                {
                    printf("FAIL\n");
                    flint_printf("Check lazy crt failure i = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();                    
                }
            }

            for (k = 0; k < moduli_count; k++)
            {
                fmpz_clear(moduli + k);
                fmpz_clear(inputs + k);
            }
            flint_free(moduli);
            flint_free(inputs);
        }

        fmpz_clear(t);
        fmpz_clear(p);
        fmpz_clear(output);
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}
