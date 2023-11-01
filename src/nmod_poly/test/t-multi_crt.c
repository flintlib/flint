/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

TEST_FUNCTION_START(nmod_poly_multi_crt, state)
{
    slong i, j, k;

    /* test internal interface */
    {
        nmod_poly_multi_crt_t P;
        nmod_poly_struct ** moduli, ** inputs, * outputs;
        mp_limb_t modulus = 1009;
        slong moduli_count = 1000;

        moduli = (nmod_poly_struct **) flint_malloc(moduli_count*sizeof(nmod_poly_struct *));
        inputs = (nmod_poly_struct **) flint_malloc(moduli_count*sizeof(nmod_poly_struct *));
        outputs = (nmod_poly_struct *) flint_malloc(moduli_count*sizeof(nmod_poly_struct));

        for (k = 0; k < moduli_count; k++)
        {
            moduli[k] = (nmod_poly_struct *) flint_malloc(sizeof(nmod_poly_struct));
            inputs[k] = (nmod_poly_struct *) flint_malloc(sizeof(nmod_poly_struct));

            nmod_poly_init(moduli[k], modulus);
            nmod_poly_init(inputs[k], modulus);
            nmod_poly_init(outputs + k, modulus);

            nmod_poly_set_coeff_ui(moduli[k], 1, 1);
            nmod_poly_set_coeff_ui(moduli[k], 0, modulus - k);

            nmod_poly_set_coeff_ui(inputs[k], 0, (k*k)^k);
        }

        nmod_poly_multi_crt_init(P);
        if (!nmod_poly_multi_crt_precompute_p(P,
                     (const nmod_poly_struct * const *) moduli, moduli_count))
        {
            printf("FAIL\n");
            flint_printf("Check simple example\n");
            fflush(stdout);
            flint_abort();
        }

        FLINT_ASSERT(_nmod_poly_multi_crt_local_size(P) <= moduli_count);

        for (k = 0; k < 1; k++)
        {
            _nmod_poly_multi_crt_run_p(outputs, P,
                                    (const nmod_poly_struct * const *) inputs);
        }

        for (k = 0; k < moduli_count; k++)
        {
            nmod_poly_clear(moduli[k]);
            nmod_poly_clear(inputs[k]);
            nmod_poly_clear(outputs + k);
            flint_free(moduli[k]);
            flint_free(inputs[k]);
        }

        flint_free(moduli);
        flint_free(inputs);
        flint_free(outputs);

        nmod_poly_multi_crt_clear(P);
    }

    /* test pointer interface */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_poly_multi_crt_t P;
        nmod_poly_t t, p;
        slong total_degree, moduli_length, moduli_count;
        nmod_poly_struct ** moduli, ** inputs;
        nmod_poly_t output;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_poly_init(t, modulus);
        nmod_poly_init(p, modulus);
        nmod_poly_init(output, modulus);

        nmod_poly_multi_crt_init(P);

        for (j = 0; j < 4; j++)
        {
            moduli_length = n_randint(state, 20) + 1;
            moduli_count = n_randint(state, 50) + 1;

            moduli = (nmod_poly_struct **) flint_malloc(moduli_count
                                                  *sizeof(nmod_poly_struct *));
            inputs = (nmod_poly_struct **) flint_malloc(moduli_count
                                                  *sizeof(nmod_poly_struct *));
            for (k = 0; k < moduli_count; k++)
            {
                moduli[k] = (nmod_poly_struct *) flint_malloc(
                                                      sizeof(nmod_poly_struct));
                inputs[k] = (nmod_poly_struct *) flint_malloc(
                                                      sizeof(nmod_poly_struct));
                nmod_poly_init(moduli[k], modulus);
                nmod_poly_init(inputs[k], modulus);
                nmod_poly_randtest(moduli[k], state, moduli_length);
                nmod_poly_randtest(inputs[k], state, moduli_length);
            }

            if (nmod_poly_multi_crt_precompute_p(P,
                      (const nmod_poly_struct * const *) moduli, moduli_count))
            {
                nmod_poly_multi_crt_precomp_p(output, P,
                                    (const nmod_poly_struct * const *) inputs);

                total_degree = 0;
                for (k = 0; k < moduli_count; k++)
                {
                    total_degree += nmod_poly_degree(moduli[k]);

                    nmod_poly_sub(t, output, inputs[k]);
                    nmod_poly_rem(t, t, moduli[k]);
                    if (!nmod_poly_is_zero(t))
                    {
                        printf("FAIL\n");
                        flint_printf("Check remainder "
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                        fflush(stdout);
                        flint_abort();
                    }
                }
                if (nmod_poly_degree(output) >= total_degree)
                {
                    printf("FAIL\n");
                    flint_printf("Check output degree "
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                    fflush(stdout);
                    flint_abort();
                }
            }
            else
            {
                /* check if it was ok to fail on these moduli */
                int ok = 0;

                nmod_poly_one(p);
                for (k = 0; !ok && k < moduli_count; k++)
                {
                    nmod_poly_mul(p, p, moduli[k]);
                    if (nmod_poly_degree(moduli[k]) < 1)
                    {
                        ok = 1;
                    }
                }

                for (k = 0; !ok && k < moduli_count; k++)
                {
                    nmod_poly_div(t, p, moduli[k]);
                    nmod_poly_gcd(t, t, moduli[k]);
                    if (nmod_poly_degree(t) > 0)
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
                nmod_poly_clear(moduli[k]);
                nmod_poly_clear(inputs[k]);
                flint_free(moduli[k]);
                flint_free(inputs[k]);
            }
            flint_free(moduli);
            flint_free(inputs);
        }

        nmod_poly_clear(t);
        nmod_poly_clear(p);
        nmod_poly_clear(output);

        nmod_poly_multi_crt_clear(P);
    }

    /* test flat interface */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_poly_multi_crt_t P;
        nmod_poly_t t, p;
        slong total_degree, moduli_length, moduli_count;
        nmod_poly_struct * moduli, * inputs;
        nmod_poly_t output;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_poly_init(t, modulus);
        nmod_poly_init(p, modulus);
        nmod_poly_init(output, modulus);

        nmod_poly_multi_crt_init(P);

        for (j = 0; j < 4; j++)
        {
            moduli_length = n_randint(state, 20) + 1;
            moduli_count = n_randint(state, 50) + 1;

            moduli = (nmod_poly_struct *) flint_malloc(moduli_count
                                                  *sizeof(nmod_poly_struct));
            inputs = (nmod_poly_struct *) flint_malloc(moduli_count
                                                  *sizeof(nmod_poly_struct));
            for (k = 0; k < moduli_count; k++)
            {
                nmod_poly_init(moduli + k, modulus);
                nmod_poly_init(inputs + k, modulus);
                nmod_poly_randtest(moduli + k, state, moduli_length);
                nmod_poly_randtest(inputs + k, state, moduli_length);
            }

            if (nmod_poly_multi_crt_precompute(P, moduli, moduli_count))
            {
                nmod_poly_multi_crt_precomp(output, P, inputs);

                total_degree = 0;
                for (k = 0; k < moduli_count; k++)
                {
                    total_degree += nmod_poly_degree(moduli + k);

                    nmod_poly_sub(t, output, inputs + k);
                    nmod_poly_rem(t, t, moduli + k);
                    if (!nmod_poly_is_zero(t))
                    {
                        printf("FAIL\n");
                        flint_printf("Check remainder flat "
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                        fflush(stdout);
                        flint_abort();
                    }
                }
                if (nmod_poly_degree(output) >= total_degree)
                {
                    printf("FAIL\n");
                    flint_printf("Check output degree flat "
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                    fflush(stdout);
                    flint_abort();
                }
            }
            else
            {
                /* check if it was ok to fail on these moduli */
                int ok = 0;

                nmod_poly_one(p);
                for (k = 0; !ok && k < moduli_count; k++)
                {
                    nmod_poly_mul(p, p, moduli + k);
                    if (nmod_poly_degree(moduli + k) < 1)
                    {
                        ok = 1;
                    }
                }

                for (k = 0; !ok && k < moduli_count; k++)
                {
                    nmod_poly_div(t, p, moduli + k);
                    nmod_poly_gcd(t, t, moduli + k);
                    if (nmod_poly_degree(t) > 0)
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
                nmod_poly_clear(moduli + k);
                nmod_poly_clear(inputs + k);
            }
            flint_free(moduli);
            flint_free(inputs);
        }

        nmod_poly_clear(t);
        nmod_poly_clear(p);
        nmod_poly_clear(output);

        nmod_poly_multi_crt_clear(P);
    }

    /* test lazy interface */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_poly_t t, p;
        slong total_degree, moduli_length, moduli_count;
        nmod_poly_struct * moduli, * inputs;
        nmod_poly_t output;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_poly_init(t, modulus);
        nmod_poly_init(p, modulus);
        nmod_poly_init(output, modulus);

        for (j = 0; j < 4; j++)
        {
            moduli_length = n_randint(state, 20) + 1;
            moduli_count = n_randint(state, 50) + 1;

            moduli = (nmod_poly_struct *) flint_malloc(moduli_count
                                                  *sizeof(nmod_poly_struct));
            inputs = (nmod_poly_struct *) flint_malloc(moduli_count
                                                  *sizeof(nmod_poly_struct));
            for (k = 0; k < moduli_count; k++)
            {
                nmod_poly_init(moduli + k, modulus);
                nmod_poly_init(inputs + k, modulus);
                nmod_poly_randtest(moduli + k, state, moduli_length);
                nmod_poly_randtest(inputs + k, state, moduli_length);
            }

            if (nmod_poly_multi_crt(output, moduli, inputs, moduli_count))
            {
                total_degree = 0;
                for (k = 0; k < moduli_count; k++)
                {
                    total_degree += nmod_poly_degree(moduli + k);

                    nmod_poly_sub(t, output, inputs + k);
                    nmod_poly_rem(t, t, moduli + k);
                    if (!nmod_poly_is_zero(t))
                    {
                        printf("FAIL\n");
                        flint_printf("Check remainder lazy "
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                        fflush(stdout);
                        flint_abort();
                    }
                }
                if (nmod_poly_degree(output) >= total_degree)
                {
                    printf("FAIL\n");
                    flint_printf("Check output degree lazy "
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
                    fflush(stdout);
                    flint_abort();
                }
            }
            else
            {
                /* check if it was ok to fail on these moduli */
                int ok = 0;

                nmod_poly_one(p);
                for (k = 0; !ok && k < moduli_count; k++)
                {
                    nmod_poly_mul(p, p, moduli + k);
                    if (nmod_poly_degree(moduli + k) < 1)
                    {
                        ok = 1;
                    }
                }

                for (k = 0; !ok && k < moduli_count; k++)
                {
                    nmod_poly_div(t, p, moduli + k);
                    nmod_poly_gcd(t, t, moduli + k);
                    if (nmod_poly_degree(t) > 0)
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
                nmod_poly_clear(moduli + k);
                nmod_poly_clear(inputs + k);
            }
            flint_free(moduli);
            flint_free(inputs);
        }

        nmod_poly_clear(t);
        nmod_poly_clear(p);
        nmod_poly_clear(output);
    }

    TEST_FUNCTION_END(state);
}
