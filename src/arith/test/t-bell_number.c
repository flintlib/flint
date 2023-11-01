/*
    Copyright (C) 2011, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "arith.h"

#define STRESS_TEST 0

#if STRESS_TEST
#define MAXN 100000
#define MAXN_VEC 10000
#else
#define MAXN 4000
#define MAXN_VEC 500
#endif

TEST_FUNCTION_START(arith_bell_number, state)
{
    {
        slong len, prev_len;
        fmpz * vb1, * vb2;
        fmpz_t b;
        mp_ptr vnb, vnr;
        slong n, iter;
        ulong nb;
        nmod_t mod;
        prev_len = 0;

        fmpz_init(b);

        for (n = 0; n < MAXN; n += n_randint(state, n / 4 + 2))
        {
#if STRESS_TEST
            flint_printf("%wd\n", n);
#endif

            arith_bell_number(b, n);

            for (iter = 0; iter < 3 + 10 * STRESS_TEST; iter++)
            {
                nmod_init(&mod, n_randtest_not_zero(state));
                nb = arith_bell_number_nmod(n, mod);

                if (nb != fmpz_fdiv_ui(b, mod.n))
                {
                    flint_printf("FAIL (vs nmod, n = %wd)\n", n);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        for (len = 0; len < MAXN_VEC; len = FLINT_MAX(len + 1, len * 1.25))
        {
            vb1 = _fmpz_vec_init(len);
            vb2 = _fmpz_vec_init(len);
            vnb = _nmod_vec_init(len);
            vnr = _nmod_vec_init(len);

            arith_bell_number_vec_recursive(vb1, len);
            arith_bell_number_vec_multi_mod(vb2, len);

            if (!_fmpz_vec_equal(vb1, vb2, len))
            {
                flint_printf("FAIL (len = %wd)\n", len);
                fflush(stdout);
                flint_abort();
            }

            for (n = prev_len; n < len; n++)
            {
#if STRESS_TEST
                flint_printf("%wd\n", n);
#endif
                if (n < 5000)
                {
                    arith_bell_number_dobinski(b, n);
                    if (!fmpz_equal(vb1 + n, b))
                    {
                        flint_printf("FAIL (dobinski, n = %wd)\n", n);
                        fflush(stdout);
                        flint_abort();
                    }
                }

                arith_bell_number_multi_mod(b, n);
                if (!fmpz_equal(vb1 + n, b))
                {
                    flint_printf("FAIL (multi_mod, n = %wd)\n", n);
                    fflush(stdout);
                    flint_abort();
                }
            }

            for (iter = 0; iter < 30; iter++)
            {
                nmod_init(&mod, n_randtest_not_zero(state));
                arith_bell_number_nmod_vec(vnb, len, mod);

                _fmpz_vec_get_nmod_vec(vnr, vb1, len, mod);

                if (!_nmod_vec_equal(vnr, vnb, len))
                {
                    flint_printf("FAIL (nmod_vec, len = %wd)\n", len);
                    fflush(stdout);
                    flint_abort();
                }

                if (len)
                {
                    n = n_randint(state, len);
                    nb = arith_bell_number_nmod(n, mod);

                    if (nb != fmpz_fdiv_ui(vb1 + n, mod.n))
                    {
                        flint_printf("FAIL (nmod n = %wd)\n", n);
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }

            _fmpz_vec_clear(vb1, len);
            _fmpz_vec_clear(vb2, len);
            _nmod_vec_clear(vnb);
            _nmod_vec_clear(vnr);

            prev_len = len;
        }

        fmpz_clear(b);
    }

    TEST_FUNCTION_END(state);
}
