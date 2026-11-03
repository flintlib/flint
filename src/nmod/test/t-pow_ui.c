/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"
#include "nmod.h"

/* nmod_pow_ui and the algorithms behind it against mpz_powm_ui, for moduli
   of every size (odd and even, around 2^(FLINT_BITS/2)) and exponents of
   every size, including those around the switches between algorithms */
TEST_FUNCTION_START(nmod_pow_ui, state)
{
    slong i;
    mpz_t a_m, n_m, r_m;

    mpz_init(a_m);
    mpz_init(n_m);
    mpz_init(r_m);

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        nmod_t mod;
        ulong n, a, exp, r, s;
        flint_bitcnt_t bits = 1 + n_randint(state, FLINT_BITS);

        switch (n_randint(state, 8))
        {
            case 0:
                n = (bits == FLINT_BITS) ? UWORD_MAX : (UWORD(1) << bits) - 1;
                break;
            case 1:
                /* around 2^(FLINT_BITS/2) and 2^(FLINT_BITS/2 - 2) */
                n = (UWORD(1) << (FLINT_BITS / 2 - 2 * n_randint(state, 2))) - 2 + n_randint(state, 5);
                break;
            default:
                n = n_randbits(state, bits);
        }

        nmod_init(&mod, n);

        a = n_randint(state, 4) ? n_randtest(state) % n : n - 1;

        switch (n_randint(state, 4))
        {
            case 0:
                exp = n_randint(state, 20);
                break;
            case 1:
                /* around 2^k, the possible switches */
                exp = (UWORD(1) << n_randint(state, FLINT_BITS - 1)) - 2 + n_randint(state, 5);
                break;
            default:
                exp = n_randtest(state);
        }

        flint_mpz_set_ui(a_m, a);
        flint_mpz_set_ui(n_m, n);
        flint_mpz_powm_ui(r_m, a_m, exp, n_m);
        s = flint_mpz_get_ui(r_m);

        r = nmod_pow_ui(a, exp, mod);
        if (r != s)
            TEST_FUNCTION_FAIL("nmod_pow_ui\nn = %wu\na = %wu\nexp = %wu\nr = %wu\ns = %wu\n",
                n, a, exp, r, s);

        r = n_powmod2_ui_preinv(a, exp, n, mod.ninv);
        if (r != s)
            TEST_FUNCTION_FAIL("n_powmod2_ui_preinv\nn = %wu\na = %wu\nexp = %wu\nr = %wu\ns = %wu\n",
                n, a, exp, r, s);

        if (exp >= 1)
        {
            r = _nmod_pow_ui_binexp(a, exp, mod);
            if (r != s)
                TEST_FUNCTION_FAIL("_nmod_pow_ui_binexp\nn = %wu\na = %wu\nexp = %wu\nr = %wu\ns = %wu\n",
                    n, a, exp, r, s);

            if (NMOD_IS_HALFLIMB(mod))
            {
                r = _nmod_pow_ui_binexp_halflimb(a, exp, mod);
                if (r != s)
                    TEST_FUNCTION_FAIL("_nmod_pow_ui_binexp_halflimb\nn = %wu\na = %wu\nexp = %wu\nr = %wu\ns = %wu\n",
                        n, a, exp, r, s);
            }

            if (n % 2 == 1 && mod.norm >= FLINT_BITS / 2 + 2)
            {
                r = _nmod_pow_ui_redc_half(a, exp, mod);
                if (r != s)
                    TEST_FUNCTION_FAIL("_nmod_pow_ui_redc_half\nn = %wu\na = %wu\nexp = %wu\nr = %wu\ns = %wu\n",
                        n, a, exp, r, s);
            }

            if (n % 2 == 1)
            {
                r = _nmod_pow_ui_redc(a, exp, mod);
                if (r != s)
                    TEST_FUNCTION_FAIL("_nmod_pow_ui_redc\nn = %wu\na = %wu\nexp = %wu\nr = %wu\ns = %wu\n",
                        n, a, exp, r, s);
            }
        }
    }

    mpz_clear(a_m);
    mpz_clear(n_m);
    mpz_clear(r_m);

    TEST_FUNCTION_END(state);
}
