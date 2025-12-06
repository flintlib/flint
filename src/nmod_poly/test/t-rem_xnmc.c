/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_rem_xnmc, state)
{
    int i, result = 1;

    for (i = 0; i < 5000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        nmod_poly_t rok, r1;
        nn_ptr rvec1;
        ulong modn, len, n, c;
        nmod_t mod;

        /* modulus */
        if (i < 128)
            modn = n_randtest_bits(state, FLINT_BITS);
        else if (i < 256)
            modn = n_randtest_bits(state, FLINT_BITS-1);
        else if (i < 384)
            modn = n_randtest_bits(state, FLINT_BITS-2);
        else
            modn = n_randtest_not_zero(state);
        nmod_init(&mod, modn);

        /* initialize polynomials */
        nmod_poly_init_mod(a, mod);
        nmod_poly_init_mod(b, mod);
        nmod_poly_init_mod(rok, mod);
        nmod_poly_init_mod(r1, mod);

        /* pick parameters */
        len = 5 + n_randint(state, 100);    /* poly length, [5, 105) */
        n = 1 + n_randint(state, 30);       /* divisor degree, [1, 30) */

        if (i % 16 == 0 && modn > 1)
            c = 1;
        else if (i % 16 == 1)
            c = modn - 1;
        else if (i % 16 == 2)
            c = 0;
        else
            c = n_randint(state, modn);

        /* random a of length <= len; rvec1 of length == n */
        nmod_poly_randtest(a, state, len);
        rvec1 = _nmod_vec_init(n);

        /* divisor xnmc:  b == x**n - c */
        nmod_poly_set_coeff_ui(b, n, 1);
        nmod_poly_set_coeff_ui(b, 0, n_negmod(c, modn));

        // correct quotient/remainder
        if (modn == 1)
            nmod_poly_zero(rok);
        else
            nmod_poly_rem(rok, a, b);

        // general interface
        nmod_poly_rem_xnmc(r1, a, n, c);

        result = (nmod_poly_equal(rok, r1));
        if (!result)
        {
            flint_printf("FAIL (poly 1):\n");
            flint_printf("a->length = %wd, len = %wu, modn = %wu\n", a->length, len, a->mod.n);
            flint_printf("n = %wu, c = %wu, poly a :\n", n, c);
            nmod_poly_print(a), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        /* testing in-place underscore version, requires length >= n */
        if (a->length >= (slong)n)
        {
            /* general function */
            _nmod_poly_rem_xnmc(rvec1, a->coeffs, a->length, n, c, a->mod);

            result = (_nmod_vec_equal(rvec1, rok->coeffs, rok->length) && _nmod_vec_is_zero(rvec1 + rok->length, n - rok->length));
            if (!result)
            {
                flint_printf("FAIL (vec):\n");
                flint_printf("a->length = %wd, len = %wu, modn = %wu\n", a->length, len, a->mod.n);
                flint_printf("n = %wu, c = %wu, input vec rvec :\n", n, c);
                _nmod_vec_print(a->coeffs, a->length, mod), flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }

            /* special cases 1, -1 */
            if (c == 1 || c == modn - 1)
            {
                if (c == 1)
                    _nmod_poly_rem_xnm1(rvec1, a->coeffs, a->length, n, a->mod.n);
                else if (c == modn-1)
                    _nmod_poly_rem_xnp1(rvec1, a->coeffs, a->length, n, a->mod.n);

                result = (_nmod_vec_equal(rvec1, rok->coeffs, rok->length) && _nmod_vec_is_zero(rvec1 + rok->length, n - rok->length));
                if (!result)
                {
                    flint_printf("FAIL (vec, +1 / - 1):\n");
                    flint_printf("a->length = %wd, len = %wu, modn = %wu\n", a->length, len, a->mod.n);
                    flint_printf("n = %wu, c = %wu, input vec rvec :\n", n, c);
                    _nmod_vec_print(a->coeffs, a->length, mod), flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            /* cases using precomp */
            if (NMOD_CAN_USE_SHOUP(mod))
            {
                ulong c_precomp = n_mulmod_precomp_shoup(c, modn);
                int lazy = 0;

#if FLINT_BITS == 64
                if (modn <= UWORD(6148914691236517205))
#else // FLINT_BITS == 32
                if (modn <= UWORD(1431655765))
#endif
                {
                    lazy = 1;
                    _nmod_poly_rem_xnmc_precomp_lazy(rvec1, a->coeffs, a->length, n, c, c_precomp, modn);
                    for (ulong ii = 0; ii < n; ii++)
                    {
                        if (rvec1[ii] >= 2*modn)
                            rvec1[ii] -= 2*modn;
                        else if (rvec1[ii] >= modn)
                            rvec1[ii] -= modn;
                    }
                }
                else
                    _nmod_poly_rem_xnmc_precomp(rvec1, a->coeffs, a->length, n, c, c_precomp, modn);

                result = (_nmod_vec_equal(rvec1, rok->coeffs, rok->length) && _nmod_vec_is_zero(rvec1 + rok->length, n - rok->length));
                if (!result)
                {
                    if (lazy)
                        flint_printf("FAIL (vec, precomp lazy):\n");
                    else
                        flint_printf("FAIL (vec, precomp):\n");
                    flint_printf("a->length = %wd, len = %wu, modn = %wu\n", a->length, len, a->mod.n);
                    flint_printf("n = %wu, c = %wu, input vec rvec :\n", n, c);
                    _nmod_vec_print(a->coeffs, a->length, mod), flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(rok);
        nmod_poly_clear(r1);
        _nmod_vec_clear(rvec1);
    }

    TEST_FUNCTION_END(state);
}
