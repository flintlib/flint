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

TEST_FUNCTION_START(nmod_poly_divrem_xnmc, state)
{
    int i, result = 1;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        nmod_poly_t qok, rok, q1, r1;
        nn_ptr qr1;
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
        nmod_poly_init_mod(qok, mod);
        nmod_poly_init_mod(rok, mod);
        nmod_poly_init_mod(q1, mod);
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

        /* random a of length <= len; qr1 with same coefficients */
        nmod_poly_randtest(a, state, len);
        qr1 = _nmod_vec_init(a->length);
        _nmod_vec_set(qr1, a->coeffs, a->length);

        /* divisor xnmc:  b == x**n - c */
        nmod_poly_set_coeff_ui(b, n, 1);
        nmod_poly_set_coeff_ui(b, 0, n_negmod(c, modn));

        // correct quotient/remainder
        nmod_poly_divrem(qok, rok, a, b);

        // general interface
        nmod_poly_divrem_xnmc(q1, r1, a, n, c);

        result = (nmod_poly_equal(qok, q1) && nmod_poly_equal(rok, r1));
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
            _nmod_poly_divrem_xnmc(qr1, qr1, a->length, n, c, a->mod);

            result = (_nmod_vec_equal(qr1, rok->coeffs, rok->length) && _nmod_vec_is_zero(qr1 + rok->length, n - rok->length)
                      && _nmod_vec_equal(qr1+n, qok->coeffs, qok->length));
            if (!result)
            {
                flint_printf("FAIL (vec):\n");
                flint_printf("a->length = %wd, len = %wu, modn = %wu\n", a->length, len, a->mod.n);
                flint_printf("n = %wu, c = %wu, input vec qr :\n", n, c);
                _nmod_vec_print(a->coeffs, a->length, mod), flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }

            /* special cases 1, -1 */
            _nmod_vec_set(qr1, a->coeffs, a->length);
            if (c == 1 || c == modn - 1)
            {
                if (c == 1)
                    _nmod_poly_divrem_xnm1(qr1, qr1, a->length, n, a->mod.n);
                else if (c == modn-1)
                    _nmod_poly_divrem_xnp1(qr1, qr1, a->length, n, a->mod.n);

                result = (_nmod_vec_equal(qr1, rok->coeffs, rok->length) && _nmod_vec_is_zero(qr1 + rok->length, n - rok->length)
                        && _nmod_vec_equal(qr1+n, qok->coeffs, qok->length));
                if (!result)
                {
                    flint_printf("FAIL (vec, +1 / - 1):\n");
                    flint_printf("a->length = %wd, len = %wu, modn = %wu\n", a->length, len, a->mod.n);
                    flint_printf("n = %wu, c = %wu, input vec qr :\n", n, c);
                    _nmod_vec_print(a->coeffs, a->length, mod), flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            /* cases using precomp */
            _nmod_vec_set(qr1, a->coeffs, a->length);
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
                    _nmod_poly_divrem_xnmc_precomp_lazy(qr1, qr1, a->length, n, c, c_precomp, modn);
                    for (long ii = 0; ii < a->length; ii++)
                    {
                        if (qr1[ii] >= 2*modn)
                            qr1[ii] -= 2*modn;
                        else if (qr1[ii] >= modn)
                            qr1[ii] -= modn;
                    }
                }
                else
                    _nmod_poly_divrem_xnmc_precomp(qr1, qr1, a->length, n, c, c_precomp, modn);

                result = (_nmod_vec_equal(qr1, rok->coeffs, rok->length) && _nmod_vec_is_zero(qr1 + rok->length, n - rok->length)
                        && _nmod_vec_equal(qr1+n, qok->coeffs, qok->length));
                if (!result)
                {
                    if (lazy)
                        flint_printf("FAIL (vec, precomp lazy):\n");
                    else
                        flint_printf("FAIL (vec, precomp):\n");
                    flint_printf("a->length = %wd, len = %wu, modn = %wu\n", a->length, len, a->mod.n);
                    flint_printf("n = %wu, c = %wu, input vec qr :\n", n, c);
                    _nmod_vec_print(a->coeffs, a->length, mod), flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(qok);
        nmod_poly_clear(rok);
        nmod_poly_clear(q1);
        nmod_poly_clear(r1);
        _nmod_vec_clear(qr1);
    }

    TEST_FUNCTION_END(state);
}
