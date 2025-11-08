/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_fft/impl_macros_tft.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "n_fft.h"
#include "n_fft/impl.h"

#define MAX_EVAL_DEPTH 10

void extend_recurrence(nn_ptr p, ulong len, nn_ptr rec, ulong d, nmod_t mod)
{
    if (len > d)
    {
        for (ulong i = 0; i < len - d; i++)
        {
            p[d+i] = 0;
            for (ulong j = 0; j < d; j++)
            {
                p[d+i] = nmod_sub(p[d+i], nmod_mul(rec[j], p[j+i], mod), mod);
            }
        }
    }
}

TEST_FUNCTION_START(n_fft_special_divrems, state)
{
    int i, result = 1;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        // take some FFT prime p with max_depth >= MAX_EVAL_DEPTH
        ulong max_depth, prime;

        // half of tests == fixed large prime, close to limit
        // 62 bits: prime = 4611686018427322369 == 2**62 - 2**16 + 1
        // 30 bits: prime = 1073479681 == 2**30 - 2**18 + 1
        if (i > 50 * flint_test_multiplier())
#if FLINT_BITS == 64
            prime = UWORD(4611686018427322369);
#else // FLINT_BITS == 32
            prime = UWORD(1073479681);
#endif
        else
        {
            max_depth = MAX_EVAL_DEPTH + n_randint(state, 6);
            prime = 1 + (UWORD(1) << max_depth);
            while (! n_is_prime(prime))
                prime += (UWORD(1) << max_depth);
        }
        max_depth = flint_ctz(prime-1);

        nmod_t mod;
        nmod_init(&mod, prime);

        // init FFT root tables
        n_fft_ctx_t F;
        n_fft_ctx_init2(F, MAX_EVAL_DEPTH, prime);
        n_fft_args_t Fargs;
        n_fft_set_args(Fargs, F->mod, F->tab_w);

        // retrieve roots, used later for building master poly
        nn_ptr roots = flint_malloc((UWORD(1) << MAX_EVAL_DEPTH) * sizeof(ulong));
        for (ulong k = 0; k < (UWORD(1) << (MAX_EVAL_DEPTH-1)); k++)
        {
            roots[2*k] = F->tab_w[2*k];
            roots[2*k+1] = prime - F->tab_w[2*k];  // < prime since F->tab_w[2*k] != 0
        }

        for (ulong depth = 1; depth < MAX_EVAL_DEPTH; depth++)
        {
            const ulong len = (UWORD(1) << depth);

            nmod_poly_t a, b;
            nmod_poly_t qok, rok;

            /* parameter d */
            ulong d = 1+len/2 + n_randint(state, len/2);  /* (len/2, len] */

            /* input vector and copy */
            ulong ilen = d+1 + n_randint(state, 5*len);  /* [d+1, ~ 5*len] */
            nn_ptr p = _nmod_vec_init(ilen);
            _nmod_vec_randtest(p, state, ilen, mod);
            while (p[ilen-1] == 0)
                p[ilen-1] = n_randint(state, prime);
            nn_ptr pp = _nmod_vec_init(ilen);
            nn_ptr pok = _nmod_vec_init(ilen);
            nn_ptr rec = _nmod_vec_init(d+1);
            /* flint_printf("depth = %wu | %wu, %wu...\n", depth, d, ilen); */

            /* initialize polynomials */
            nmod_poly_init_mod(a, mod);
            nmod_poly_init_mod(b, mod);
            nmod_poly_init_mod(qok, mod);
            nmod_poly_init_mod(rok, mod);

            /* fill a with p */
            nmod_poly_fit_length(a, ilen);
            _nmod_poly_set_length(a, ilen);
            _nmod_vec_set(a->coeffs, p, ilen);

            /* divrem circulant, general, v0 */
            {
                ulong c = n_randint(state, prime);
                ulong c_precomp = n_mulmod_precomp_shoup(c, prime);

                /* divisor b == x**d - c */
                nmod_poly_zero(b);
                nmod_poly_set_coeff_ui(b, d, 1);
                nmod_poly_set_coeff_ui(b, 0, n_negmod(c, prime));

                nmod_poly_divrem(qok, rok, a, b);

                _nmod_vec_set(pp, p, ilen);
                _nmod_poly_divrem_circulant_lazy_4_4_v0(pp, ilen, d, c, c_precomp, prime, 2*prime);

                for (ulong i = 0; i < ilen; i++)
                {
                    if (pp[i] >= 2*prime)
                        pp[i] -= 2*prime;
                    if (pp[i] >= prime)
                        pp[i] -= prime;
                }

                result = (_nmod_vec_equal(pp, rok->coeffs, rok->length) && _nmod_vec_is_zero(pp + rok->length, d - rok->length)
                        && _nmod_vec_equal(pp+d, qok->coeffs, qok->length));
                if (!result)
                {
                    flint_printf("FAIL (_nmod_poly_divrem_circulant_lazy_4_4_v0):\n");
                    flint_printf("a->length = %wd, ilen = %wu, modulus = %wu\n", a->length, ilen, prime);
                    flint_printf("d = %wu, c = %wu, input vec :\n", d, c);
                    _nmod_vec_print(p, ilen, mod), flint_printf("\n\n");
                    _nmod_vec_print(rok->coeffs, rok->length, mod), flint_printf("\n\n");
                    _nmod_vec_print(qok->coeffs, qok->length, mod), flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            /* divrem circulant, general, v1 */
            {
                ulong c = n_randint(state, prime);
                ulong c_precomp = n_mulmod_precomp_shoup(c, prime);

                /* divisor b == x**d - c */
                nmod_poly_zero(b);
                nmod_poly_set_coeff_ui(b, d, 1);
                nmod_poly_set_coeff_ui(b, 0, n_negmod(c, prime));

                nmod_poly_divrem(qok, rok, a, b);

                _nmod_vec_set(pp, p, ilen);
                _nmod_poly_divrem_circulant_lazy_4_4(pp, ilen, d, c, c_precomp, prime, 2*prime);

                for (ulong i = 0; i < ilen; i++)
                {
                    if (pp[i] >= 2*prime)
                        pp[i] -= 2*prime;
                    if (pp[i] >= prime)
                        pp[i] -= prime;
                }

                result = (_nmod_vec_equal(pp, rok->coeffs, rok->length) && _nmod_vec_is_zero(pp + rok->length, d - rok->length)
                        && _nmod_vec_equal(pp+d, qok->coeffs, qok->length));
                if (!result)
                {
                    flint_printf("FAIL (_nmod_poly_divrem_circulant_lazy_4_4):\n");
                    flint_printf("a->length = %wd, ilen = %wu, modulus = %wu\n", a->length, ilen, prime);
                    flint_printf("d = %wu, c = %wu, input vec :\n", d, c);
                    _nmod_vec_print(p, ilen, mod), flint_printf("\n\n");
                    _nmod_vec_print(rok->coeffs, rok->length, mod), flint_printf("\n\n");
                    _nmod_vec_print(qok->coeffs, qok->length, mod), flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            /* divrem circulant, 1, v0 */
            {
                ulong c = 1;

                /* divisor b == x**d - 1 */
                nmod_poly_zero(b);
                nmod_poly_set_coeff_ui(b, d, 1);
                nmod_poly_set_coeff_ui(b, 0, prime - 1);

                nmod_poly_divrem(qok, rok, a, b);

                _nmod_vec_set(pp, p, ilen);
                _nmod_poly_divrem_circulant1(pp, ilen, d, prime);

                result = (_nmod_vec_equal(pp, rok->coeffs, rok->length) && _nmod_vec_is_zero(pp + rok->length, d - rok->length)
                        && _nmod_vec_equal(pp+d, qok->coeffs, qok->length));
                if (!result)
                {
                    flint_printf("FAIL (_nmod_poly_divrem_circulant1):\n");
                    flint_printf("a->length = %wd, ilen = %wu, modulus = %wu\n", a->length, ilen, prime);
                    flint_printf("d = %wu, c = %wu, input vec :\n", d, c);
                    _nmod_vec_print(p, ilen, mod), flint_printf("\n\n");
                    _nmod_vec_print(rok->coeffs, rok->length, mod), flint_printf("\n\n");
                    _nmod_vec_print(qok->coeffs, qok->length, mod), flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            /* divrem circulant, 1, v1 */
            {
                ulong c = 1;

                /* divisor b == x**d - 1 */
                nmod_poly_zero(b);
                nmod_poly_set_coeff_ui(b, d, 1);
                nmod_poly_set_coeff_ui(b, 0, prime - 1);

                nmod_poly_rem(rok, a, b);

                _nmod_vec_set(pp, p, ilen);
                _nmod_poly_divrem_circulant1_v1(pp, ilen, d, prime);

                result = (_nmod_vec_equal(pp, rok->coeffs, rok->length) && _nmod_vec_is_zero(pp + rok->length, d - rok->length)
                        && _nmod_vec_equal(pp+d, qok->coeffs, qok->length));
                if (!result)
                {
                    flint_printf("FAIL (_nmod_poly_divrem_circulant1_v1):\n");
                    flint_printf("a->length = %wd, ilen = %wu, modulus = %wu\n", a->length, ilen, prime);
                    flint_printf("d = %wu, c = %wu, input vec :\n", d, c);
                    _nmod_vec_print(p, ilen, mod), flint_printf("\n\n");
                    _nmod_vec_print(rok->coeffs, rok->length, mod), flint_printf("\n\n");
                    _nmod_vec_print(qok->coeffs, qok->length, mod), flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            /* extend circulant, general */
            {
                ulong c = n_randint(state, prime);
                ulong c_precomp = n_mulmod_precomp_shoup(c, prime);

                _nmod_vec_zero(rec, d+1);
                rec[0] = n_negmod(c, prime);
                rec[d] = 1;
                _nmod_vec_set(pok, p, d);
                extend_recurrence(pok, ilen, rec, d, mod);

                _nmod_vec_set(pp, p, d);
                _nmod_poly_divrem_circulant_lazy_4_2_t(pp, ilen, d, c, c_precomp, prime);

                for (ulong i = d; i < ilen; i++)
                {
                    if (pp[i] >= prime)
                        pp[i] -= prime;
                }

                result = (_nmod_vec_equal(pp, pok, ilen));
                if (!result)
                {
                    flint_printf("FAIL (_nmod_poly_divrem_circulant_lazy_4_2_t):\n");
                    flint_printf("a->length = %wd, ilen = %wu, modulus = %wu\n", a->length, ilen, prime);
                    flint_printf("d = %wu, c = %wu, input vec :\n", d, c);
                    _nmod_vec_print(p, ilen, mod), flint_printf("\n\n");
                    _nmod_vec_print(rok->coeffs, rok->length, mod), flint_printf("\n\n");
                    _nmod_vec_print(qok->coeffs, qok->length, mod), flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            /* extend circulant, 1 */
            {
                ulong c = 1;

                _nmod_vec_zero(rec, d+1);
                rec[0] = n_negmod(c, prime);
                rec[d] = 1;
                _nmod_vec_set(pok, p, d);
                extend_recurrence(pok, ilen, rec, d, mod);

                _nmod_vec_set(pp, p, d);
                _nmod_poly_divrem_circulant1_t(pp, ilen, d);

                result = _nmod_vec_equal(pp, pok, ilen);
                if (!result)
                {
                    flint_printf("FAIL (_nmod_poly_divrem_circulant1_t):\n");
                    flint_printf("a->length = %wd, ilen = %wu, modulus = %wu\n", a->length, ilen, prime);
                    flint_printf("d = %wu, c = %wu, input vec :\n", d, c);
                    _nmod_vec_print(p, ilen, mod), flint_printf("\n\n");
                    _nmod_vec_print(rok->coeffs, rok->length, mod), flint_printf("\n\n");
                    _nmod_vec_print(qok->coeffs, qok->length, mod), flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            /* rem mod prod and extend mod prod work with even d */
            if (d % 2 == 0)
            {
                // number of roots in F->tab_w at least mod_node * 2**mod_depth + d
                ulong mod_depth = depth + n_randint(state, 3);  /* d <= 2**mod_depth */
                mod_depth = FLINT_MIN(mod_depth, MAX_EVAL_DEPTH-1);
                ulong mod_node = n_randint(state, UWORD(1) << (MAX_EVAL_DEPTH-1 - mod_depth));

                /* compute modulus polynomial b == prod((x-w_k) (x+w_k) for */
                /* w_k == F->tab_w[2**depth * node + 2*k] for k in range(d/2)) */
                /* this is also prod(x - w_k), w_k = roots[2**depth * node + k] for k in range(d) */
                nmod_poly_product_roots_nmod_vec(b, roots + (mod_node << mod_depth), d);

                /* rem mod prod */
                nmod_poly_rem(rok, a, b);

                _nmod_vec_set(pp, p, ilen);
                _nmod_poly_rem_prod_root1_lazy_4_4(pp, ilen, d, mod_depth, mod_node, Fargs);
                for (ulong i = 0; i < d; i++)
                {
                    if (pp[i] >= 2*prime)
                        pp[i] -= 2*prime;
                    if (pp[i] >= prime)
                        pp[i] -= prime;
                }

                result = (_nmod_vec_equal(pp, rok->coeffs, rok->length) && _nmod_vec_is_zero(pp+rok->length, d - rok->length));
                if (!result)
                {
                    flint_printf("FAIL (_nmod_poly_rem_prod_root1_lazy_4_4):\n");
                    flint_printf("a->length = %wd, ilen = %wu, modulus = %wu\n", a->length, ilen, prime);
                    flint_printf("d = %wu, mod_depth = %wu, mod_node = %wu, input vec :\n", d, mod_depth, mod_node);
                    _nmod_vec_print(p, ilen, mod), flint_printf("\n\n");
                    _nmod_vec_print(rok->coeffs, d, mod); flint_printf("\n\n");
                    _nmod_vec_print(pp, d, mod); flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }

                /* extend mod prod */
                _nmod_vec_set(pok, p, d);
                extend_recurrence(pok, ilen, b->coeffs, d, mod);

                _nmod_vec_set(pp, p, d);
                _nmod_poly_rem_prod_root1_t_lazy_4_4(pp, ilen, d, mod_depth, mod_node, Fargs);

                for (ulong i = 0; i < ilen; i++)
                {
                    if (pp[i] >= 2*prime)
                        pp[i] -= 2*prime;
                    if (pp[i] >= prime)
                        pp[i] -= prime;
                }

                result = _nmod_vec_equal(pp, pok, ilen);
                if (!result)
                {
                    flint_printf("FAIL (_nmod_poly_rem_prod_root1_t_lazy_4_4):\n");
                    flint_printf("a->length = %wd, ilen = %wu, modulus = %wu\n", a->length, ilen, prime);
                    flint_printf("d = %wu, mod_depth = %wu, mod_node = %wu, input vec :\n", d, mod_depth, mod_node);
                    _nmod_vec_print(p, d, mod), flint_printf("\n\n");
                    _nmod_vec_print(pok, ilen, mod); flint_printf("\n\n");
                    _nmod_vec_print(pp, ilen, mod); flint_printf("\n\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            _nmod_vec_clear(p);
            _nmod_vec_clear(pp);
            _nmod_vec_clear(pok);
            _nmod_vec_clear(rec);
            nmod_poly_clear(a);
            nmod_poly_clear(b);
            nmod_poly_clear(qok);
            nmod_poly_clear(rok);
        }

    }

    TEST_FUNCTION_END(state);
}

#undef MAX_EVAL_DEPTH
