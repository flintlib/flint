/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"
#include "ulong_extras.h"
#include "fmpz_poly/impl.h"

/*
    Tests whether the polynomial is suitably normalised for the
    result of a GCD operation, that is, whether it's leading
    coefficient is non-negative.
 */
#ifndef _t_gcd_is_canonical
#define _t_gcd_is_canonical _t_gcd_is_canonical
static
int _t_gcd_is_canonical(const fmpz_poly_t poly)
{
    return fmpz_poly_is_zero(poly) || (fmpz_sgn(fmpz_poly_lead(poly)) > 0);
}
#endif

TEST_FUNCTION_START(fmpz_poly_gcd_modular, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(c, state, n_randint(state, 40), 80);

        fmpz_poly_gcd_modular(a, b, c);
        fmpz_poly_gcd_modular(b, b, c);

        result = (fmpz_poly_equal(a, b) && _t_gcd_is_canonical(a));
        if (!result)
        {
            flint_printf("FAIL (aliasing a and b):\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(c, state, n_randint(state, 40), 80);

        fmpz_poly_gcd_modular(a, b, c);
        fmpz_poly_gcd_modular(c, b, c);

        result = (fmpz_poly_equal(a, c) && _t_gcd_is_canonical(a));
        if (!result)
        {
            flint_printf("FAIL (aliasing a and c):\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("c = "), fmpz_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check that a divides GCD(af, ag) */
    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, d, f, g, q, r;

        fmpz_poly_init(a);
        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_randtest_not_zero(a, state, n_randint(state, 100) + 1, 40);
        fmpz_poly_randtest(f, state, n_randint(state, 100), 40);
        fmpz_poly_randtest(g, state, n_randint(state, 100), 40);

        fmpz_poly_mul(f, a, f);
        fmpz_poly_mul(g, a, g);
        fmpz_poly_gcd_modular(d, f, g);

        fmpz_poly_divrem_divconquer(q, r, d, a);

        result = fmpz_poly_is_zero(r) && _t_gcd_is_canonical(d);
        if (!result)
        {
           flint_printf("FAIL (check a | gcd(af, ag)):\n");
           flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n");
           flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n");
           flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n");
           flint_printf("d = "), fmpz_poly_print(d), flint_printf("\n");
           fflush(stdout);
           flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
    }

    /* Check that a == GCD(af, ag) when GCD(f, g) = 1 */
    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, d, f, g, q, r;

        fmpz_poly_init(a);
        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_randtest_not_zero(a, state, n_randint(state, 100) + 1, 200);
        do {
           fmpz_poly_randtest(f, state, n_randint(state, 100), 200);
           fmpz_poly_randtest(g, state, n_randint(state, 100), 200);
           fmpz_poly_gcd_heuristic(d, f, g);
        } while (!(d->length == 1 && fmpz_is_one(d->coeffs)));

        fmpz_poly_mul(f, a, f);
        fmpz_poly_mul(g, a, g);
        fmpz_poly_gcd_modular(d, f, g);

        if (!_t_gcd_is_canonical(a)) fmpz_poly_neg(a, a);

        result = fmpz_poly_equal(d, a) && _t_gcd_is_canonical(d);
        if (!result)
        {
           flint_printf("FAIL (check a == gcd(af, ag) when gcd(f, g) = 1):\n");
           flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n");
           flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n");
           flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n");
           flint_printf("d = "), fmpz_poly_print(d), flint_printf("\n");
           fflush(stdout);
           flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
    }

    /* Sebastian's test case */
    {
        fmpz_poly_t a, b, d;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(d);

        fmpz_poly_set_coeff_ui(b, 2, 1);
        fmpz_poly_set_coeff_si(a, 0, -32);
        fmpz_poly_set_coeff_si(a, 1, 24);

        fmpz_poly_gcd_modular(d, a, b);

        result = (d->length == 1 && fmpz_is_one(d->coeffs));
        if (!result)
        {
            flint_printf("FAIL (check 1 == gcd(x^2, 24*x - 32):\n");
            fmpz_poly_print(d); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(d);
    }

    /* another test case */
    {
        fmpz_poly_t a, b, d, e;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(d);
        fmpz_poly_init(e);

        fmpz_poly_set_str(a, "12  0 0 0 0 0 0 0 0 0 8582594367 -9297159048333985579007 33822867456");
        fmpz_poly_set_str(b, "8  0 0 -258272396248218664896 0 -2762 -549690802047 -3771028 8796059467776");
        fmpz_poly_set_str(e, "3  0 0 1");

        fmpz_poly_gcd_modular(d, a, b);

        result = fmpz_poly_equal(d, e);
        if (!result)
        {
            flint_printf("FAIL (check special #2):\n");
            fmpz_poly_print(d); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(d);
        fmpz_poly_clear(e);
    }

    /* Compare with the subresultant algorithm, using small primes to
       exercise unlucky and bad primes */
    for (i = 0; i < 2000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, f, g, d, e;
        ulong first_prime;
        slong len;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(d);
        fmpz_poly_init(e);

        len = 1 + n_randint(state, 20);

        fmpz_poly_randtest_not_zero(a, state, 1 + n_randint(state, len), 1 + n_randint(state, 100));
        fmpz_poly_randtest(f, state, n_randint(state, len), 1 + n_randint(state, 100));
        fmpz_poly_randtest(g, state, n_randint(state, len), 1 + n_randint(state, 100));

        if (n_randint(state, 2))
        {
            fmpz_poly_mul(f, a, f);
            fmpz_poly_mul(g, a, g);
        }

        /* structured leading coefficients */
        if (n_randint(state, 4) == 0 && f->length > 0 && g->length > 0)
        {
            fmpz_t c;
            fmpz_init(c);
            fmpz_randtest_not_zero(c, state, 1 + n_randint(state, 100));
            fmpz_mul(f->coeffs + f->length - 1, f->coeffs + f->length - 1, c);
            fmpz_mul(g->coeffs + g->length - 1, g->coeffs + g->length - 1, c);
            fmpz_clear(c);
        }

        /* structured constant coefficients */
        if (n_randint(state, 4) == 0 && f->length > 0 && g->length > 0)
        {
            fmpz_t c;
            fmpz_init(c);
            fmpz_randtest_not_zero(c, state, 1 + n_randint(state, 100));
            fmpz_mul(f->coeffs, f->coeffs, c);
            fmpz_mul(g->coeffs, g->coeffs, c);
            _fmpz_poly_normalise(f);
            _fmpz_poly_normalise(g);
            fmpz_clear(c);
        }

        if (n_randint(state, 8) == 0)
            fmpz_poly_mul(f, f, g);

        switch (n_randint(state, 4))
        {
            case 0: first_prime = 0; break;
            case 1: first_prime = 1; break;
            case 2: first_prime = n_randint(state, 100); break;
            default: first_prime = n_randtest(state) >> (1 + n_randint(state, FLINT_BITS - 1));
        }

        if (f->length < g->length)
            fmpz_poly_swap(f, g);

        if (g->length == 0)
        {
            fmpz_poly_gcd(e, f, g);
            fmpz_poly_set(d, e);
        }
        else
        {
            fmpz_poly_fit_length(d, g->length);
            _fmpz_poly_gcd_modular_primes(d->coeffs, f->coeffs, f->length, g->coeffs, g->length, first_prime);
            _fmpz_poly_set_length(d, g->length);
            _fmpz_poly_normalise(d);
            fmpz_poly_gcd_subresultant(e, f, g);
        }

        if (!fmpz_poly_equal(d, e))
        {
            flint_printf("FAIL (comparison with subresultant):\n");
            flint_printf("first_prime = %wu\n", first_prime);
            flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n");
            flint_printf("d = "), fmpz_poly_print(d), flint_printf("\n");
            flint_printf("e = "), fmpz_poly_print(e), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(d);
        fmpz_poly_clear(e);
    }

    /* Large input (exercising many primes and FFT primes): check that
       gcd(a u, a v) = a when u and v are coprime modulo some prime */
    for (i = 0; i < 4 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, u, v, f, g, d;
        nmod_poly_t um, vm, gm;
        slong lena, lenu, lenv, bits, ubits;
        ulong p;

        fmpz_poly_init(a);
        fmpz_poly_init(u);
        fmpz_poly_init(v);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(d);

        if (i % 2 == 0)
        {
            /* long polynomials (using FFT primes); sometimes with a large
               enough gcd to need more than the FFT primes of the default
               fft_small context */
            lena = 1000 + n_randint(state, 2000);
            lenu = 1500 + n_randint(state, 1000);
            lenv = 1500 + n_randint(state, 1000);
            if (i % 4 == 0)
                bits = 1 + n_randint(state, 100);
            else
                bits = 500 + n_randint(state, 100);
            ubits = 1 + n_randint(state, 100);
        }
        else if (i % 4 == 1)
        {
            /* large coefficients */
            lena = 1 + n_randint(state, 10);
            lenu = 1 + n_randint(state, 10);
            lenv = 1 + n_randint(state, 10);
            bits = 4000 + n_randint(state, 2000);
            ubits = 1 + n_randint(state, 100);
        }
        else
        {
            /* the quotient a u / a has much larger coefficients than the
               gcd (so it is not computed multimodularly) */
            lena = 10 + n_randint(state, 10);
            lenu = 10 + n_randint(state, 10);
            lenv = 10 + n_randint(state, 10);
            bits = 1 + n_randint(state, 20);
            ubits = 4000 + n_randint(state, 2000);
        }

        fmpz_poly_randtest_not_zero(a, state, lena, bits);
        fmpz_poly_randtest_not_zero(u, state, lenu, ubits);
        fmpz_poly_randtest_not_zero(v, state, lenv, 1 + n_randint(state, 100));

        /* make the lengths exact */
        fmpz_poly_set_coeff_ui(a, lena - 1, 1 + n_randint(state, 100));
        fmpz_poly_set_coeff_ui(u, lenu - 1, 1 + n_randint(state, 100));
        fmpz_poly_set_coeff_ui(v, lenv - 1, 1 + n_randint(state, 100));
        fmpz_poly_set_coeff_si(a, 0, -(slong) n_randint(state, 1000) - 1);
        fmpz_randbits(a->coeffs + a->length / 2, state, bits);

        p = n_randprime(state, FLINT_BITS - 1, 0);
        nmod_poly_init(um, p);
        nmod_poly_init(vm, p);
        nmod_poly_init(gm, p);
        fmpz_poly_get_nmod_poly(um, u);
        fmpz_poly_get_nmod_poly(vm, v);
        nmod_poly_gcd(gm, um, vm);

        if (um->length == u->length && vm->length == v->length && gm->length == 1)
        {
            fmpz_poly_mul(f, a, u);
            fmpz_poly_mul(g, a, v);
            fmpz_poly_gcd_modular(d, f, g);

            /* the expected result is a times the gcd of the contents of u, v,
               with positive leading coefficient */
            {
                fmpz_t cu, cv;
                fmpz_init(cu);
                fmpz_init(cv);
                fmpz_poly_content(cu, u);
                fmpz_poly_content(cv, v);
                fmpz_gcd(cu, cu, cv);
                fmpz_poly_scalar_mul_fmpz(a, a, cu);
                if (fmpz_sgn(fmpz_poly_lead(a)) < 0)
                    fmpz_poly_neg(a, a);
                fmpz_clear(cu);
                fmpz_clear(cv);
            }

            if (!fmpz_poly_equal(d, a))
            {
                flint_printf("FAIL (large input):\n");
                flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n");
                flint_printf("d = "), fmpz_poly_print(d), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_poly_clear(um);
        nmod_poly_clear(vm);
        nmod_poly_clear(gm);
        fmpz_poly_clear(a);
        fmpz_poly_clear(u);
        fmpz_poly_clear(v);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
