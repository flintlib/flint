/*
    Copyright (C) 2015 FLINT authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_factor.h"

void _fmpz_factor_randtest(fmpz_factor_t f, flint_rand_t state,
        slong num, flint_bitcnt_t bits);
void _fmpz_factor_set(fmpz_factor_t z, const fmpz_factor_t x);
int _fmpz_factor_equal(const fmpz_factor_t x, const fmpz_factor_t y);

void
_fmpz_factor_randtest(fmpz_factor_t f, flint_rand_t state,
        slong num, flint_bitcnt_t bits)
{
    slong i;
    ulong n;
    int s;
    ulong bases_not_zero, exp_not_zero;

    bases_not_zero = n_randint(state, 2);
    exp_not_zero = n_randint(state, 2);

    /* sign is -1 or 1 or rarely 0 */
    s = 0;
    n = n_randint(state, 11);
    if (n)
    {
        s = (n % 2) ? 1 : -1;
    }
    f->sign = s;

    _fmpz_factor_fit_length(f, num);
    _fmpz_factor_set_length(f, num);
    for (i = 0; i < num; i++)
    {
        if (bases_not_zero)
        {
            fmpz_randtest_not_zero(f->p+i, state, bits);
        }
        else
        {
            fmpz_randtest(f->p+i, state, bits);
        }
        if (exp_not_zero)
        {
            f->exp[i] = n_randint(state, 3) + 1;
        }
        else
        {
            f->exp[i] = n_randint(state, 4);
        }
    }
}

void
_fmpz_factor_set(fmpz_factor_t z, const fmpz_factor_t x)
{
    if (z != x)
    {
        slong i, len;
        len = x->num;
        z->sign = x->sign;
        _fmpz_factor_fit_length(z, len);
        _fmpz_factor_set_length(z, len);
        _fmpz_vec_set(z->p, x->p, len);
        for (i = 0; i < len; i++)
        {
            z->exp[i] = x->exp[i];
        }
    }
}

int
_fmpz_factor_equal(const fmpz_factor_t x, const fmpz_factor_t y)
{
    /* this is not equivalence, but equality in the entrywise sense */
    slong len, i;

    if (x->sign != y->sign || x->num != y->num)
    {
        return 0;
    }
    len = x->num;
    if (!_fmpz_vec_equal(x->p, y->p, len))
    {
        return 0;
    }
    for (i = 0; i < len; i++)
    {
        if (x->exp[i] != y->exp[i])
        {
            return 0;
        }
    }
    return 1;
}

TEST_FUNCTION_START(fmpz_factor_refine, state)
{
    int iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        int i;
        fmpz_factor_t f, g;
        slong num;
        flint_bitcnt_t bits;

        bits = n_randint(state, 80) + 2;
        num = n_randint(state, 10);

        /* sample a factor structure that is probably not in canonical form */
        fmpz_factor_init(f);
        _fmpz_factor_randtest(f, state, num, bits);

        /* compute the factor refinement */
        fmpz_factor_init(g);
        fmpz_factor_refine(g, f);

        /* each base must not be less than 2 */
        for (i = 0; i < g->num; i++)
        {
            if (fmpz_cmp_ui(g->p+i, 2) < 0)
            {
                flint_printf("FAIL (base minimum)\n");
                fmpz_print(g->p+i); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        /* bases must be increasing */
        for (i = 0; i < g->num-1; i++)
        {
            if (fmpz_cmp(g->p+i, g->p+i+1) >= 0)
            {
                flint_printf("FAIL (base sorting)\n");
                fmpz_print(g->p+i); flint_printf(" ");
                fmpz_print(g->p+i+1); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        /* each exponent must not be less than 1 */
        for (i = 0; i < g->num; i++)
        {
            if (g->exp[i] < 1)
            {
                flint_printf("FAIL (exponent minimum)\n");
                flint_printf("%wd\n", g->exp[i]);
                fflush(stdout);
                flint_abort();
            }
        }

        /* bases must be coprime */
        {
            slong u, v;
            fmpz_t x;
            fmpz_init(x);
            for (u = 0; u < g->num; u++)
            {
                for (v = 0; v < u; v++)
                {
                    fmpz_gcd(x, g->p+u, g->p+v);
                    if (!fmpz_is_one(x))
                    {
                        flint_printf("FAIL (coprime bases)\n");
                        fmpz_print(g->p+u); flint_printf(" ");
                        fmpz_print(g->p+v); flint_printf("\n");
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }
            fmpz_clear(x);
        }

        /* pre- and post- refinement products must be equal */
        {
            fmpz_t x, y;

            fmpz_init(x);
            fmpz_init(y);

            fmpz_factor_expand(x, f);
            fmpz_factor_expand(y, g);

            if (!fmpz_equal(x, y))
            {
                flint_printf("FAIL (product equality)\n");
                fmpz_factor_print(f); flint_printf(" : ");
                fmpz_print(x); flint_printf("\n");
                fmpz_factor_print(g); flint_printf(" : ");
                fmpz_print(y); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_clear(x);
            fmpz_clear(y);
        }

        /*
         * If the product is not zero, each pre-refinement base with a positive
         * exponent must be representable as a product of powers of
         * post-refinement bases.
         */
        if (g->sign)
        {
            slong i, j;
            fmpz_t x;
            fmpz_init(x);
            for (i = 0; i < f->num; i++)
            {
                if (f->exp[i] != WORD(0))
                {
                    fmpz_abs(x, f->p+i);
                    for (j = 0; j < g->num; j++)
                    {
                        while (fmpz_divisible(x, g->p+j))
                        {
                            fmpz_divexact(x, x, g->p+j);
                        }
                    }
                    if (!fmpz_is_one(x))
                    {
                        flint_printf("FAIL (representation)\n");
                        fmpz_factor_print(f); flint_printf("\n");
                        fmpz_factor_print(g); flint_printf("\n");
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }
            fmpz_clear(x);
        }

        /* check idempotence */
        {
            fmpz_factor_t h;
            fmpz_factor_init(h);
            fmpz_factor_refine(h, g);
            if (!_fmpz_factor_equal(g, h))
            {
                flint_printf("FAIL (idempotence)\n");
                fmpz_factor_print(f); flint_printf("\n");
                fmpz_factor_print(g); flint_printf("\n");
                fmpz_factor_print(h); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
            fmpz_factor_clear(h);
        }

        /* check aliasing */
        {
            fmpz_factor_t h;
            fmpz_factor_init(h);
            _fmpz_factor_set(h, f);
            fmpz_factor_refine(h, h);
            if (!_fmpz_factor_equal(g, h))
            {
                flint_printf("FAIL (aliasing)\n");
                fmpz_factor_print(f); flint_printf("\n");
                fmpz_factor_print(g); flint_printf("\n");
                fmpz_factor_print(h); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
            fmpz_factor_clear(h);
        }

        fmpz_factor_clear(f);
        fmpz_factor_clear(g);
    }

    /*
     * check that factor refinement of (1, 2, 3, ..., n) gives the same
     * result as factoring factorial(n).
     */
    {
        slong i, j;
        for (i = 0; i < 100; i++)
        {
            fmpz_t x, z;
            mpz_t y;
            fmpz_factor_t f, g, h;

            fmpz_init(x);
            fmpz_init(z);
            mpz_init(y);
            fmpz_factor_init(f);
            fmpz_factor_init(g);
            fmpz_factor_init(h);

            flint_mpz_fac_ui(y, i);
            fmpz_set_mpz(x, y);
            f->sign = 1;
            _fmpz_factor_append(f, x, 1);

            g->sign = 1;
            for (j = 1; j < i+1; j++)
            {
                _fmpz_factor_append_ui(g, j, 1);
            }

            fmpz_factor(f, x);
            fmpz_factor_refine(h, g);

            if (!_fmpz_factor_equal(f, h))
            {
                flint_printf("FAIL (factorial)\n");
                flint_printf("%wd\n", i);
                fflush(stdout);
                flint_abort();
            }

            fmpz_clear(x);
            fmpz_clear(z);
            mpz_clear(y);
            fmpz_factor_clear(f);
            fmpz_factor_clear(g);
            fmpz_factor_clear(h);
        }
    }

    TEST_FUNCTION_END(state);
}
