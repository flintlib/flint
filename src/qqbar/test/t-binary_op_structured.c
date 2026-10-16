/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "qqbar.h"

/* Test arithmetic on inputs with structure (elements of a common cyclotomic
   field, conjugates) which is exploited by qqbar_binary_op. */

static void
_rand_cyclotomic(qqbar_t x, flint_rand_t state, slong n, int real)
{
    qqbar_t t;
    slong k, terms, e, c;

    qqbar_init(t);

    do
    {
        qqbar_zero(x);
        terms = 1 + n_randint(state, 3);

        for (k = 0; k < terms; k++)
        {
            e = n_randint(state, n);
            c = (slong) n_randint(state, 5) - 2;
            if (c == 0)
                c = 1;

            if (real)
                qqbar_cos_pi(t, 2 * e, n);
            else
                qqbar_root_of_unity(t, e, n);

            qqbar_mul_si(t, t, c);
            qqbar_add(x, x, t);
        }
    }
    while (qqbar_degree(x) < 2);

    qqbar_clear(t);
}

static int
_check_result(const qqbar_t res, const acb_t expected)
{
    acb_t z;
    fmpz_poly_factor_t fac;
    int ok;

    acb_init(z);
    fmpz_poly_factor_init(fac);

    qqbar_get_acb(z, res, 128);
    ok = acb_overlaps(z, expected);

    fmpz_poly_factor(fac, QQBAR_POLY(res));
    ok = ok && (fac->num == 1) && (fac->exp[0] == 1) && fmpz_is_one(&fac->c);

    acb_clear(z);
    fmpz_poly_factor_clear(fac);

    return ok;
}

/* random element of a tower of quadratic extensions:
   sqrt(a + b sqrt(c)) + e sqrt(c) + f */
static void
_rand_tower(qqbar_t x, flint_rand_t state)
{
    qqbar_t s, t;
    qqbar_init(s);
    qqbar_init(t);

    do
    {
        qqbar_set_ui(s, 2 + n_randint(state, 5));
        qqbar_sqrt(s, s);
        qqbar_mul_si(t, s, (slong) n_randint(state, 5) - 2);
        qqbar_add_si(t, t, 3 + n_randint(state, 5));
        qqbar_sqrt(x, t);
        if (n_randint(state, 2))
        {
            qqbar_mul_si(t, s, (slong) n_randint(state, 5) - 2);
            qqbar_add(x, x, t);
        }
        qqbar_add_si(x, x, (slong) n_randint(state, 5) - 2);
    }
    while (qqbar_degree(x) < 2);

    qqbar_clear(s);
    qqbar_clear(t);
}

/* random quadratic a + b sqrt(c) */
static void
_rand_quadratic(qqbar_t x, flint_rand_t state)
{
    do
    {
        qqbar_set_si(x, (slong) n_randint(state, 30) - 10);
        qqbar_sqrt(x, x);
        qqbar_mul_si(x, x, 1 + n_randint(state, 3));
        qqbar_add_si(x, x, (slong) n_randint(state, 5) - 2);
    }
    while (qqbar_degree(x) != 2);
}

TEST_FUNCTION_START(qqbar_binary_op_structured, state)
{
    slong iter;
    const slong ns[] = { 13, 16, 20, 21, 24, 28, 36, 44 };

    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y, z, w;
        acb_t a, b, c;
        slong n, i;
        int op, real, kind;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);
        qqbar_init(w);
        acb_init(a);
        acb_init(b);
        acb_init(c);

        n = ns[n_randint(state, sizeof(ns) / sizeof(slong))];
        real = n_randint(state, 2);
        kind = n_randint(state, 6);
        op = n_randint(state, 4);

        _rand_cyclotomic(x, state, n, real);

        if (kind == 3)
        {
            /* quadratic operand */
            if (n_randint(state, 2))
                _rand_tower(x, state);
            _rand_quadratic(y, state);
        }
        else if (kind == 4)
        {
            /* towers of quadratic extensions */
            _rand_tower(x, state);
            _rand_tower(y, state);
            if (n_randint(state, 2))
            {
                qqbar_mul(x, x, y);
                _rand_tower(y, state);
            }
        }
        else if (kind == 5)
        {
            /* two quadratics */
            _rand_quadratic(x, state);
            _rand_quadratic(y, state);
        }
        else if (kind == 0)
        {
            /* a conjugate of x */
            slong d = qqbar_degree(x);
            qqbar_ptr v = _qqbar_vec_init(d);
            qqbar_conjugates(v, x);
            qqbar_set(y, v + n_randint(state, d));
            _qqbar_vec_clear(v, d);
        }
        else if (kind == 1)
        {
            /* element of the same field */
            _rand_cyclotomic(y, state, n, real);
        }
        else
        {
            /* element of a subfield */
            slong m = (n % 4 == 0) ? n / 4 : ((n % 2 == 0) ? n / 2 : n);
            /* the real subfield of Q(zeta_m) must not be Q */
            if (m < 7)
                m = n;
            _rand_cyclotomic(y, state, m, n_randint(state, 2));
        }

        if (n_randint(state, 2))
            qqbar_swap(x, y);

        if (op == 3 && qqbar_is_zero(y))
            op = 2;

        qqbar_get_acb(a, x, 128);
        qqbar_get_acb(b, y, 128);

        if (op == 0)
        {
            qqbar_add(z, x, y);
            acb_add(c, a, b, 128);
            qqbar_sub(w, z, y);
        }
        else if (op == 1)
        {
            qqbar_sub(z, x, y);
            acb_sub(c, a, b, 128);
            qqbar_add(w, z, y);
        }
        else if (op == 2)
        {
            qqbar_mul(z, x, y);
            acb_mul(c, a, b, 128);
            if (qqbar_is_zero(y))
                qqbar_set(w, x);
            else
                qqbar_div(w, z, y);
        }
        else
        {
            qqbar_div(z, x, y);
            acb_div(c, a, b, 128);
            qqbar_mul(w, z, y);
        }

        if (!_check_result(z, c) || !qqbar_equal(w, x))
        {
            flint_printf("FAIL (binary op)\n");
            flint_printf("op = %d, n = %wd\n", op, n);
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
            flint_printf("w = "); qqbar_print(w); flint_printf("\n\n");
            flint_abort();
        }

        /* unary operations which combine x with its complex conjugate */
        for (i = 0; i < 4; i++)
        {
            if (i == 0)
            {
                qqbar_re(z, x);
                acb_set_arb(c, acb_realref(a));
            }
            else if (i == 1)
            {
                qqbar_im(z, x);
                acb_set_arb(c, acb_imagref(a));
            }
            else if (i == 2)
            {
                qqbar_abs(z, x);
                acb_abs(acb_realref(c), a, 128);
                arb_zero(acb_imagref(c));
            }
            else
            {
                qqbar_abs2(z, x);
                acb_abs(acb_realref(c), a, 128);
                arb_sqr(acb_realref(c), acb_realref(c), 128);
                arb_zero(acb_imagref(c));
            }

            if (!_check_result(z, c))
            {
                flint_printf("FAIL (unary op %wd)\n", i);
                flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
                flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
                flint_abort();
            }
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
        qqbar_clear(w);
        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
    }

    TEST_FUNCTION_END(state);
}
