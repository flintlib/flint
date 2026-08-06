/*
    Copyright (C) 2019 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpq_vec.h"

/* check that a (= approximation of polynomial root) is in between */
/* c 2^k and (c+1) 2^k (c and k are produced by root isolation).   */
static int check_isolation(fmpz * c, slong k, fmpq_t a)
{
    fmpq_t r1, r2;

    fmpq_init(r1);
    fmpq_init(r2);

    fmpz_set(fmpq_numref(r1), c);
    fmpz_one(fmpq_denref(r1));
    fmpq_set(r2, r1);
    fmpq_add_si(r2, r2, 1);
    if (k > 0)
    {
        fmpq_mul_2exp(r1, r1, (ulong)k);
        fmpq_mul_2exp(r2, r2, (ulong)k);
    }
    else if (k < 0)
    {
        fmpq_div_2exp(r1, r1, (ulong)-k);
        fmpq_div_2exp(r2, r2, (ulong)-k);
    }
    if (fmpq_cmp(r1, a) >= 0 || fmpq_cmp(r2, a) <= 0)
    {
        fprintf(stderr, "a  = "); fmpq_fprint(stderr, a); fprintf(stderr, "\n");
        fprintf(stderr, "r1 = "); fmpq_fprint(stderr, r1); fprintf(stderr, "\n");
        fprintf(stderr, "r2 = "); fmpq_fprint(stderr, r2); fprintf(stderr, "\n");
        return 1;
    }

    fmpq_clear(r1);
    fmpq_clear(r2);

    return 0;
}

static void _slong_vec_print(const slong * vec, slong len)
{
    slong i;
    flint_printf("%wd", len);
    for (i = 0; i < len; ++i)
        flint_printf(" %wd", vec[i]);
}


static void check_intervals(
      fmpq * vec, slong len,
      fmpq * exact, slong n_exact,
      fmpz * c_array, slong * k_array, slong n_interval)
{
    slong i,j,k;
    fmpq_t x,y;

    if (n_exact + n_interval != len)
    {
        flint_printf("ERROR:\n");
        flint_printf("found n_exact = %wd and n_interval = %wd but n = %wd\n", n_exact, n_interval, len);
        flint_printf("vec = "); _fmpq_vec_print(vec, len);
        flint_printf("\n");
        flint_abort();
    }

    fmpq_init(x);
    fmpq_init(y);
    i = j = k = 0;
    while ((j < n_exact) || (k < n_interval))
    {
        if (k < n_interval)
        {
            fmpz_set(fmpq_numref(x), c_array + k);
            fmpz_set(fmpq_numref(y), c_array + k);
            fmpz_one(fmpq_denref(x));
            fmpz_one(fmpq_denref(y));
            fmpz_add_ui(fmpq_numref(y), fmpq_numref(y), 1);
            if (k_array[k] > 0)
            {
                fmpq_mul_2exp(x, x, (ulong)k_array[k]);
                fmpq_mul_2exp(y, y, (ulong)k_array[k]);
            }
            else if (k_array[k] < 0)
            {
                fmpq_div_2exp(x, x, (ulong)-k_array[k]);
                fmpq_div_2exp(y, y, (ulong)-k_array[k]);
            }
        }


        if ((j < n_exact) && (k < n_interval))
        {
            if ((fmpq_cmp(exact + j, x) > 0) && (fmpq_cmp(y, exact + j) > 0))
            {
                flint_printf("ERROR:\n");
                flint_printf("the %wd-th exact root ", j);
                fmpq_print(exact + j); flint_printf("\n");
                flint_printf("belongs to the %wd-th interval ", k);
                flint_printf("("); fmpq_print(x); flint_printf(","); fmpq_print(y); flint_printf(")");
                flint_printf("\n");
                flint_printf("vec   = "); _fmpq_vec_print(vec, len);
                flint_printf("\nc     = "); _fmpz_vec_print(c_array, n_interval);
                flint_printf("\nk     = "); _slong_vec_print(k_array, n_interval);
                flint_printf("\nexact = "); _fmpq_vec_print(exact, n_exact);
                flint_printf("\n");
                flint_abort();
            }
        }

        if ((j == n_exact) || ((k < n_interval) && (fmpq_cmp(exact + j, y) >= 0)))
        {
            /* check interval */
            if ((fmpq_cmp(x, vec + i) > 0) ||
                 fmpq_cmp(vec + i, y) > 0)
            {
                flint_printf("ERROR:\n");
                flint_printf(" the %wd-th root is vec[%wd] = ", i, i);
                fmpq_print(vec + i);
                flint_printf(" but got %wd-th interval ", k);
                flint_printf("("); fmpq_print(x); flint_printf(","); fmpq_print(y); flint_printf(")");
                flint_printf("\n");
                flint_printf("vec   = "); _fmpq_vec_print(vec, len);
                flint_printf("\nc     = "); _fmpz_vec_print(c_array, n_interval);
                flint_printf("\nk     = "); _slong_vec_print(k_array, n_interval);
                flint_printf("\nexact = "); _fmpq_vec_print(exact, n_exact);
                flint_printf("\n");
                flint_abort();
            }
            k += 1;
            i += 1;
        }
        else
        {
            /* check root */
            if (!fmpq_equal(vec + i, exact + j))
            {
                flint_printf("ERROR:\n");
                flint_printf("the %wd-th root is vec[%wd] = ", i, i);
                fmpq_print(vec + i);
                flint_printf(" but got %wd-th exact root ", j);
                fmpq_print(exact + j);
                flint_printf("\n");
                flint_printf("vec   = "); _fmpq_vec_print(vec, len);
                flint_printf("\nc     = "); _fmpz_vec_print(c_array, n_interval);
                flint_printf("\nk     = "); _slong_vec_print(k_array, n_interval);
                flint_printf("\nexact = "); _fmpq_vec_print(exact, n_exact);
                flint_printf("\n");
                flint_abort();
            }
            j += 1;
            i += 1;
        }
    }

    fmpq_clear(x);
    fmpq_clear(y);
}

static void fmpz_poly_from_fmpq_roots(fmpz_poly_t p, const fmpq * vec, slong n)
{
    fmpz_poly_t q;
    slong i;

    fmpz_poly_init(q);
    fmpz_poly_one(p);
    for (i = 0; i < n; i++)
    {
        if (fmpq_is_zero(vec + i))
        {
            fmpz_poly_set_coeff_si(q, 0, 0);
            fmpz_poly_set_coeff_si(q, 1, 1);
        }
        else
        {
            fmpz_poly_set_coeff_fmpz(q, 0, fmpq_numref(vec + i));
            fmpz_neg(fmpq_poly_numref(q), fmpq_poly_numref(q));
            fmpz_poly_set_coeff_fmpz(q, 1, fmpq_denref(vec + i));
        }
        fmpz_poly_mul(p, p, q);
    }
    fmpz_poly_clear(q);
}

TEST_FUNCTION_START(fmpz_poly_isolate_real_roots, state)
{
    {
        fmpz_poly_t p;
        fmpq * exact_roots;
        fmpz * c_array;
        slong * k_array;
        slong n_exact, n_interval;
        fmpq_t a;

        fmpq_init(a);
        fmpz_poly_init(p);
        exact_roots = _fmpq_vec_init(5);
        c_array = _fmpz_vec_init(5);
        k_array = (slong *) flint_malloc(5 * sizeof(slong));

        /* -1705*x^2 - 7650*x - 3297 */
        /* roots: -4.0038 and -0.4829 */
        fmpz_poly_set_coeff_si(p, 0, -3297);
        fmpz_poly_set_coeff_si(p, 1, -7650);
        fmpz_poly_set_coeff_si(p, 2, -1705);

        fmpz_poly_isolate_real_roots(exact_roots, &n_exact,
                c_array, k_array, &n_interval, p);

        if (n_exact != 0 || n_interval != 2)
        {
            fprintf(stderr, "wrong number of isolated roots\n");
            flint_abort();
        }

        fmpq_set_si(a, WORD(-1933157935), WORD(482826508)); /* approx of root 1 */
        if (check_isolation(c_array, k_array[0], a))
        {
            fprintf(stderr, "Failed root1 of poly1\n");
            flint_abort();
        }

        fmpq_set_si(a, WORD(-151354505), WORD(313384144)); /* approx of root 2 */
        if (check_isolation(c_array + 1, k_array[1], a))
        {
            fprintf(stderr, "Failed root2 of poly1\n");
            flint_abort();
        }

        /* x^2 - 7650*x - 13297 */
        /* roots: -1.7377 and 7651.7377 */
        fmpz_poly_set_coeff_si(p, 0, -13297);
        fmpz_poly_set_coeff_si(p, 1, -7650);
        fmpz_poly_set_coeff_si(p, 2, 1);

        n_exact = n_interval = 0;
        fmpz_poly_isolate_real_roots(exact_roots, &n_exact,
                c_array, k_array, &n_interval, p);

        if (n_exact != 0 || n_interval != 2)
        {
            fprintf(stderr, "wrong number of isolated roots\n");
            flint_abort();
        }

        fmpq_set_si(a, WORD(-833585025), WORD(479685194));
        if (check_isolation(c_array, k_array[0], a))
        {
            fprintf(stderr, "Failed root1 of poly2\n");
            flint_abort();
        }
        fmpq_set_si(a, WORD(1283601967), WORD(167753));
        if (check_isolation(c_array + 1, k_array[1], a))
        {
            fprintf(stderr, "Failed root2 of poly2\n");
            flint_abort();
        }

        /* x^2 - 1505*x + 566255 */
        /* roots: 751.381, 753.618 */
        fmpz_poly_set_coeff_si(p, 0, 566255);
        fmpz_poly_set_coeff_si(p, 1, -1505);
        fmpz_poly_set_coeff_si(p, 2, 1);

        n_exact = n_interval = 0;
        fmpz_poly_isolate_real_roots(exact_roots, &n_exact,
                c_array, k_array, &n_interval, p);

        if (n_exact != 0 || n_interval != 2)
        {
            fprintf(stderr, "wrong number of isolated roots\n");
            flint_abort();
        }

        fmpq_set_si(a, WORD(1011562248), WORD(1346269));
        if (check_isolation(c_array, k_array[0], a))
        {
            fprintf(stderr, "Failed root1 of poly2\n");
            flint_abort();
        }
        fmpq_set_si(a, WORD(148024147), WORD(196418));
        if (check_isolation(c_array + 1, k_array[1], a))
        {
            fprintf(stderr, "Failed root2 of poly2\n");
            flint_abort();
        }

        /* 146434129 * x^2 - 134751 * x + 31 */
        /* roots: 0.0004601002 and 0.0004601155 */
        fmpz_poly_set_coeff_si(p, 0, 31);
        fmpz_poly_set_coeff_si(p, 1, -134751);
        fmpz_poly_set_coeff_si(p, 2, 146434129);

        n_exact = n_interval = 0;
        fmpz_poly_isolate_real_roots(exact_roots, &n_exact,
                c_array, k_array, &n_interval, p);

        if (n_exact != 0 || n_interval != 2)
        {
            fprintf(stderr, "wrong number of isolated roots\n");
            flint_abort();
        }

        fmpq_set_si(a, WORD(381890), WORD(830014731));
        if (check_isolation(c_array, k_array[0], a))
        {
            fprintf(stderr, "Failed root1 of poly2\n");
            flint_abort();
        }
        fmpq_set_si(a, WORD(456915), WORD(993044056));
        if (check_isolation(c_array + 1, k_array[1], a))
        {
            fprintf(stderr, "Failed root2 of poly2\n");
            flint_abort();
        }

        _fmpq_vec_clear(exact_roots, 5);
        _fmpz_vec_clear(c_array, 5);
        flint_free(k_array);
        fmpq_clear(a);
        fmpz_poly_clear(p);
    }

    {
        int iter;

        for (iter = 0; iter < 500; iter++)
        {
            fmpq vec[30];
            fmpz c_array[30];
            slong k_array[30];
            fmpq exact_array[30];
            fmpz_poly_t p,q;

            slong n = (slong)n_randint(state, 30);      /* real roots            */
            slong nc = 1 + (slong)n_randint(state, 30); /* complex roots */
            slong i;
            slong n_exact, n_interval;

            if (n + nc == 0) continue;

            for(i = 0; i < 30; ++i)
            {
                fmpq_init(vec + i);
                fmpz_init(c_array + i);
                fmpq_init(exact_array + i);
            }

            _fmpq_vec_randtest_uniq_sorted(vec, state, n, 30);

            fmpz_poly_init(p);
            fmpz_poly_from_fmpq_roots(p, vec, n);

            fmpz_poly_init(q);
            fmpz_poly_randtest_no_real_root(q, state, nc, 100);
            fmpz_poly_mul(p, p, q);

            fmpz_poly_isolate_real_roots(exact_array, &n_exact, c_array, k_array, &n_interval, p);

            check_intervals(vec, n, exact_array, n_exact, c_array, k_array, n_interval);

            fmpz_poly_clear(p);
            fmpz_poly_clear(q);
            for(i = 0; i < 30; ++i)
            {
                fmpq_clear(vec + i);
                fmpz_clear(c_array + i);
                fmpq_clear(exact_array + i);
            }
        }
    }

    TEST_FUNCTION_END(state);
}

