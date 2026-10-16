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
#include "fmpz.h"
#include "fmpq_vec.h"
#include "fmpz_poly/impl.h"

/* check that (c 2^k, (c+1) 2^k) contains exactly one root of p */
static int
check_one_root(const fmpz_poly_t p, const fmpz_t c, slong k)
{
    fmpz_poly_t r;
    fmpz_t one;
    slong n;

    fmpz_poly_init(r);
    fmpz_init(one);
    fmpz_one(one);
    fmpz_poly_set(r, p);
    /* r(x) = p(2^k x), then r(x + c) */
    _fmpz_poly_scale_2exp(r->coeffs, r->length, k);
    fmpz_poly_taylor_shift(r, r, c);
    n = fmpz_poly_num_real_roots_0_1_vca(r);
    fmpz_poly_clear(r);
    fmpz_clear(one);
    return n == 1;
}

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

    /* Sturm isolation with rational roots */
    {
        int iter;

        for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
        {
            fmpq vec[30];
            fmpz c_array[30];
            slong k_array[30];
            fmpq exact_array[30];
            fmpz_poly_t p, q;
            /* small, since the Sturm sequence grows quickly for such input */
            slong n = (slong) n_randint(state, 16);
            slong nc = 1 + (slong) n_randint(state, 8);
            slong i, n_exact, n_interval;

            for (i = 0; i < 30; ++i)
            {
                fmpq_init(vec + i);
                fmpz_init(c_array + i);
                fmpq_init(exact_array + i);
            }

            _fmpq_vec_randtest_uniq_sorted(vec, state, n, 16);
            fmpz_poly_init(p);
            fmpz_poly_from_fmpq_roots(p, vec, n);
            fmpz_poly_init(q);
            fmpz_poly_randtest_no_real_root(q, state, nc, 20);
            fmpz_poly_mul(p, p, q);

            if (_fmpz_poly_isolate_real_roots_sturm(exact_array, &n_exact,
                c_array, k_array, &n_interval, p->coeffs, p->length, 0, WORD_MAX))
            {
                check_intervals(vec, n, exact_array, n_exact, c_array, k_array, n_interval);
            }

            else if (fmpz_poly_is_squarefree(p))
            {
                flint_printf("FAIL: Sturm isolation failed without size bound\n");
                flint_printf("p = %{fmpz_poly}\n", p);
                flint_abort();
            }

            /* the sign change method must either fail or be correct
               (a small budget suffices to exercise both outcomes) */
            if (_fmpz_poly_isolate_real_roots_signs(exact_array, &n_exact,
                c_array, k_array, &n_interval, p->coeffs, p->length, 0, 100))
            {
                check_intervals(vec, n, exact_array, n_exact, c_array, k_array, n_interval);
            }

            fmpz_poly_clear(p);
            fmpz_poly_clear(q);
            for (i = 0; i < 30; ++i)
            {
                fmpq_clear(vec + i);
                fmpz_clear(c_array + i);
                fmpq_clear(exact_array + i);
            }
        }
    }

    /* sign change method for real-rooted polynomials (products of linear
       factors), which must succeed given a large enough budget; we use
       well-separated roots (the number of evaluations grows like the
       inverse of the root separation) */
    {
        int iter;

        for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
        {
            fmpq vec[30];
            fmpz c_array[30];
            slong k_array[30];
            fmpq exact_array[30];
            fmpz_poly_t p;
            slong n = (slong) n_randint(state, 30);
            slong i, n_exact, n_interval;

            for (i = 0; i < 30; ++i)
            {
                fmpq_init(vec + i);
                fmpz_init(c_array + i);
                fmpq_init(exact_array + i);
            }

            /* roots +/- m 2^e with 1 <= m <= 20, -20 <= e <= 20 */
            for (i = 0; i < n; i++)
            {
                fmpq_set_si(vec + i, 1 + n_randint(state, 20), 1);
                if (n_randint(state, 2))
                    fmpq_neg(vec + i, vec + i);
                slong e = (slong) n_randint(state, 41) - 20;
                if (e >= 0)
                    fmpq_mul_2exp(vec + i, vec + i, e);
                else
                    fmpq_div_2exp(vec + i, vec + i, -e);
            }
            _fmpq_vec_sort(vec, n);
            for (i = 0; i + 1 < n; i++)
            {
                if (fmpq_equal(vec + i, vec + i + 1))
                {
                    slong j;
                    for (j = i + 1; j + 1 < n; j++)
                        fmpq_swap(vec + j, vec + j + 1);
                    n--;
                    i--;
                }
            }
            fmpz_poly_init(p);
            fmpz_poly_from_fmpq_roots(p, vec, n);

            if (!_fmpz_poly_isolate_real_roots_signs(exact_array, &n_exact,
                c_array, k_array, &n_interval, p->coeffs, p->length, 0, 1000000))
            {
                flint_printf("FAIL: sign change method failed for a real-rooted polynomial\n");
                flint_printf("p = %{fmpz_poly}\n", p);
                flint_abort();
            }

            check_intervals(vec, n, exact_array, n_exact, c_array, k_array, n_interval);

            fmpz_poly_clear(p);
            for (i = 0; i < 30; ++i)
            {
                fmpq_clear(vec + i);
                fmpz_clear(c_array + i);
                fmpq_clear(exact_array + i);
            }
        }
    }

    /* the default budget for the sign change method must suffice for
       Eulerian polynomials (which are real-rooted, with roots spread over
       many orders of magnitude, and need about len^2 / 4 evaluations) */
    {
        int iter;

        for (iter = 0; iter < 3 * flint_test_multiplier(); iter++)
        {
            fmpz_poly_t p;
            fmpq * ex;
            fmpz * c_array;
            slong * k_array;
            slong i, len, n_exact, n_interval;

            fmpz_poly_init(p);
            fmpz_poly_eulerian_polynomial(p, 2 + n_randint(state, 80));
            len = p->length;

            ex = _fmpq_vec_init(len);
            c_array = _fmpz_vec_init(len);
            k_array = flint_malloc(len * sizeof(slong));

            if (!_fmpz_poly_isolate_real_roots_signs(ex, &n_exact, c_array,
                    k_array, &n_interval, p->coeffs, len, 0,
                    _fmpz_poly_isolate_real_roots_signs_budget(len)) ||
                n_exact + n_interval != len - 1)
            {
                flint_printf("FAIL: sign change method with the default budget (Eulerian)\n");
                flint_printf("p = %{fmpz_poly}\n", p);
                flint_abort();
            }

            for (i = 0; i < n_interval; i++)
            {
                if (!check_one_root(p, c_array + i, k_array[i]))
                {
                    flint_printf("FAIL: interval (Eulerian)\n");
                    flint_printf("p = %{fmpz_poly}\n", p);
                    flint_abort();
                }
            }

            _fmpq_vec_clear(ex, len);
            _fmpz_vec_clear(c_array, len);
            flint_free(k_array);
            fmpz_poly_clear(p);
        }
    }

    /* regression test: exact positive roots when there is a root at zero */
    {
        fmpz_poly_t p;
        fmpq * ex = _fmpq_vec_init(3);
        fmpz * c = _fmpz_vec_init(3);
        slong k[3], ne, ni, i;

        fmpz_poly_init(p);
        fmpz_poly_set_str(p, "4  0 2 -3 1");   /* x (x - 1) (x - 2) */
        fmpz_poly_isolate_positive_roots(ex, &ne, c, k, &ni, p);
        if (ne + ni != 2)
        {
            flint_printf("FAIL: positive roots with a root at zero\n");
            flint_abort();
        }
        for (i = 0; i < ne; i++)
        {
            fmpq_t y;
            fmpq_init(y);
            fmpz_poly_evaluate_fmpq(y, p, ex + i);
            if (!fmpq_is_zero(y))
            {
                flint_printf("FAIL: exact positive root with a root at zero\n");
                flint_abort();
            }
            fmpq_clear(y);
        }
        fmpz_poly_clear(p);
        _fmpq_vec_clear(ex, 3);
        _fmpz_vec_clear(c, 3);
    }

    /* the sign change method must fail when Descartes' rule is not exact:
       x^4 - 3x^2 + 9 = (x^2 - 3x + 3)(x^2 + 3x + 3) has no real roots but
       two sign variations for x and -x */
    {
        fmpz_poly_t p;
        slong ne, ni;
        fmpz_poly_init(p);
        fmpz_poly_set_str(p, "5  9 0 -3 0 1");
        if (_fmpz_poly_isolate_real_roots_signs(NULL, &ne, NULL, NULL, &ni,
                p->coeffs, p->length, 0, 1000))
        {
            flint_printf("FAIL: sign change method on polynomial without real roots\n");
            flint_abort();
        }
        fmpz_poly_clear(p);
    }

    /* polynomials with irrational roots where the Sturm method is used */
    {
        int iter;

        for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
        {
            fmpz_poly_t p, q;
            fmpq * ex;
            fmpz * c_array;
            slong * k_array;
            slong i, n, n_exact, n_interval, count;
            int positive_only = n_randint(state, 4) == 0;

            fmpz_poly_init(p);
            fmpz_poly_init(q);

            switch (n_randint(state, 5))
            {
                case 3:
                    /* real roots spread over many orders of magnitude
                       (plain VCA, used for checking, is slow here) */
                    fmpz_poly_eulerian_polynomial(p, 1 + n_randint(state, 40));
                    break;
                case 4:
                    /* prod (x + 2^e_i + 1) for distinct e_i in [-25, 25],
                       with some roots made positive */
                    fmpz_poly_one(p);
                    for (n = -25; n <= 25; n++)
                    {
                        if (n_randint(state, 2) != 0)
                            continue;
                        fmpz_poly_zero(q);
                        fmpz_poly_set_coeff_si(q, 0, 1);
                        fmpz_poly_set_coeff_si(q, 1, 1);
                        if (n >= 0)
                            fmpz_mul_2exp(q->coeffs, q->coeffs, n);
                        else
                            fmpz_mul_2exp(q->coeffs + 1, q->coeffs + 1, -n);
                        fmpz_add(q->coeffs, q->coeffs, q->coeffs + 1);
                        if (n_randint(state, 4) == 0)
                            fmpz_neg(q->coeffs, q->coeffs);
                        fmpz_poly_mul(p, p, q);
                    }
                    fmpz_poly_zero(q);
                    break;
                case 0:
                    fmpz_poly_chebyshev_t(p, 1 + n_randint(state, 100));
                    break;
                case 1:
                    /* x^n - 2 (a x - 1)^2 */
                    /* (plain VCA, used for checking, is slow here) */
                    n = 5 + n_randint(state, 60);
                    fmpz_poly_set_coeff_si(q, 1, 3 + n_randint(state, 50));
                    fmpz_poly_set_coeff_si(q, 0, -1);
                    fmpz_poly_mul(p, q, q);
                    fmpz_poly_scalar_mul_si(p, p, -2);
                    fmpz_poly_set_coeff_si(q, 0, 0);
                    fmpz_poly_set_coeff_si(q, 1, 0);
                    fmpz_poly_set_coeff_si(q, n, 1);
                    fmpz_poly_add(p, p, q);
                    break;
                default:
                {
                    fmpq_poly_t r;
                    fmpq_poly_init(r);
                    fmpq_poly_laguerre_l(r, 1 + n_randint(state, 80));
                    fmpq_poly_get_numerator(p, r);
                    fmpq_poly_clear(r);
                }
            }

            /* random scaling and shift of the roots */
            _fmpz_poly_scale_2exp(p->coeffs, p->length, (slong) n_randint(state, 5) - 2);
            if (n_randint(state, 2))
            {
                fmpz_t c;
                fmpz_init(c);
                fmpz_set_si(c, (slong) n_randint(state, 5) - 2);
                fmpz_poly_taylor_shift(p, p, c);
                fmpz_clear(c);
            }

            /* times a factor without real roots */
            if (n_randint(state, 2))
            {
                fmpz_poly_randtest_no_real_root(q, state, 1 + n_randint(state, 10), 10);
                if (!fmpz_poly_is_zero(q))
                    fmpz_poly_mul(p, p, q);
            }

            if (fmpz_poly_is_zero(p) || !fmpz_poly_is_squarefree(p))
            {
                fmpz_poly_clear(p);
                fmpz_poly_clear(q);
                continue;
            }

            n = p->length;
            ex = _fmpq_vec_init(n);
            c_array = _fmpz_vec_init(n);
            k_array = flint_malloc(n * sizeof(slong));

            if (positive_only)
                fmpz_poly_isolate_positive_roots(ex, &n_exact, c_array, k_array, &n_interval, p);
            else
                fmpz_poly_isolate_real_roots(ex, &n_exact, c_array, k_array, &n_interval, p);

            /* expected count, computed with VCA */
            if (positive_only)
            {
                slong ne2, ni2;
                _fmpz_poly_isolate_real_roots_vca(NULL, &ne2, NULL, NULL, &ni2, p, 1);
                count = ne2 + ni2;
            }
            else
                count = fmpz_poly_num_real_roots_vca(p);

            if (n_exact + n_interval != count)
            {
                flint_printf("FAIL: wrong number of roots\n");
                flint_printf("p = %{fmpz_poly}\n", p);
                flint_printf("%wd + %wd, expected %wd\n", n_exact, n_interval, count);
                flint_abort();
            }

            for (i = 0; i < n_exact; i++)
            {
                fmpq_t y;
                fmpq_init(y);
                fmpz_poly_evaluate_fmpq(y, p, ex + i);
                if (!fmpq_is_zero(y) || (i > 0 && fmpq_cmp(ex + i - 1, ex + i) >= 0) ||
                    (positive_only && fmpq_sgn(ex + i) <= 0))
                {
                    flint_printf("FAIL: exact root\n");
                    flint_printf("p = %{fmpz_poly}\n", p);
                    flint_abort();
                }
                fmpq_clear(y);
            }

            for (i = 0; i < n_interval; i++)
            {
                if (!check_one_root(p, c_array + i, k_array[i]) ||
                    (positive_only && fmpz_sgn(c_array + i) < 0))
                {
                    flint_printf("FAIL: interval\n");
                    flint_printf("p = %{fmpz_poly}\n", p);
                    flint_printf("c = %{fmpz}, k = %wd\n", c_array + i, k_array[i]);
                    flint_abort();
                }

                /* intervals are sorted and disjoint */
                if (i > 0)
                {
                    fmpq_t a, b;
                    fmpq_init(a);
                    fmpq_init(b);
                    fmpz_add_ui(fmpq_numref(a), c_array + i - 1, 1);
                    fmpq_mul_2exp(a, a, 0);
                    if (k_array[i - 1] >= 0) fmpq_mul_2exp(a, a, k_array[i - 1]); else fmpq_div_2exp(a, a, -k_array[i - 1]);
                    fmpz_set(fmpq_numref(b), c_array + i);
                    if (k_array[i] >= 0) fmpq_mul_2exp(b, b, k_array[i]); else fmpq_div_2exp(b, b, -k_array[i]);
                    if (fmpq_cmp(a, b) > 0)
                    {
                        flint_printf("FAIL: intervals not sorted\n");
                        flint_printf("p = %{fmpz_poly}\n", p);
                        flint_abort();
                    }
                    fmpq_clear(a);
                    fmpq_clear(b);
                }
            }

            _fmpq_vec_clear(ex, n);
            _fmpz_vec_clear(c_array, n);
            flint_free(k_array);
            fmpz_poly_clear(p);
            fmpz_poly_clear(q);
        }
    }

    TEST_FUNCTION_END(state);
}

