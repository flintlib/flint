/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "acb.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

/*
    Factorization of polynomials over towers (Trager's method): random
    products over random radical towers are factored into monic
    irreducible factors whose product is the input (each factor factors
    no further), known factorizations over algebraic and transcendental
    towers, the lazy fields (linear factors, and quadratic factors for
    the pairs of nonreal roots over the real field), and the refinement
    of a reducible dynamic step over a transcendental tower.
*/

/* a random element: small integer combination of 1 and the generators */
static int
_small_elem(gr_ptr x, flint_rand_t state, gr_tower_t T)
{
    gr_ctx_struct * top = gr_tower_field(T);
    gr_ptr t;
    slong d;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, top);
    status |= gr_set_si(x, (slong) n_randint(state, 7) - 3, top);
    for (d = 0; d < T->num_gens; d++)
    {
        if (n_randint(state, 2))
        {
            status |= gr_tower_gen_get(t, T, d);
            status |= gr_mul_si(t, t, (slong) n_randint(state, 5) - 2, top);
            status |= gr_add(x, x, t, top);
        }
    }
    GR_TMP_CLEAR(t, top);
    return status;
}

static int
_small_poly(gr_poly_t f, slong deg, flint_rand_t state, gr_tower_t T)
{
    gr_ctx_struct * top = gr_tower_field(T);
    slong i;
    int status = GR_SUCCESS;

    gr_poly_fit_length(f, deg + 1, top);
    for (i = 0; i < deg; i++)
        status |= _small_elem(gr_poly_coeff_ptr(f, i, top), state, T);
    status |= gr_set_si(gr_poly_coeff_ptr(f, deg, top), 1 + n_randint(state, 2), top);
    _gr_poly_set_length(f, deg + 1, top);
    _gr_poly_normalise(f, top);
    return status;
}

/* checks c prod fac^e == f, the factors monic of degree >= 1; returns
   the number of factors with multiplicity */
static slong
_check_factorization(gr_srcptr c, gr_vec_t fac, const fmpz_vec_t e, const gr_poly_t f, gr_ctx_t F, const char * what)
{
    gr_ctx_t P;
    gr_poly_t prod;
    slong i, j, num = 0;

    gr_ctx_init_gr_poly(P, F);
    gr_poly_init(prod, F);
    GR_MUST_SUCCEED(gr_poly_set_scalar(prod, c, F));

    for (i = 0; i < fac->length; i++)
    {
        const gr_poly_struct * h = gr_vec_entry_ptr(fac, i, P);
        if (h->length < 2 || gr_is_one(gr_poly_coeff_srcptr(h, h->length - 1, F), F) != T_TRUE)
        {
            flint_printf("FAIL (%s): factor not monic of positive degree\n", what);
            flint_abort();
        }
        for (j = 0; j < fmpz_get_si(e->entries + i); j++)
            GR_MUST_SUCCEED(gr_poly_mul(prod, prod, h, F));
        num += fmpz_get_si(e->entries + i);
    }

    if (gr_poly_equal(prod, f, F) != T_TRUE)
    {
        flint_printf("FAIL (%s): product of the factors\n", what);
        flint_printf("f = "); gr_poly_print(f, F); flint_printf("\n");
        flint_printf("prod = "); gr_poly_print(prod, F); flint_printf("\n");
        flint_abort();
    }

    gr_poly_clear(prod, F);
    gr_ctx_clear(P);
    return num;
}

static void
_expect_factors(gr_tower_t T, const gr_poly_t f, slong num_expected, const char * what)
{
    gr_ctx_struct * F = gr_tower_field(T);
    gr_ctx_t P;
    gr_vec_t fac;
    fmpz_vec_t e;
    gr_ptr c;
    slong num;

    gr_ctx_init_gr_poly(P, F);
    gr_vec_init(fac, 0, P);
    fmpz_vec_init(e, 0);
    GR_TMP_INIT(c, F);

    if (gr_tower_poly_factor(c, fac, e, f, T->length, T) != GR_SUCCESS)
    {
        flint_printf("FAIL (%s): not factored\n", what);
        flint_abort();
    }
    num = _check_factorization(c, fac, e, f, F, what);
    if (num != num_expected)
    {
        flint_printf("FAIL (%s): %wd factors, expected %wd\n", what, num, num_expected);
        gr_tower_print(T);
        flint_abort();
    }

    GR_TMP_CLEAR(c, F);
    fmpz_vec_clear(e);
    gr_vec_clear(fac, P);
    gr_ctx_clear(P);
}

/* x^n + c0 (c0 an element of the top field) */
static void
_binomial(gr_poly_t f, slong n, gr_srcptr c0, gr_tower_t T)
{
    gr_ctx_struct * F = gr_tower_field(T);
    GR_MUST_SUCCEED(gr_poly_zero(f, F));
    GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, n, 1, F));
    GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(f, 0, c0, F));
}

TEST_FUNCTION_START(gr_tower_factor, state)
{
    gr_ctx_t QQ;
    slong iter;

    gr_ctx_init_fmpq(QQ);

    /* random products over random radical towers */
    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        gr_tower_t T;
        gr_ctx_struct * F;
        gr_ctx_t P;
        gr_poly_t f, g, h;
        gr_vec_t fac, fac2;
        fmpz_vec_t e, e2;
        gr_ptr c, c2;
        slong steps = n_randint(state, 3), k, i, num, e_h;
        int status = GR_SUCCESS;

        gr_tower_init(T, QQ);
        for (k = 0; k < steps && status == GR_SUCCESS; k++)
        {
            gr_ctx_struct * top = gr_tower_field(T);
            gr_ptr x;
            GR_TMP_INIT(x, top);
            status |= gr_set_ui(x, 2 + n_randint(state, 10), top);
            status |= gr_tower_adjoin_root_ui(T, x, 2 + n_randint(state, 2), NULL);
            GR_TMP_CLEAR(x, top);
        }

        if (status != GR_SUCCESS || gr_tower_degree(T) > 6)
        {
            gr_tower_clear(T);
            continue;
        }

        F = gr_tower_field(T);
        gr_ctx_init_gr_poly(P, F);
        gr_poly_init(f, F);
        gr_poly_init(g, F);
        gr_poly_init(h, F);
        gr_vec_init(fac, 0, P);
        gr_vec_init(fac2, 0, P);
        fmpz_vec_init(e, 0);
        fmpz_vec_init(e2, 0);
        GR_TMP_INIT2(c, c2, F);

        /* f = g h^e_h */
        status |= _small_poly(g, 1 + n_randint(state, 3), state, T);
        status |= _small_poly(h, 1 + n_randint(state, 2), state, T);
        e_h = 1 + n_randint(state, 2);
        status |= gr_poly_pow_ui(f, h, e_h, F);
        status |= gr_poly_mul(f, f, g, F);

        if (status == GR_SUCCESS)
        {
            status = gr_tower_poly_factor(c, fac, e, f, T->length, T);
            if (status != GR_SUCCESS)
            {
                flint_printf("FAIL: not factored (status %d)\n", status);
                gr_tower_print(T);
                flint_printf("f = "); gr_poly_print(f, F); flint_printf("\n");
                flint_abort();
            }

            num = _check_factorization(c, fac, e, f, F, "random");
            if (num < 1 + e_h)
            {
                flint_printf("FAIL: too few factors\n");
                flint_abort();
            }

            /* the factors are irreducible */
            for (i = 0; i < fac->length; i++)
            {
                GR_MUST_SUCCEED(gr_tower_poly_factor(c2, fac2, e2, gr_vec_entry_ptr(fac, i, P), T->length, T));
                if (fac2->length != 1 || !fmpz_is_one(e2->entries))
                {
                    flint_printf("FAIL: a factor factors further\n");
                    gr_tower_print(T);
                    flint_printf("factor = "); gr_poly_print(gr_vec_entry_ptr(fac, i, P), F); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        GR_TMP_CLEAR2(c, c2, F);
        fmpz_vec_clear(e);
        fmpz_vec_clear(e2);
        gr_vec_clear(fac, P);
        gr_vec_clear(fac2, P);
        gr_poly_clear(f, F);
        gr_poly_clear(g, F);
        gr_poly_clear(h, F);
        gr_ctx_clear(P);
        gr_tower_clear(T);
    }

    /* known factorizations over Q(sqrt 2, sqrt 3) */
    {
        gr_tower_t T;
        gr_ctx_struct * F;
        gr_poly_t f;
        gr_ptr x;
        slong k;

        gr_tower_init(T, QQ);
        for (k = 0; k < 2; k++)
        {
            gr_ctx_struct * top = gr_tower_field(T);
            GR_TMP_INIT(x, top);
            GR_MUST_SUCCEED(gr_set_ui(x, 2 + k, top));
            GR_MUST_SUCCEED(gr_tower_adjoin_sqrt(T, x, NULL));
            GR_TMP_CLEAR(x, top);
        }
        F = gr_tower_field(T);
        gr_poly_init(f, F);
        GR_TMP_INIT(x, F);

        /* x^4 + 1 = (x^2 + sqrt2 x + 1)(x^2 - sqrt2 x + 1) */
        GR_MUST_SUCCEED(gr_one(x, F));
        _binomial(f, 4, x, T);
        _expect_factors(T, f, 2, "x^4 + 1");

        /* x^4 - 10 x^2 + 1: the four conjugates of sqrt 2 + sqrt 3 */
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 2, -10, F));
        _expect_factors(T, f, 4, "x^4 - 10 x^2 + 1");

        /* x^3 - 2: irreducible */
        GR_MUST_SUCCEED(gr_set_si(x, -2, F));
        _binomial(f, 3, x, T);
        _expect_factors(T, f, 1, "x^3 - 2");

        /* (x^2 - 6)^3 */
        GR_MUST_SUCCEED(gr_set_si(x, -6, F));
        _binomial(f, 2, x, T);
        GR_MUST_SUCCEED(gr_poly_pow_ui(f, f, 3, F));
        _expect_factors(T, f, 6, "(x^2 - 6)^3");

        GR_TMP_CLEAR(x, F);
        gr_poly_clear(f, F);
        gr_tower_clear(T);
    }

    /* over Q(pi)(sqrt 2) */
    {
        gr_tower_t T;
        gr_ctx_struct * F;
        gr_poly_t f;
        gr_ptr x, pi;

        gr_tower_init(T, QQ);
        GR_MUST_SUCCEED(gr_tower_adjoin_pi(T, NULL));
        {
            gr_ctx_struct * top = gr_tower_field(T);
            GR_TMP_INIT(x, top);
            GR_MUST_SUCCEED(gr_set_ui(x, 2, top));
            GR_MUST_SUCCEED(gr_tower_adjoin_sqrt(T, x, NULL));
            GR_TMP_CLEAR(x, top);
        }
        F = gr_tower_field(T);
        gr_poly_init(f, F);
        GR_TMP_INIT2(x, pi, F);
        GR_MUST_SUCCEED(gr_tower_gen_get(pi, T, 0));

        /* x^4 + pi^4 = (x^2 + sqrt2 pi x + pi^2)(x^2 - sqrt2 pi x + pi^2) */
        GR_MUST_SUCCEED(gr_pow_ui(x, pi, 4, F));
        _binomial(f, 4, x, T);
        _expect_factors(T, f, 2, "x^4 + pi^4");

        /* x^2 - 2 pi^2 = (x - sqrt2 pi)(x + sqrt2 pi) */
        GR_MUST_SUCCEED(gr_sqr(x, pi, F));
        GR_MUST_SUCCEED(gr_mul_si(x, x, -2, F));
        _binomial(f, 2, x, T);
        _expect_factors(T, f, 2, "x^2 - 2 pi^2");

        /* x^3 - pi: irreducible */
        GR_MUST_SUCCEED(gr_neg(x, pi, F));
        _binomial(f, 3, x, T);
        _expect_factors(T, f, 1, "x^3 - pi");

        /* x^4 + pi^4 adjoined as a dynamic step: refined by Trager's
           method to a proven quadratic step */
        {
            acb_t z;
            gr_poly_t m;
            acb_init(z);
            gr_poly_init(m, F);
            GR_MUST_SUCCEED(gr_pow_ui(x, pi, 4, F));
            _binomial(m, 4, x, T);
            /* the root pi exp(i pi / 4) */
            acb_set_d_d(z, 2.221441469079183, 2.221441469079183);
            GR_MUST_SUCCEED(gr_tower_adjoin_algebraic(T, m, z, GR_TOWER_STATUS_DYNAMIC, NULL));
            if (!gr_tower_prove_step_trager(T, T->length, GR_TOWER_OPTION(T, GR_TOWER_OPT_TRAGER_DEGREE_LIMIT)) ||
                gr_tower_step_degree(T, T->length) != 2)
            {
                flint_printf("FAIL: x^4 + pi^4 over Q(pi, sqrt 2) not refined\n");
                gr_tower_print(T);
                flint_abort();
            }
            gr_poly_clear(m, F);
            acb_clear(z);
        }

        GR_TMP_CLEAR2(x, pi, F);
        gr_poly_clear(f, F);
        gr_tower_clear(T);
    }

    /* the lazy fields */
    {
        gr_ctx_t C, R, PC, PR;
        gr_poly_t f;
        gr_vec_t fac;
        fmpz_vec_t e;
        gr_poly_t c;
        slong num, i, lin;

        gr_ctx_init_tower_lazy(C, QQ, GR_TOWER_MERGE_EXPRESS);
        gr_ctx_init_tower_lazy(R, QQ, GR_TOWER_MERGE_EXPRESS | GR_TOWER_LAZY_REAL);
        gr_ctx_init_gr_poly(PC, C);
        gr_ctx_init_gr_poly(PR, R);

        /* x^4 + pi x + 1 over the complex field: four linear factors */
        gr_poly_init(f, C);
        gr_poly_init(c, C);
        gr_vec_init(fac, 0, PC);
        fmpz_vec_init(e, 0);
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 4, 1, C));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 0, 1, C));
        {
            gr_ptr pi;
            GR_TMP_INIT(pi, C);
            GR_MUST_SUCCEED(gr_pi(pi, C));
            GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(f, 1, pi, C));
            GR_TMP_CLEAR(pi, C);
        }
        GR_MUST_SUCCEED(gr_factor(c, fac, e, f, 0, PC));
        num = _check_factorization(gr_poly_coeff_srcptr(c, 0, C), fac, e, f, C, "lazy complex");
        if (num != 4 || fac->length != 4)
        {
            flint_printf("FAIL: x^4 + pi x + 1 over the complex field\n");
            flint_abort();
        }
        gr_vec_clear(fac, PC);
        gr_poly_clear(f, C);
        gr_poly_clear(c, C);

        /* over the real field: x^4 + 1 (two quadratic factors), and
           (x^3 - 2)(x - 1)^2 (a linear, a quadratic and a double factor) */
        gr_poly_init(f, R);
        gr_poly_init(c, R);
        gr_vec_init(fac, 0, PR);
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 4, 1, R));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 0, 1, R));
        GR_MUST_SUCCEED(gr_factor(c, fac, e, f, 0, PR));
        num = _check_factorization(gr_poly_coeff_srcptr(c, 0, R), fac, e, f, R, "lazy real x^4 + 1");
        if (fac->length != 2 || ((gr_poly_struct *) gr_vec_entry_ptr(fac, 0, PR))->length != 3)
        {
            flint_printf("FAIL: x^4 + 1 over the real field\n");
            flint_abort();
        }

        {
            gr_poly_t g;
            gr_poly_init(g, R);
            GR_MUST_SUCCEED(gr_poly_zero(f, R));
            GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 3, 1, R));
            GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 0, -2, R));
            GR_MUST_SUCCEED(gr_poly_set_coeff_si(g, 1, 1, R));
            GR_MUST_SUCCEED(gr_poly_set_coeff_si(g, 0, -1, R));
            GR_MUST_SUCCEED(gr_poly_mul(g, g, g, R));
            GR_MUST_SUCCEED(gr_poly_mul(f, f, g, R));
            gr_poly_clear(g, R);
        }
        GR_MUST_SUCCEED(gr_factor(c, fac, e, f, 0, PR));
        num = _check_factorization(gr_poly_coeff_srcptr(c, 0, R), fac, e, f, R, "lazy real (x^3 - 2)(x - 1)^2");
        for (lin = 0, i = 0; i < fac->length; i++)
            lin += (((gr_poly_struct *) gr_vec_entry_ptr(fac, i, PR))->length == 2);
        if (fac->length != 3 || lin != 2 || num != 4)
        {
            flint_printf("FAIL: (x^3 - 2)(x - 1)^2 over the real field\n");
            flint_abort();
        }

        gr_vec_clear(fac, PR);
        fmpz_vec_clear(e);
        gr_poly_clear(f, R);
        gr_poly_clear(c, R);
        gr_ctx_clear(PC);
        gr_ctx_clear(PR);
        gr_ctx_clear(C);
        gr_ctx_clear(R);
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
