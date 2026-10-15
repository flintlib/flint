/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpq.h"
#include "gr.h"
#include "gr_poly.h"

/* x^2 - 2 over QQ, then y^2 - 3 over QQ(x), then z^2 - x*y over that */
static void
_sqrt_tower(gr_ctx_t K1, gr_ctx_t K2, gr_ctx_t K3, gr_ctx_t QQ)
{
    gr_poly_t m;

    gr_poly_init(m, QQ);
    GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 2, 1, QQ));
    GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 0, -2, QQ));
    gr_ctx_init_gr_poly_quotient(K1, QQ, m);
    GR_MUST_SUCCEED(gr_ctx_set_gen_name(K1, "a"));
    GR_MUST_SUCCEED(gr_ctx_set_is_field(K1, T_TRUE));
    gr_poly_clear(m, QQ);

    gr_poly_init(m, K1);
    GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 2, 1, K1));
    GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 0, -3, K1));
    gr_ctx_init_gr_poly_quotient(K2, K1, m);
    GR_MUST_SUCCEED(gr_ctx_set_gen_name(K2, "b"));
    GR_MUST_SUCCEED(gr_ctx_set_is_field(K2, T_TRUE));
    gr_poly_clear(m, K1);

    {
        gr_ptr a, b, ab;
        GR_TMP_INIT3(a, b, ab, K2);
        GR_MUST_SUCCEED(gr_gen(b, K2));
        {
            gr_ptr a1;
            GR_TMP_INIT(a1, K1);
            GR_MUST_SUCCEED(gr_gen(a1, K1));
            GR_MUST_SUCCEED(gr_set_other(a, a1, K1, K2));
            GR_TMP_CLEAR(a1, K1);
        }
        GR_MUST_SUCCEED(gr_mul(ab, a, b, K2));
        gr_poly_init(m, K2);
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 2, 1, K2));
        GR_MUST_SUCCEED(gr_neg(ab, ab, K2));
        GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(m, 0, ab, K2));
        gr_ctx_init_gr_poly_quotient(K3, K2, m);
        GR_MUST_SUCCEED(gr_ctx_set_gen_name(K3, "c"));
        GR_MUST_SUCCEED(gr_ctx_set_is_field(K3, T_TRUE));
        gr_poly_clear(m, K2);
        GR_TMP_CLEAR3(a, b, ab, K2);
    }
}

TEST_FUNCTION_START(gr_poly_quotient, state)
{
    gr_ctx_t QQ, ZZ, K1, K2, K3, R;
    gr_poly_t m;
    int reps = 10 * flint_test_multiplier();
    slong iter;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_fmpz(ZZ);

    /* Random monic moduli over QQ and ZZ */
    for (iter = 0; iter < 20; iter++)
    {
        gr_ctx_struct * base = (iter % 2) ? QQ : ZZ;
        slong deg = 1 + n_randint(state, 5);

        gr_poly_init(m, base);
        GR_MUST_SUCCEED(gr_poly_randtest(m, state, deg, base));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, deg, 1, base));
        gr_ctx_init_gr_poly_quotient(R, base, m);
        gr_test_ring(R, reps, 0);
        gr_ctx_clear(R);
        gr_poly_clear(m, base);
    }

    /* Random monic moduli over GF(p): fields when the modulus is
       irreducible, other rings otherwise; also pretending to be fields */
    for (iter = 0; iter < 30; iter++)
    {
        gr_ctx_t Fp;
        slong deg = 1 + n_randint(state, 4);
        truth_t irred;
        int pretend = n_randint(state, 2);

        GR_MUST_SUCCEED(gr_ctx_init_nmod(Fp, n_randprime(state, 2 + n_randint(state, 6), 1)));
        GR_MUST_SUCCEED(gr_ctx_set_is_field(Fp, T_TRUE));

        gr_poly_init(m, Fp);
        GR_MUST_SUCCEED(gr_poly_randtest(m, state, deg, Fp));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, deg, 1, Fp));
        irred = gr_poly_is_irreducible(m, Fp);

        gr_ctx_init_gr_poly_quotient(R, Fp, m);
        if (irred != T_UNKNOWN)
            GR_MUST_SUCCEED(gr_ctx_set_is_field(R, irred));
        if (pretend)
            GR_MUST_SUCCEED(gr_ctx_set_is_pretend_field(R, T_TRUE));

        if (gr_ctx_is_field(R) != irred ||
            gr_ctx_is_pretend_field(R) != ((pretend || irred == T_TRUE) ? T_TRUE : T_FALSE))
        {
            flint_printf("FAIL: field predicates\n");
            gr_ctx_println(R);
            flint_abort();
        }

        gr_test_ring(R, reps, (irred == T_TRUE) ? GR_TEST_ALWAYS_ABLE : 0);

        gr_ctx_clear(R);
        gr_poly_clear(m, Fp);
        gr_ctx_clear(Fp);
    }

    /* Long and sparse moduli over GF(p), exercising the different
       representations of the precomputed modulus: products and powers
       against plain polynomial arithmetic */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t Fp;
        gr_poly_t a, b, c, r;
        gr_ptr x, y, z;
        slong i, deg = 1 + n_randint(state, n_randint(state, 4) ? 50 : 400);
        ulong e = n_randint(state, 100);
        int sparse = n_randint(state, 2);

        GR_MUST_SUCCEED(gr_ctx_init_nmod(Fp, n_randprime(state, 2 + n_randint(state, FLINT_BITS - 2), 1)));
        GR_MUST_SUCCEED(gr_ctx_set_is_field(Fp, T_TRUE));

        gr_poly_init(m, Fp);
        if (sparse)
        {
            for (i = 0; i < 4; i++)
                GR_MUST_SUCCEED(gr_poly_set_coeff_ui(m, n_randint(state, deg), n_randlimb(state), Fp));
        }
        else
        {
            GR_MUST_SUCCEED(gr_poly_randtest(m, state, deg, Fp));
        }
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, deg, 1, Fp));

        gr_ctx_init_gr_poly_quotient(R, Fp, m);
        gr_poly_init(a, Fp);
        gr_poly_init(b, Fp);
        gr_poly_init(c, Fp);
        gr_poly_init(r, Fp);
        x = gr_heap_init(R);
        y = gr_heap_init(R);
        z = gr_heap_init(R);

        /* (unreduced inputs too) */
        GR_MUST_SUCCEED(gr_poly_randtest(a, state, 1 + n_randint(state, 2 * deg), Fp));
        GR_MUST_SUCCEED(gr_poly_randtest(b, state, 1 + n_randint(state, 2 * deg), Fp));
        GR_MUST_SUCCEED(gr_poly_set(x, a, Fp));
        GR_MUST_SUCCEED(gr_poly_set(y, b, Fp));

        GR_MUST_SUCCEED(gr_mul(z, x, y, R));
        GR_MUST_SUCCEED(gr_poly_mul(c, a, b, Fp));
        GR_MUST_SUCCEED(gr_poly_rem(c, c, m, Fp));
        GR_MUST_SUCCEED(gr_poly_quotient_get_poly(r, z, R));
        if (gr_poly_equal(r, c, Fp) != T_TRUE)
        {
            flint_printf("FAIL: mul (deg %wd)\n", deg);
            flint_abort();
        }

        /* aliased squaring and powers */
        GR_MUST_SUCCEED(gr_set(z, x, R));
        GR_MUST_SUCCEED(gr_sqr(z, z, R));
        GR_MUST_SUCCEED(gr_poly_mul(c, a, a, Fp));
        GR_MUST_SUCCEED(gr_poly_rem(c, c, m, Fp));
        GR_MUST_SUCCEED(gr_poly_quotient_get_poly(r, z, R));
        if (gr_poly_equal(r, c, Fp) != T_TRUE)
        {
            flint_printf("FAIL: sqr (deg %wd)\n", deg);
            flint_abort();
        }

        GR_MUST_SUCCEED(gr_set(z, x, R));
        GR_MUST_SUCCEED(gr_pow_ui(z, z, e, R));
        GR_MUST_SUCCEED(gr_poly_one(c, Fp));
        for (i = 0; i < (slong) e; i++)
        {
            GR_MUST_SUCCEED(gr_poly_mul(c, c, a, Fp));
            GR_MUST_SUCCEED(gr_poly_rem(c, c, m, Fp));
        }
        GR_MUST_SUCCEED(gr_poly_rem(c, c, m, Fp));
        GR_MUST_SUCCEED(gr_poly_quotient_get_poly(r, z, R));
        if (gr_poly_equal(r, c, Fp) != T_TRUE)
        {
            flint_printf("FAIL: pow_ui (deg %wd, e %wu)\n", deg, e);
            flint_abort();
        }

        /* negative powers: x^-e x^e = 1 when x is invertible */
        if (gr_pow_si(z, x, -(slong) e, R) == GR_SUCCESS)
        {
            GR_MUST_SUCCEED(gr_pow_ui(y, x, e, R));
            GR_MUST_SUCCEED(gr_mul(z, z, y, R));
            if (gr_is_one(z, R) != T_TRUE)
            {
                flint_printf("FAIL: pow_si (deg %wd, e %wu)\n", deg, e);
                flint_abort();
            }
        }

        gr_heap_clear(x, R);
        gr_heap_clear(y, R);
        gr_heap_clear(z, R);
        gr_poly_clear(a, Fp);
        gr_poly_clear(b, Fp);
        gr_poly_clear(c, Fp);
        gr_poly_clear(r, Fp);
        gr_ctx_clear(R);
        gr_poly_clear(m, Fp);
        gr_ctx_clear(Fp);
    }

    /* Tower of square roots */
    _sqrt_tower(K1, K2, K3, QQ);
    /* gr_test_ring renames the generator, so restore distinct names
       between runs to keep string parsing unambiguous */
    gr_test_ring(K1, reps, GR_TEST_ALWAYS_ABLE);
    GR_MUST_SUCCEED(gr_ctx_set_gen_name(K1, "a"));
    gr_test_ring(K2, reps, GR_TEST_ALWAYS_ABLE);
    GR_MUST_SUCCEED(gr_ctx_set_gen_name(K2, "b"));
    gr_test_ring(K3, reps, GR_TEST_ALWAYS_ABLE);
    GR_MUST_SUCCEED(gr_ctx_set_gen_name(K3, "c"));

    /* c^2 = a*b, so (c - ...)? check a few identities: (a*b)^2 = 6, c^4 = 6 */
    {
        gr_ptr c, t;
        GR_TMP_INIT2(c, t, K3);
        GR_MUST_SUCCEED(gr_gen(c, K3));
        GR_MUST_SUCCEED(gr_pow_ui(t, c, 4, K3));
        if (gr_sub_si(t, t, 6, K3) != GR_SUCCESS || gr_is_zero(t, K3) != T_TRUE)
        {
            flint_printf("FAIL: c^4 != 6\n");
            gr_println(t, K3);
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_inv(t, c, K3));
        GR_MUST_SUCCEED(gr_mul(t, t, c, K3));
        if (gr_is_one(t, K3) != T_TRUE)
        {
            flint_printf("FAIL: c * c^-1 != 1\n");
            flint_abort();
        }
        GR_TMP_CLEAR2(c, t, K3);
    }

    gr_ctx_clear(K3);
    gr_ctx_clear(K2);
    gr_ctx_clear(K1);

    /* Dynamic evaluation: QQ[x]/(x^2 - 1), inverting x - 1 finds a zero divisor;
       refining to x - 1 makes x = 1. */
    {
        gr_ptr x, y, t;
        const gr_poly_struct * g;
        int status;

        gr_poly_init(m, QQ);
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 2, 1, QQ));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 0, -1, QQ));
        gr_ctx_init_gr_poly_quotient(R, QQ, m);
        gr_poly_clear(m, QQ);

        GR_TMP_INIT3(x, y, t, R);
        GR_MUST_SUCCEED(gr_gen(x, R));
        GR_MUST_SUCCEED(gr_sub_si(y, x, 1, R));    /* y = x - 1 */

        /* as a ring: x - 1 is simply not invertible */
        if (gr_ctx_is_field(R) != T_UNKNOWN || gr_ctx_is_pretend_field(R) != T_FALSE ||
            gr_inv(t, y, R) != GR_DOMAIN || gr_is_invertible(y, R) != T_FALSE ||
            gr_poly_quotient_ctx_num_zero_divisors(R) != 0)
        {
            flint_printf("FAIL: zero divisor in a ring which does not pretend to be a field\n");
            flint_abort();
        }

        /* D5 assumption: the zero divisor contradicts the pretense */
        GR_MUST_SUCCEED(gr_ctx_set_is_pretend_field(R, T_TRUE));
        if (gr_ctx_is_field(R) != T_UNKNOWN || gr_ctx_is_pretend_field(R) != T_TRUE)
        {
            flint_printf("FAIL: is_pretend_field\n");
            flint_abort();
        }

        status = gr_inv(t, y, R);
        if (status != GR_UNABLE || gr_poly_quotient_ctx_num_zero_divisors(R) != 1 ||
            gr_is_invertible(y, R) != T_UNKNOWN)
        {
            flint_printf("FAIL: expected GR_UNABLE and one zero divisor, got status %d, %wd\n",
                status, gr_poly_quotient_ctx_num_zero_divisors(R));
            flint_abort();
        }

        /* the zero divisor is also available as an element */
        GR_MUST_SUCCEED(gr_ctx_recover_zero_divisor(t, R));
        if (gr_is_zero(t, R) != T_FALSE || gr_mul(t, t, y, R) != GR_SUCCESS)
        {
            flint_printf("FAIL: recover_zero_divisor\n");
            flint_abort();
        }

        g = gr_poly_quotient_ctx_zero_divisor(R, 0);
        /* g must be x - 1 or x + 1 */
        if (g->length != 2)
        {
            flint_printf("FAIL: zero divisor has wrong degree\n");
            flint_abort();
        }

        /* refine to x - 1 (choose the factor vanishing at "our" root x = 1) */
        gr_poly_init(m, QQ);
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 1, 1, QQ));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 0, -1, QQ));
        GR_MUST_SUCCEED(gr_poly_quotient_ctx_refine(R, m));
        gr_poly_clear(m, QQ);

        /* stale element x should now equal 1 */
        if (gr_is_one(x, R) != T_TRUE || gr_is_zero(y, R) != T_TRUE)
        {
            flint_printf("FAIL: stale elements not reduced correctly after refinement\n");
            gr_println(x, R);
            gr_println(y, R);
            flint_abort();
        }

        GR_MUST_SUCCEED(gr_add(t, x, x, R));
        if (gr_sub_si(t, t, 2, R) != GR_SUCCESS || gr_is_zero(t, R) != T_TRUE)
        {
            flint_printf("FAIL: x + x != 2 after refinement\n");
            flint_abort();
        }

        gr_test_ring(R, reps, GR_TEST_ALWAYS_ABLE);

        GR_TMP_CLEAR3(x, y, t, R);
        gr_ctx_clear(R);
    }

    /* Norms: N(M1 M2) = N(M1) N(M2) over random quotient rings, and
       known values in Q(sqrt(2)) */
    for (iter = 0; iter < 30; iter++)
    {
        gr_ctx_struct * base = QQ;
        slong deg = 1 + n_randint(state, 4);
        gr_poly_t M1, M2, M3, N1, N2, N3;
        gr_ptr x, y, nx, ny, nxy;
        int status = GR_SUCCESS;

        gr_poly_init(m, base);
        GR_MUST_SUCCEED(gr_poly_randtest(m, state, deg, base));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, deg, 1, base));
        gr_ctx_init_gr_poly_quotient(R, base, m);

        gr_poly_init(M1, R); gr_poly_init(M2, R); gr_poly_init(M3, R);
        gr_poly_init(N1, base); gr_poly_init(N2, base); gr_poly_init(N3, base);
        GR_TMP_INIT2(x, y, R);
        GR_TMP_INIT3(nx, ny, nxy, base);

        status |= gr_poly_randtest(M1, state, 1 + n_randint(state, 4), R);
        status |= gr_poly_randtest(M2, state, 1 + n_randint(state, 4), R);
        status |= gr_poly_mul(M3, M1, M2, R);
        status |= gr_poly_quotient_norm_poly(N1, M1, R);
        status |= gr_poly_quotient_norm_poly(N2, M2, R);
        status |= gr_poly_quotient_norm_poly(N3, M3, R);
        status |= gr_poly_mul(N1, N1, N2, base);

        if (status == GR_SUCCESS && gr_poly_equal(N1, N3, base) == T_FALSE)
        {
            flint_printf("FAIL: norm_poly not multiplicative\n");
            flint_printf("m = "); gr_poly_print(m, base); flint_printf("\n");
            flint_printf("M1 = "); gr_poly_print(M1, R); flint_printf("\n");
            flint_printf("M2 = "); gr_poly_print(M2, R); flint_printf("\n");
            flint_abort();
        }

        status |= gr_randtest(x, state, R);
        status |= gr_randtest(y, state, R);
        status |= gr_poly_quotient_norm(nx, x, R);
        status |= gr_poly_quotient_norm(ny, y, R);
        status |= gr_mul(x, x, y, R);
        status |= gr_poly_quotient_norm(nxy, x, R);
        status |= gr_mul(nx, nx, ny, base);

        if (status == GR_SUCCESS && gr_equal(nx, nxy, base) == T_FALSE)
        {
            flint_printf("FAIL: norm not multiplicative\n");
            flint_printf("m = "); gr_poly_print(m, base); flint_printf("\n");
            flint_abort();
        }

        gr_poly_clear(M1, R); gr_poly_clear(M2, R); gr_poly_clear(M3, R);
        gr_poly_clear(N1, base); gr_poly_clear(N2, base); gr_poly_clear(N3, base);
        GR_TMP_CLEAR2(x, y, R);
        GR_TMP_CLEAR3(nx, ny, nxy, base);
        gr_ctx_clear(R);
        gr_poly_clear(m, base);
    }

    /* the norm of 1 + sqrt(2) is -1, that of x - sqrt(2) is x^2 - 2 */
    {
        gr_poly_t M, N;
        gr_ptr x, n;
        int status = GR_SUCCESS;

        gr_poly_init(m, QQ);
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 2, 1, QQ));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, 0, -2, QQ));
        gr_ctx_init_gr_poly_quotient(R, QQ, m);
        gr_poly_init(M, R);
        gr_poly_init(N, QQ);
        GR_TMP_INIT(x, R);
        GR_TMP_INIT(n, QQ);

        status |= gr_gen(x, R);
        status |= gr_add_si(x, x, 1, R);
        status |= gr_poly_quotient_norm(n, x, R);
        if (status != GR_SUCCESS || gr_is_neg_one(n, QQ) != T_TRUE)
        {
            flint_printf("FAIL: N(1 + sqrt(2)) != -1\n");
            gr_println(n, QQ);
            flint_abort();
        }

        status |= gr_gen(x, R);
        status |= gr_poly_set_coeff_si(M, 1, 1, R);
        status |= gr_neg(x, x, R);
        status |= gr_poly_set_coeff_scalar(M, 0, x, R);
        status |= gr_poly_quotient_norm_poly(N, M, R);
        status |= gr_poly_set(m, N, QQ);
        if (status != GR_SUCCESS || gr_poly_equal(N, m, QQ) != T_TRUE || N->length != 3)
        {
            flint_printf("FAIL: N(x - sqrt(2)) != x^2 - 2\n");
            gr_poly_print(N, QQ); flint_printf("\n");
            flint_abort();
        }

        gr_poly_clear(M, R);
        gr_poly_clear(N, QQ);
        GR_TMP_CLEAR(x, R);
        GR_TMP_CLEAR(n, QQ);
        gr_ctx_clear(R);
        gr_poly_clear(m, QQ);
    }

    gr_ctx_clear(QQ);
    gr_ctx_clear(ZZ);

    TEST_FUNCTION_END(state);
}
