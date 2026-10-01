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
#include "fmpz.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_mat.h"
#include "mpn_mod.h"

/*
    Rings Z/nZ (nmod, fmpz_mod, mpn_mod) with n composite, and quotient
    rings (Z/nZ)[x]/(m) with m possibly reducible, pretending to be fields: every operation either succeeds with a correct result, or
    fails with GR_UNABLE having recorded a zero divisor, which
    gr_ctx_recover_zero_divisor returns (a nontrivial factor of n); a
    domain error only means what it means in a field.
*/

#define CHECK(cond, msg) do { if (!(cond)) { flint_printf("FAIL: %s\n", msg); gr_ctx_println(ctx); flint_abort(); } } while (0)

/* whether a zero divisor has been recorded: a nontrivial factor of n
   for Z/nZ; a nonzero non-unit for other rings (a quotient ring may
   also have met the zero divisor in its base ring) */
static int
_check_zero_divisor_generic(gr_ctx_t ctx)
{
    gr_ptr r, t;
    int ok;

    r = gr_heap_init(ctx);
    t = gr_heap_init(ctx);
    ok = (gr_ctx_recover_zero_divisor(r, ctx) == GR_SUCCESS) &&
         (gr_is_zero(r, ctx) == T_FALSE) &&
         (gr_inv(t, r, ctx) != GR_SUCCESS);
    gr_heap_clear(r, ctx);
    gr_heap_clear(t, ctx);
    return ok;
}

static int
_check_zero_divisor_modular(gr_ctx_t ctx, const fmpz_t n)
{
    gr_ptr r;
    fmpz_t c, g;
    int ok;

    r = gr_heap_init(ctx);
    fmpz_init(c);
    fmpz_init(g);
    ok = (gr_ctx_recover_zero_divisor(r, ctx) == GR_SUCCESS) &&
         (gr_get_fmpz(c, r, ctx) == GR_SUCCESS);
    if (ok)
    {
        fmpz_gcd(g, c, n);
        ok = !fmpz_is_one(g) && !fmpz_equal(g, n);
    }
    fmpz_clear(c);
    fmpz_clear(g);
    gr_heap_clear(r, ctx);
    return ok;
}

static int
_check_zero_divisor(gr_ctx_t ctx, gr_ctx_t base, const fmpz_t n)
{
    if (base == NULL)
        return _check_zero_divisor_modular(ctx, n);
    return _check_zero_divisor_generic(ctx) || _check_zero_divisor_modular(base, n);
}

static void
_rand_poly(gr_poly_t f, flint_rand_t state, slong len, gr_ctx_t ctx)
{
    GR_MUST_SUCCEED(gr_poly_randtest(f, state, len, ctx));
}

TEST_FUNCTION_START(gr_pretend_field, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, base_ctx;
        gr_ctx_struct * base = NULL;
        fmpz_t n, p, q;
        int which = n_randint(state, 4), status;
        ulong small[] = {2, 3, 5, 7, 11};

        fmpz_init(n);
        fmpz_init(p);
        fmpz_init(q);

        /* n = p q with a small prime p, so that zero divisors occur
           often; sometimes n prime */
        fmpz_set_ui(p, small[n_randint(state, 5)]);
        if (which == 0 || which == 3)
            fmpz_set_ui(q, n_randprime(state, 2 + n_randint(state, which == 3 ? 10 : 40), 1));
        else
        {
            fmpz_randprime(q, state, (which == 2 ? 2 * FLINT_BITS : 2) + n_randint(state, 100), 1);
        }
        if (n_randint(state, 8) == 0)
            fmpz_one(p);
        fmpz_mul(n, p, q);

        if (which == 3)
        {
            /* (Z/nZ)[x]/(m) for random monic m (often reducible), over
               a base ring which pretends to be a field too */
            gr_poly_t m;
            slong d = 1 + n_randint(state, 3);
            GR_MUST_SUCCEED(gr_ctx_init_nmod(base_ctx, fmpz_get_ui(n)));
            GR_MUST_SUCCEED(gr_ctx_set_is_pretend_field(base_ctx, T_TRUE));
            base = base_ctx;
            gr_poly_init(m, base);
            GR_MUST_SUCCEED(gr_poly_randtest(m, state, d + 1, base));
            GR_MUST_SUCCEED(gr_poly_set_coeff_si(m, d, 1, base));
            gr_ctx_init_gr_poly_quotient(ctx, base, m);
            gr_poly_clear(m, base);
        }
        else if (which == 0)
            GR_MUST_SUCCEED(gr_ctx_init_nmod(ctx, fmpz_get_ui(n)));
        else if (which == 1)
            gr_ctx_init_fmpz_mod(ctx, n);
        else
        {
            if (fmpz_size(n) < 2 || gr_ctx_init_mpn_mod(ctx, n) != GR_SUCCESS)
            {
                fmpz_clear(n);
                fmpz_clear(p);
                fmpz_clear(q);
                continue;
            }
        }

        /* not pretending: domain errors for non-units, nothing recorded */
        CHECK(gr_ctx_is_pretend_field(ctx) != T_TRUE || (fmpz_is_one(p) && base == NULL), "pretend by default");
        {
            gr_ptr x, y;
            x = gr_heap_init(ctx);
            y = gr_heap_init(ctx);
            if (!fmpz_is_one(p) && base == NULL)
            {
                GR_MUST_SUCCEED(gr_set_fmpz(x, p, ctx));
                CHECK(gr_inv(y, x, ctx) == GR_DOMAIN, "inverse of a non-unit, not pretending");
                CHECK(gr_ctx_recover_zero_divisor(y, ctx) == GR_UNABLE, "recorded without pretending");
            }
            gr_heap_clear(x, ctx);
            gr_heap_clear(y, ctx);
        }

        GR_MUST_SUCCEED(gr_ctx_set_is_pretend_field(ctx, T_TRUE));
        CHECK(gr_ctx_is_pretend_field(ctx) == T_TRUE, "set_is_pretend_field");
        CHECK(gr_ctx_is_field(ctx) != T_TRUE || (fmpz_is_one(p) && base == NULL), "is_field for a composite modulus");

        /* inverses */
        {
            gr_ptr x, y, z;
            slong j;
            x = gr_heap_init(ctx);
            y = gr_heap_init(ctx);
            z = gr_heap_init(ctx);
            for (j = 0; j < 10; j++)
            {
                GR_MUST_SUCCEED(gr_randtest(x, state, ctx));
                status = gr_inv(y, x, ctx);
                if (status == GR_SUCCESS)
                {
                    GR_MUST_SUCCEED(gr_mul(z, x, y, ctx));
                    CHECK(gr_is_one(z, ctx) == T_TRUE, "inverse");
                }
                else if (status == GR_DOMAIN)
                    CHECK(gr_is_zero(x, ctx) == T_TRUE, "domain error for a nonzero element");
                else
                    CHECK(status == GR_UNABLE && _check_zero_divisor(ctx, base, n), "zero divisor (inverse)");
            }
            gr_heap_clear(x, ctx);
            gr_heap_clear(y, ctx);
            gr_heap_clear(z, ctx);
        }

        /* polynomial gcd and xgcd */
        {
            gr_poly_t A, B, G, S, T, U, V;
            gr_poly_init(A, ctx);
            gr_poly_init(B, ctx);
            gr_poly_init(G, ctx);
            gr_poly_init(S, ctx);
            gr_poly_init(T, ctx);
            gr_poly_init(U, ctx);
            gr_poly_init(V, ctx);

            _rand_poly(A, state, 1 + n_randint(state, 8), ctx);
            _rand_poly(B, state, 1 + n_randint(state, 8), ctx);
            /* (a common factor, sometimes) */
            if (n_randint(state, 2))
            {
                _rand_poly(U, state, 1 + n_randint(state, 4), ctx);
                GR_MUST_SUCCEED(gr_poly_mul(A, A, U, ctx));
                GR_MUST_SUCCEED(gr_poly_mul(B, B, U, ctx));
            }

            status = gr_poly_xgcd(G, S, T, A, B, ctx);
            if (status == GR_SUCCESS)
            {
                GR_MUST_SUCCEED(gr_poly_mul(U, S, A, ctx));
                GR_MUST_SUCCEED(gr_poly_mul(V, T, B, ctx));
                GR_MUST_SUCCEED(gr_poly_add(U, U, V, ctx));
                CHECK(gr_poly_equal(U, G, ctx) == T_TRUE, "xgcd: G = S A + T B");
            }
            else
                CHECK(status == GR_UNABLE && _check_zero_divisor(ctx, base, n), "zero divisor (xgcd)");

            status = gr_poly_gcd(G, A, B, ctx);
            if (status == GR_SUCCESS && G->length > 0)
            {
                /* G is monic: the remainders are well defined */
                GR_MUST_SUCCEED(gr_poly_divrem(U, V, A, G, ctx));
                CHECK(gr_poly_is_zero(V, ctx) == T_TRUE, "gcd divides A");
                GR_MUST_SUCCEED(gr_poly_divrem(U, V, B, G, ctx));
                CHECK(gr_poly_is_zero(V, ctx) == T_TRUE, "gcd divides B");
            }
            else if (status != GR_SUCCESS)
                CHECK(status == GR_UNABLE && _check_zero_divisor(ctx, base, n), "zero divisor (gcd)");

            /* resultants: against the Sylvester determinant */
            if (A->length >= 1 && B->length >= 1)
            {
                gr_ptr r1, r2;
                r1 = gr_heap_init(ctx);
                r2 = gr_heap_init(ctx);
                status = gr_poly_resultant(r1, A, B, ctx);
                if (status == GR_SUCCESS)
                {
                    GR_MUST_SUCCEED(gr_poly_resultant_sylvester(r2, A, B, ctx));
                    CHECK(gr_equal(r1, r2, ctx) == T_TRUE, "resultant");
                }
                else
                    CHECK(status == GR_UNABLE && _check_zero_divisor(ctx, base, n), "zero divisor (resultant)");
                gr_heap_clear(r1, ctx);
                gr_heap_clear(r2, ctx);
            }

            gr_poly_clear(A, ctx);
            gr_poly_clear(B, ctx);
            gr_poly_clear(G, ctx);
            gr_poly_clear(S, ctx);
            gr_poly_clear(T, ctx);
            gr_poly_clear(U, ctx);
            gr_poly_clear(V, ctx);
        }

        /* linear systems and ranks */
        {
            gr_mat_t A, X, B, AX;
            gr_ptr d;
            slong m = n_randint(state, 5), r;

            gr_mat_init(A, m, m, ctx);
            gr_mat_init(X, m, 2, ctx);
            gr_mat_init(B, m, 2, ctx);
            gr_mat_init(AX, m, 2, ctx);
            d = gr_heap_init(ctx);

            GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));
            GR_MUST_SUCCEED(gr_mat_randtest(B, state, ctx));

            status = gr_mat_nonsingular_solve(X, A, B, ctx);
            if (status == GR_SUCCESS)
            {
                GR_MUST_SUCCEED(gr_mat_mul(AX, A, X, ctx));
                CHECK(gr_mat_equal(AX, B, ctx) == T_TRUE, "solve: A X = B");
            }
            else if (status == GR_DOMAIN)
            {
                /* singular: with unit pivots, the determinant is zero */
                GR_MUST_SUCCEED(gr_mat_det_berkowitz(d, A, ctx));
                CHECK(gr_is_zero(d, ctx) == T_TRUE, "solve: domain error for a nonsingular matrix");
            }
            else
                CHECK(status == GR_UNABLE && _check_zero_divisor(ctx, base, n), "zero divisor (solve)");

            status = gr_mat_rank(&r, A, ctx);
            if (status != GR_SUCCESS)
                CHECK(status == GR_UNABLE && _check_zero_divisor(ctx, base, n), "zero divisor (rank)");

            gr_mat_clear(A, ctx);
            gr_mat_clear(X, ctx);
            gr_mat_clear(B, ctx);
            gr_mat_clear(AX, ctx);
            gr_heap_clear(d, ctx);
        }

        /* turning the pretense off forgets the zero divisor */
        GR_MUST_SUCCEED(gr_ctx_set_is_pretend_field(ctx, T_FALSE));
        CHECK(gr_ctx_is_pretend_field(ctx) != T_TRUE || (fmpz_is_one(p) && base == NULL), "pretense turned off");
        {
            gr_ptr r = gr_heap_init(ctx);
            CHECK(gr_ctx_recover_zero_divisor(r, ctx) == GR_UNABLE, "zero divisor kept after the pretense");
            gr_heap_clear(r, ctx);
        }

        gr_ctx_clear(ctx);
        if (base != NULL)
            gr_ctx_clear(base);
        fmpz_clear(n);
        fmpz_clear(p);
        fmpz_clear(q);
    }

    TEST_FUNCTION_END(state);
}
