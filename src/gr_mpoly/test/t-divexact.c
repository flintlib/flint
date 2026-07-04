/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "gr_mpoly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_mpoly_divexact, state)
{
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t f, g, h, q1, q2, q3, q4;
        slong lenf, leng;
        flint_bitcnt_t exp_bits;
        int status, s1, s2, s3, s4;
        int monomial_case = (n_randint(state, 4) == 0);
        int big = (!monomial_case && n_randint(state, 8) == 0);

        switch (n_randint(state, 4))
        {
            case 0:
                gr_ctx_init_random_finite_field(cctx, state);
                break;
            case 1:
                gr_ctx_init_fmpz(cctx);
                break;
            case 2:
                gr_ctx_init_fmpq(cctx);
                break;
            default:
                gr_ctx_init_nmod(cctx, n_randtest_prime(state, 1));
                break;
        }

        gr_mpoly_ctx_init_rand(ctx, state, cctx, 4);

        gr_mpoly_init(f, ctx);
        gr_mpoly_init(g, ctx);
        gr_mpoly_init(h, ctx);
        gr_mpoly_init(q1, ctx);
        gr_mpoly_init(q2, ctx);
        gr_mpoly_init(q3, ctx);
        gr_mpoly_init(q4, ctx);

        flint_set_num_threads(1 + n_randint(state, 4));

        /* occasionally big enough to actually exercise the threaded
           dispatch inside gr_mpoly_divexact (A->length > 500, B->length
           > 2), and always call gr_mpoly_divexact_heap_threaded directly
           too, regardless of size, to exercise it even when small inputs
           would otherwise make gr_mpoly_divexact stay serial */
        lenf = n_randint(state, big ? 400 : 20) + 1;
        leng = monomial_case ? 1 : (n_randint(state, big ? 80 : 10) + 3);
        exp_bits = n_randint(state, big ? 12 : 4) + 2;

        status = GR_SUCCESS;
        status |= gr_mpoly_randtest_bits(f, state, lenf, exp_bits, ctx);
        status |= gr_mpoly_randtest_bits(g, state, leng, exp_bits, ctx);
        if (gr_mpoly_is_zero(g, ctx) != T_FALSE)
            status |= gr_mpoly_one(g, ctx);

        if (status == GR_SUCCESS)
        {
            /* h = f * g, so g divides h exactly with quotient f */
            status = gr_mpoly_mul(h, f, g, ctx);

            if (status == GR_SUCCESS)
            {
                s1 = gr_mpoly_divides(q1, h, g, ctx);
                s2 = gr_mpoly_div(q2, h, g, ctx);
                s3 = gr_mpoly_divexact(q3, h, g, ctx);
                if (monomial_case)
                {
                    s4 = s3;
                    GR_IGNORE(gr_mpoly_set(q4, q3, ctx));
                }
                else
                {
                    s4 = gr_mpoly_divexact_heap_threaded(q4, h, g, ctx);
                }

                if (s1 == GR_DOMAIN)
                {
                    /* g always divides h = f*g exactly */
                    flint_printf("FAIL: exact division reported GR_DOMAIN "
                                 "by gr_mpoly_divides\n");
                    flint_printf("i = %wd\n", i);
                    gr_ctx_println(cctx);
                    fflush(stdout);
                    flint_abort();
                }

                if (s1 == GR_SUCCESS &&
                    (s2 != GR_SUCCESS || s3 != GR_SUCCESS || s4 != GR_SUCCESS))
                {
                    flint_printf("FAIL: exact division but "
                                 "div = %d, divexact = %d, "
                                 "divexact_heap_threaded = %d\n", s2, s3, s4);
                    flint_printf("i = %wd\n", i);
                    gr_ctx_println(cctx);
                    fflush(stdout);
                    flint_abort();
                }

                if (s1 == GR_SUCCESS)
                {
                    gr_mpoly_assert_canonical(q1, ctx);
                    gr_mpoly_assert_canonical(q2, ctx);
                    gr_mpoly_assert_canonical(q3, ctx);
                    gr_mpoly_assert_canonical(q4, ctx);

                    if (gr_mpoly_equal(q1, q2, ctx) == T_FALSE ||
                        gr_mpoly_equal(q1, q3, ctx) == T_FALSE ||
                        gr_mpoly_equal(q1, q4, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: quotient mismatch "
                                     "(monomial_case = %d, big = %d)\n",
                                     monomial_case, big);
                        flint_printf("i = %wd\n", i);
                        gr_ctx_println(cctx);
                        flint_printf("divides               = "); gr_mpoly_print_pretty(q1, ctx); flint_printf("\n");
                        flint_printf("div                   = "); gr_mpoly_print_pretty(q2, ctx); flint_printf("\n");
                        flint_printf("divexact              = "); gr_mpoly_print_pretty(q3, ctx); flint_printf("\n");
                        flint_printf("divexact_heap_threaded = "); gr_mpoly_print_pretty(q4, ctx); flint_printf("\n");
                        fflush(stdout);
                        flint_abort();
                    }

                    if (gr_mpoly_equal(q1, f, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: quotient != f\n");
                        flint_printf("i = %wd\n", i);
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }
        }

        gr_mpoly_clear(f, ctx);
        gr_mpoly_clear(g, ctx);
        gr_mpoly_clear(h, ctx);
        gr_mpoly_clear(q1, ctx);
        gr_mpoly_clear(q2, ctx);
        gr_mpoly_clear(q3, ctx);
        gr_mpoly_clear(q4, ctx);
        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    /*
        Regression test: divexact against a monomial divisor must not
        confuse "not divisible by the monomial" with "not exact". If B is
        a single term x^beta, and A is genuinely divisible by B, then
        *every* term of A must be divisible by x^beta -- a term that is
        not can only be a spurious near-zero artifact of an inexact ring
        (there is no possibility of cancellation between distinct terms
        of the quotient when multiplying by a single monomial). Such a
        term must be silently dropped, not reported as GR_DOMAIN or
        GR_UNABLE, exactly like the discarded remainder of the general
        (non-monomial) case.
    */
    {
        gr_ctx_t R;
        gr_mpoly_ctx_t Rx;
        gr_mpoly_t A, B, Q;
        ulong exp[2];
        gr_ptr c;
        int status;

        gr_ctx_init_real_arb(R, 64);
        gr_mpoly_ctx_init(Rx, R, 2, ORD_LEX);

        gr_mpoly_init(A, Rx);
        gr_mpoly_init(B, Rx);
        gr_mpoly_init(Q, Rx);

        /* A = ([+/- eps]) * x^3 + y^2, a ball straddling zero but not
           provably zero */
        GR_TMP_INIT(c, R);

        exp[0] = 3; exp[1] = 0;
        status = gr_zero(c, R);
        arb_add_error_2exp_si((arb_struct *) c, -3 - n_randint(state, 20));
        status |= gr_mpoly_push_term_scalar_ui(A, c, exp, Rx);

        exp[0] = 0; exp[1] = 2;
        status |= gr_one(c, R);
        status |= gr_mpoly_push_term_scalar_ui(A, c, exp, Rx);

        GR_TMP_CLEAR(c, R);

        gr_mpoly_sort_terms(A, Rx);
        status |= gr_mpoly_combine_like_terms(A, Rx);

        /* B = y */
        status |= gr_mpoly_gen(B, 1, Rx);

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: unexpected failure setting up ball test\n");
            fflush(stdout);
            flint_abort();
        }

        status = gr_mpoly_divexact(Q, A, B, Rx);

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: gr_mpoly_divexact on ball input with a "
                         "spurious non-dividing term returned %d "
                         "(expected GR_SUCCESS)\n", status);
            flint_printf("A = "); gr_mpoly_print_pretty(A, Rx); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        if (gr_mpoly_equal(Q, B, Rx) == T_FALSE)
        {
            flint_printf("FAIL: expected quotient y, got ");
            gr_mpoly_print_pretty(Q, Rx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        gr_mpoly_clear(A, Rx);
        gr_mpoly_clear(B, Rx);
        gr_mpoly_clear(Q, Rx);
        gr_mpoly_ctx_clear(Rx);
        gr_ctx_clear(R);
    }

    /*
        Regression test: an "unknown" gr_is_zero verdict on the running
        accumulator, while reducing a term whose monomial *does* divide
        that of the (non-monomial) divisor, must not abort gr_mpoly_divexact
        or gr_mpoly_divrem with GR_UNABLE. We only need the coefficient
        division itself to succeed; an accumulator that cannot be proven
        zero is simply kept as a (possibly redundant, wide) term instead
        of being dropped, exactly like gr_poly_divrem already does for
        univariate polynomials.

        Concretely: a = y, b = (1/3)*x + 1, c = a*b = (1/3)*x*y + y.
        Dividing c by b (whose leading coefficient 1/3 is a unit but not
        exactly representable) is expected to succeed with quotient y,
        even though intermediate arb accumulators along the way straddle
        zero without being provably zero.
    */
    {
        gr_ctx_t R;
        gr_mpoly_ctx_t Rx;
        gr_mpoly_t a, b, c, q, q2, r2;
        int status;

        gr_ctx_init_real_arb(R, 64);
        gr_mpoly_ctx_init(Rx, R, 2, ORD_LEX);

        {
            const char * ss[] = {"x", "y"};
            GR_MUST_SUCCEED(gr_ctx_set_gen_names(Rx, ss));
        }

        gr_mpoly_init(a, Rx);
        gr_mpoly_init(b, Rx);
        gr_mpoly_init(c, Rx);
        gr_mpoly_init(q, Rx);
        gr_mpoly_init(q2, Rx);
        gr_mpoly_init(r2, Rx);

        status = GR_SUCCESS;
        status |= gr_mpoly_gen(a, 1, Rx);                    /* a = y */
        status |= gr_set_str(b, "(1/3)*x + 1", Rx);          /* b = (1/3)x + 1 */
        status |= gr_mpoly_mul(c, a, b, Rx);                 /* c = a*b */

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: unexpected failure setting up divexact/arb test\n");
            fflush(stdout);
            flint_abort();
        }

        status = gr_mpoly_divexact(q, c, b, Rx);
        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: gr_mpoly_divexact(a*b, b) over arb returned "
                         "%d (expected GR_SUCCESS)\n", status);
            fflush(stdout);
            flint_abort();
        }
        if (gr_mpoly_equal(q, a, Rx) == T_FALSE)
        {
            flint_printf("FAIL: gr_mpoly_divexact(a*b, b) over arb gave a "
                         "quotient inconsistent with a = y: ");
            gr_mpoly_print_pretty(q, Rx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        status = gr_mpoly_divrem(q2, r2, c, b, Rx);
        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: gr_mpoly_divrem(a*b, b) over arb returned "
                         "%d (expected GR_SUCCESS)\n", status);
            fflush(stdout);
            flint_abort();
        }
        if (gr_mpoly_equal(q2, a, Rx) == T_FALSE)
        {
            flint_printf("FAIL: gr_mpoly_divrem(a*b, b) over arb gave a "
                         "quotient inconsistent with a = y: ");
            gr_mpoly_print_pretty(q2, Rx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        /* by contrast, gr_mpoly_divides is allowed (expected) to report
           GR_UNABLE here, since exactness genuinely cannot be certified
           over an inexact ring; this is a deliberate difference in
           semantics from gr_mpoly_div(exact)/gr_mpoly_divrem, not a bug */
        {
            gr_mpoly_t qd;
            gr_mpoly_init(qd, Rx);
            status = gr_mpoly_divides(qd, c, b, Rx);
            if (status == GR_DOMAIN)
            {
                flint_printf("FAIL: gr_mpoly_divides(a*b, b) over arb "
                             "incorrectly reported GR_DOMAIN\n");
                fflush(stdout);
                flint_abort();
            }
            gr_mpoly_clear(qd, Rx);
        }

        gr_mpoly_clear(a, Rx);
        gr_mpoly_clear(b, Rx);
        gr_mpoly_clear(c, Rx);
        gr_mpoly_clear(q, Rx);
        gr_mpoly_clear(q2, Rx);
        gr_mpoly_clear(r2, Rx);
        gr_mpoly_ctx_clear(Rx);
        gr_ctx_clear(R);
    }

    TEST_FUNCTION_END(state);
}
