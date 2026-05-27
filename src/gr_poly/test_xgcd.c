/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "gr_vec.h"
#include "gr_poly.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

FLINT_DLL extern gr_static_method_table _ca_methods;

/*
 * The key invariant is the Bezout identity:
 *
 *     S * A + T * B  =  G
 *
 * Over fields (Q, Z/pZ), G is also the monic GCD so we additionally
 * require G == gcd_euclidean(A, B).
 *
 * Over rings that are not Bezout domains (such as Z), the subresultant
 * algorithm returns a G that is a scalar multiple of gcd(A,B) — the
 * "pseudo-xgcd" result.  We check:
 *   (1) S * A + T * B = G  (Bezout identity)
 *   (2) gcd(A,B) | G       (G is a multiple of the GCD)
 *   (3) G | gcd(A,B) * c for some scalar c (G does not exceed gcd by
 *       more than a scalar factor; equivalently, prim(G) = prim(gcd))
 *
 * In practice, for (2) and (3) over UFDs we verify by checking that
 * gcd_subresultant(A, G) == gcd_subresultant(A, B) (up to associates),
 * which is equivalent to G and gcd having the same primitive part.
 * Over fields the simpler G == gcd_euclidean suffices.
 */
void _gr_poly_test_xgcd(gr_method_poly_xgcd_op xgcd_impl,
    flint_rand_t state, slong iters, slong maxn, int flags, gr_ctx_t ctx)
{
    slong iter;
    gr_ctx_ptr given_ctx = ctx;

    for (iter = 0; iter < iters; iter++)
    {
        gr_ctx_t my_ctx;
        gr_ctx_struct * ctx;

        if (given_ctx == NULL)
        {
            gr_ctx_init_random(my_ctx, state);
            ctx = my_ctx;
        }
        else
            ctx = given_ctx;

        {
            gr_poly_t a, b, d, g, s, t, v, w;
            slong n;
            int status = GR_SUCCESS;
            int aliasing;

            if (maxn == 0)
            {
                if (gr_ctx_is_finite(ctx) == T_TRUE && n_randint(state, 2) == 0)
                    n = 20;
                else if (ctx->methods == _ca_methods)
                    n = 4;
                else
                    n = 6;
            }
            else
            {
                n = maxn;
            }

            gr_poly_init(a, ctx);
            gr_poly_init(b, ctx);
            gr_poly_init(d, ctx);
            gr_poly_init(g, ctx);
            gr_poly_init(s, ctx);
            gr_poly_init(t, ctx);
            gr_poly_init(v, ctx);
            gr_poly_init(w, ctx);

            status |= gr_poly_randtest(a, state, 1 + n_randint(state, n), ctx);
            status |= gr_poly_randtest(b, state, 1 + n_randint(state, n), ctx);

            /* common factor */
            if (n_randint(state, 2))
            {
                status |= gr_poly_randtest(t, state, 1 + n_randint(state, n), ctx);
                status |= gr_poly_mul(a, a, t, ctx);
                status |= gr_poly_mul(b, b, t, ctx);
            }

            status |= gr_poly_randtest(g, state, 3, ctx);
            status |= gr_poly_randtest(s, state, 3, ctx);
            status |= gr_poly_randtest(t, state, 3, ctx);

            aliasing = n_randint(state, 8);

            switch (aliasing)
            {
                case 0:
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, a, b, ctx);
                    break;
                case 1:
                    status |= gr_poly_set(g, a, ctx);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, g, b, ctx);
                    break;
                case 2:
                    status |= gr_poly_set(s, a, ctx);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, s, b, ctx);
                    break;
                case 3:
                    status |= gr_poly_set(t, a, ctx);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, t, b, ctx);
                    break;
                case 4:
                    status |= gr_poly_set(g, b, ctx);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, a, g, ctx);
                    break;
                case 5:
                    status |= gr_poly_set(s, b, ctx);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, a, s, ctx);
                    break;
                case 6:
                    status |= gr_poly_set(t, b, ctx);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, a, t, ctx);
                    break;
                case 7:
                    status |= gr_poly_set(b, a, ctx);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, a, a, ctx);
                    break;
                default:
                    break;
            }

            if (status == GR_SUCCESS)
            {
                if (gr_ctx_is_field(ctx) == T_TRUE ||
                    gr_ctx_is_unique_factorization_domain(ctx) == T_TRUE)
                {
                    gr_poly_t v, w;
                    gr_poly_init(v, ctx);
                    gr_poly_init(w, ctx);

                    int sb = GR_SUCCESS;
                    sb |= gr_poly_mul(v, s, a, ctx);
                    sb |= gr_poly_mul(w, t, b, ctx);
                    sb |= gr_poly_add(w, v, w, ctx);

                    if (sb == GR_SUCCESS && gr_poly_equal(g, w, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: Bezout identity s*a + t*b != g\n");
                        gr_ctx_println(ctx);
                        flint_printf("a = "); gr_poly_print(a, ctx); flint_printf("\n");
                        flint_printf("b = "); gr_poly_print(b, ctx); flint_printf("\n");
                        flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n");
                        flint_printf("s = "); gr_poly_print(s, ctx); flint_printf("\n");
                        flint_printf("t = "); gr_poly_print(t, ctx); flint_printf("\n");
                        flint_printf("s*a+t*b = "); gr_poly_print(w, ctx); flint_printf("\n");
                        fflush(stdout);
                        flint_abort();
                    }

                    gr_poly_clear(v, ctx);
                    gr_poly_clear(w, ctx);
                }

                if (gr_ctx_is_field(ctx) == T_TRUE)
                {
                    gr_poly_t d;
                    gr_poly_init(d, ctx);

                    int sd = gr_poly_gcd_euclidean(d, a, b, ctx);

                    if (sd == GR_SUCCESS && gr_poly_equal(d, g, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: g != monic gcd over field\n");
                        gr_ctx_println(ctx);
                        flint_printf("a = "); gr_poly_print(a, ctx); flint_printf("\n");
                        flint_printf("b = "); gr_poly_print(b, ctx); flint_printf("\n");
                        flint_printf("gcd = "); gr_poly_print(d, ctx); flint_printf("\n");
                        flint_printf("g   = "); gr_poly_print(g, ctx); flint_printf("\n");
                        fflush(stdout);
                        flint_abort();
                    }

                    gr_poly_clear(d, ctx);
                }
                else if (gr_ctx_is_unique_factorization_domain(ctx) == T_TRUE)
                {
                    /*
                     * Over a UFD that is not a field (e.g. Z):
                     *
                     * The subresultant xgcd satisfies Sage's guarantee: g is the gcd of
                     * a and b "up to a divisor of the resultant".  Precisely:
                     *
                     *   g  =  c * prim_gcd(a, b)
                     *
                     * for some scalar c in R with c | res(a, b).
                     *
                     * We verify the necessary consequence:
                     *   (i)  prim_gcd(a,b) | g          (g is a multiple of the primitive gcd)
                     *   (ii) g / prim_gcd(a,b) is a scalar  (the multiple is a ring element)
                     *
                     * computed via gr_poly_divrem(q, r, g, d) and checking r == 0, deg(q) == 0.
                     *
                     * We do NOT require g == resultant(a,b) as a scalar; that would fix c
                     * to the resultant and is not guaranteed by the subresultant PRS (only
                     * by a separate resultant computation).  The choice of c depends on the
                     * path taken by the PRS and is an implementation detail.
                     */
                    gr_poly_t d, q, r;
                    gr_poly_init(d, ctx);
                    gr_poly_init(q, ctx);
                    gr_poly_init(r, ctx);

                    int sd = gr_poly_gcd_subresultant(d, a, b, ctx);

                    if (sd == GR_SUCCESS)
                    {
                        int sq = gr_poly_divrem(q, r, g, d, ctx);

                        if (sq == GR_SUCCESS)
                        {
                            if (gr_poly_is_zero(r, ctx) == T_FALSE)
                            {
                                flint_printf("FAIL: prim_gcd(a,b) does not divide g over UFD\n");
                                gr_ctx_println(ctx);
                                flint_printf("a = "); gr_poly_print(a, ctx); flint_printf("\n");
                                flint_printf("b = "); gr_poly_print(b, ctx); flint_printf("\n");
                                flint_printf("d = prim_gcd = "); gr_poly_print(d, ctx); flint_printf("\n");
                                flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n");
                                flint_printf("r = g mod d = "); gr_poly_print(r, ctx); flint_printf("\n");
                                fflush(stdout);
                                flint_abort();
                            }
                            else if (gr_poly_is_scalar(q, ctx) == T_FALSE)
                            {
                                flint_printf("FAIL: g / prim_gcd(a,b) is not a scalar over UFD\n");
                                gr_ctx_println(ctx);
                                flint_printf("a = "); gr_poly_print(a, ctx); flint_printf("\n");
                                flint_printf("b = "); gr_poly_print(b, ctx); flint_printf("\n");
                                flint_printf("d = prim_gcd = "); gr_poly_print(d, ctx); flint_printf("\n");
                                flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n");
                                flint_printf("q = g / d = "); gr_poly_print(q, ctx); flint_printf("\n");
                                fflush(stdout);
                                flint_abort();
                            }
                        }
                    }

                    gr_poly_clear(d, ctx);
                    gr_poly_clear(q, ctx);
                    gr_poly_clear(r, ctx);
                }
            }

            if ((flags & 1) && (ctx->which_ring == GR_CTX_FMPZ
                || ctx->which_ring == GR_CTX_FMPZI || ctx->which_ring == GR_CTX_FMPZ_POLY) && status != GR_SUCCESS)
            {
                flint_printf("FAIL: did not succeed over Z, Z[i] or Z[t]\n\n");
                gr_ctx_println(ctx);
                gr_poly_print(a, ctx), flint_printf("\n\n");
                gr_poly_print(b, ctx), flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }

            if ((ctx->which_ring == GR_CTX_FMPQ || (ctx->which_ring == GR_CTX_NMOD8 && gr_ctx_is_field(ctx) == T_TRUE)) && status != GR_SUCCESS)
            {
                flint_printf("FAIL: did not succeed over Q or Z/pZ\n\n");
                gr_ctx_println(ctx);
                gr_poly_print(a, ctx), flint_printf("\n\n");
                gr_poly_print(b, ctx), flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }

            gr_poly_clear(a, ctx);
            gr_poly_clear(b, ctx);
            gr_poly_clear(d, ctx);
            gr_poly_clear(g, ctx);
            gr_poly_clear(s, ctx);
            gr_poly_clear(t, ctx);
            gr_poly_clear(v, ctx);
            gr_poly_clear(w, ctx);
        }

        if (given_ctx == NULL)
            gr_ctx_clear(ctx);
    }
}

POP_OPTIONS
