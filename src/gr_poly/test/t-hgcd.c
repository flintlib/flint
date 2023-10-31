/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_poly.h"

#define __swap(u, l, v, m) \
  do {                                                  \
    { gr_ptr _; _ = (u), (u) = (v), (v) = _;}           \
    { slong _; _ = (l), (l) = (m), (m) = _;}			\
  } while (0)

#define __mul(C, lenC, A, lenA, B, lenB)                        \
do {                                                            \
    if ((lenA) != 0 && (lenB) != 0)                             \
    {                                                           \
        if ((lenA) >= (lenB))                                   \
            status |= _gr_poly_mul((C), (A), (lenA), (B), (lenB), ctx); \
        else                                                    \
            status |= _gr_poly_mul((C), (B), (lenB), (A), (lenA), ctx); \
        (lenC) = (lenA) + (lenB) - 1;                           \
    }                                                           \
    else                                                        \
    {                                                           \
        (lenC) = 0;                                             \
    }                                                           \
} while (0)

TEST_FUNCTION_START(gr_poly_hgcd, state)
{
    int i, result;

    /*
       Find coprime polys, multiply by another poly
       and check the GCD is that poly
     */
    for (i = 0; i < 1000; i++)
    {
        gr_poly_t a, b, c, d, c1, d1, s, t;
        gr_ctx_t ctx;
        gr_ptr M[4];
        slong lenM[4];
        slong sgnM;
        slong cutoff;
        int status = GR_SUCCESS;
        slong sz;

        while (1)
        {
            gr_ctx_init_random(ctx, state);

            if (gr_ctx_is_finite(ctx) == T_TRUE)
                break;
            else
                gr_ctx_clear(ctx);
        }

        sz = ctx->sizeof_elem;

        gr_poly_init(a, ctx);
        gr_poly_init(b, ctx);
        gr_poly_init(c, ctx);
        gr_poly_init(d, ctx);
        gr_poly_init(c1, ctx);
        gr_poly_init(d1, ctx);
        gr_poly_init(s, ctx);
        gr_poly_init(t, ctx);

        do
        {
            status |= gr_poly_randtest(a, state, n_randint(state, 100) + 1, ctx);
            status |= gr_poly_randtest(b, state, n_randint(state, 100) + 1, ctx);

            /* escape the zero ring */
            GR_IGNORE(gr_poly_one(c, ctx));
            if (gr_poly_is_zero(c, ctx) != T_FALSE)
                break;

        } while (a->length == 0 || b->length == 0 || a->length == b->length);

        cutoff = n_randint(state, 100);

        if (a->length < b->length)
            gr_poly_swap(a, b, ctx);

        M[0] = flint_malloc(a->length * sz);
        M[1] = flint_malloc(a->length * sz);
        M[2] = flint_malloc(a->length * sz);
        M[3] = flint_malloc(a->length * sz);

        _gr_vec_init(M[0], a->length, ctx);
        _gr_vec_init(M[1], a->length, ctx);
        _gr_vec_init(M[2], a->length, ctx);
        _gr_vec_init(M[3], a->length, ctx);

        gr_poly_fit_length(c, a->length, ctx);
        gr_poly_fit_length(d, b->length, ctx);

        status |= _gr_poly_hgcd(NULL, &sgnM, M, lenM,
                                        c->coeffs, &(c->length), d->coeffs,
                                        &(d->length), a->coeffs, a->length,
                                        b->coeffs, b->length, cutoff, ctx);

        gr_poly_fit_length(s, 2 * a->length, ctx);
        gr_poly_fit_length(t, 2 * a->length, ctx);

        /* [c1,d1] := sgnM * M^{-1} [a,b] */
        {
            __swap(M[0], lenM[0], M[3], lenM[3]);
            status |= _gr_vec_neg(M[1], M[1], lenM[1], ctx);
            status |= _gr_vec_neg(M[2], M[2], lenM[2], ctx);

            __mul(s->coeffs, s->length, M[0], lenM[0], a->coeffs, a->length);
            __mul(t->coeffs, t->length, M[1], lenM[1], b->coeffs, b->length);
            status |= gr_poly_add(c1, s, t, ctx);
            __mul(s->coeffs, s->length, M[2], lenM[2], a->coeffs, a->length);
            __mul(t->coeffs, t->length, M[3], lenM[3], b->coeffs, b->length);
            status |= gr_poly_add(d1, s, t, ctx);
        }

        if (sgnM < 0)
        {
            status |= gr_poly_neg(c1, c1, ctx);
            status |= gr_poly_neg(d1, d1, ctx);
        }

        if (status != GR_SUCCESS && gr_ctx_is_field(ctx) == T_TRUE)
        {
            flint_printf("FAIL: unexpected status %d\n", status);
            gr_ctx_println(ctx);
            flint_printf("a  = "), gr_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b  = "), gr_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("c  = "), gr_poly_print(c, ctx), flint_printf("\n\n");
            flint_printf("d  = "), gr_poly_print(d, ctx), flint_printf("\n\n");
            flint_printf("c1  = "), gr_poly_print(c1, ctx), flint_printf("\n\n");
            flint_printf("d1  = "), gr_poly_print(d1, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        if (status == GR_SUCCESS)
        {
            result = (gr_poly_equal(c, c1, ctx) == T_TRUE && gr_poly_equal(d, d1, ctx) == T_TRUE);
            if (!result)
            {
                flint_printf("FAIL:\n");
                flint_printf("a  = "), gr_poly_print(a, ctx), flint_printf("\n\n");
                flint_printf("b  = "), gr_poly_print(b, ctx), flint_printf("\n\n");
                flint_printf("c  = "), gr_poly_print(c, ctx), flint_printf("\n\n");
                flint_printf("d  = "), gr_poly_print(d, ctx), flint_printf("\n\n");
                flint_printf("c1  = "), gr_poly_print(c1, ctx), flint_printf("\n\n");
                flint_printf("d1  = "), gr_poly_print(d1, ctx), flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }

        _gr_vec_clear(M[0], a->length, ctx);
        _gr_vec_clear(M[1], a->length, ctx);
        _gr_vec_clear(M[2], a->length, ctx);
        _gr_vec_clear(M[3], a->length, ctx);

        flint_free(M[0]);
        flint_free(M[1]);
        flint_free(M[2]);
        flint_free(M[3]);

        gr_poly_clear(a, ctx);
        gr_poly_clear(b, ctx);
        gr_poly_clear(c, ctx);
        gr_poly_clear(d, ctx);
        gr_poly_clear(c1, ctx);
        gr_poly_clear(d1, ctx);
        gr_poly_clear(s, ctx);
        gr_poly_clear(t, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
