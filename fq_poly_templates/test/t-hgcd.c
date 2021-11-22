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

#ifdef T

#include "templates.h"

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "mpn_extras.h"

#define __mul(C, lenC, A, lenA, B, lenB)                        \
do {                                                            \
    if ((lenA) != 0 && (lenB) != 0)                             \
    {                                                           \
        if ((lenA) >= (lenB))                                   \
            _TEMPLATE(T, poly_mul)((C), (A), (lenA), (B), (lenB), ctx); \
        else                                                    \
            _TEMPLATE(T, poly_mul)((C), (B), (lenB), (A), (lenA), ctx); \
        (lenC) = (lenA) + (lenB) - 1;                           \
    }                                                           \
    else                                                        \
    {                                                           \
        (lenC) = 0;                                             \
    }                                                           \
} while (0)


#define __swap(u, l, v, m) \
  do {								\
    { TEMPLATE(T, struct)* _; _ = (u), (u) = (v), (v) = _;}     \
    { slong _; _ = (l), (l) = (m), (m) = _;}			\
  } while (0)

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    flint_printf("hgcd....");
    fflush(stdout);

    /* 
       Find coprime polys, multiply by another poly 
       and check the GCD is that poly 
     */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) a, b, c, d, c1, d1, s, t;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, struct) * M[4];
        slong lenM[4];
        slong sgnM;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (d, ctx);
        TEMPLATE(T, poly_init) (c1, ctx);
        TEMPLATE(T, poly_init) (d1, ctx);
        TEMPLATE(T, poly_init) (s, ctx);
        TEMPLATE(T, poly_init) (t, ctx);

        do
        {
            TEMPLATE(T, poly_randtest_not_zero) (a, state,
                                                 n_randint(state, 800) + 1,
                                                 ctx);
            TEMPLATE(T, poly_randtest_not_zero) (b, state,
                                                 n_randint(state, 800) + 1,
                                                 ctx);
        } while (a->length == b->length);

        if (a->length < b->length)
            TEMPLATE(T, poly_swap) (a, b, ctx);

        M[0] = _TEMPLATE(T, vec_init) (a->length, ctx);
        M[1] = _TEMPLATE(T, vec_init) (a->length, ctx);
        M[2] = _TEMPLATE(T, vec_init) (a->length, ctx);
        M[3] = _TEMPLATE(T, vec_init) (a->length, ctx);

        TEMPLATE(T, poly_fit_length) (c, a->length, ctx);
        TEMPLATE(T, poly_fit_length) (d, b->length, ctx);

        sgnM = _TEMPLATE(T, poly_hgcd) (M, lenM,
                                        c->coeffs, &(c->length), d->coeffs,
                                        &(d->length), a->coeffs, a->length,
                                        b->coeffs, b->length, ctx);

        TEMPLATE(T, poly_fit_length) (s, 2 * a->length, ctx);
        TEMPLATE(T, poly_fit_length) (t, 2 * a->length, ctx);

        /* [c1,d1] := sgnM * M^{-1} [a,b] */
        {
            __swap(M[0], lenM[0], M[3], lenM[3]);
            _TEMPLATE(T, vec_neg) (M[1], M[1], lenM[1], ctx);
            _TEMPLATE(T, vec_neg) (M[2], M[2], lenM[2], ctx);

            __mul(s->coeffs, s->length, M[0], lenM[0], a->coeffs, a->length);
            __mul(t->coeffs, t->length, M[1], lenM[1], b->coeffs, b->length);
            TEMPLATE(T, poly_add) (c1, s, t, ctx);
            __mul(s->coeffs, s->length, M[2], lenM[2], a->coeffs, a->length);
            __mul(t->coeffs, t->length, M[3], lenM[3], b->coeffs, b->length);
            TEMPLATE(T, poly_add) (d1, s, t, ctx);
        }

        if (sgnM < 0)
        {
            TEMPLATE(T, poly_neg) (c1, c1, ctx);
            TEMPLATE(T, poly_neg) (d1, d1, ctx);
        }

        result = (TEMPLATE(T, poly_equal) (c, c1, ctx)
                  && TEMPLATE(T, poly_equal) (d, d1, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a  = "), TEMPLATE(T, poly_print_pretty) (a, "x",
                                                                   ctx),
                flint_printf("\n\n");
            flint_printf("b  = "), TEMPLATE(T, poly_print_pretty) (b, "x",
                                                                   ctx),
                flint_printf("\n\n");
            flint_printf("c  = "), TEMPLATE(T, poly_print_pretty) (c, "x",
                                                                   ctx),
                flint_printf("\n\n");
            flint_printf("d  = "), TEMPLATE(T, poly_print_pretty) (d, "x",
                                                                   ctx),
                flint_printf("\n\n");
            flint_printf("c1 = "), TEMPLATE(T, poly_print_pretty) (c1, "x",
                                                                   ctx),
                flint_printf("\n\n");
            flint_printf("d1 = "), TEMPLATE(T, poly_print_pretty) (d1, "x",
                                                                   ctx),
                flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (d, ctx);
        TEMPLATE(T, poly_clear) (c1, ctx);
        TEMPLATE(T, poly_clear) (d1, ctx);
        TEMPLATE(T, poly_clear) (s, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);

        _TEMPLATE(T, vec_clear) (M[0], a->length, ctx);
        _TEMPLATE(T, vec_clear) (M[1], a->length, ctx);
        _TEMPLATE(T, vec_clear) (M[2], a->length, ctx);
        _TEMPLATE(T, vec_clear) (M[3], a->length, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    flint_randclear(state);

    flint_printf("PASS\n");
    return 0;
}

#undef __mul
#undef __swap

#endif
