/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_EC_TEST_HELPERS_H
#define GR_EC_TEST_HELPERS_H

#include "ulong_extras.h"
#include "fmpz.h"
#include "gr.h"
#include "gr_ec.h"

/* A base ring over which elliptic curves make sense. Includes rings that
   are not fields, and rings of characteristic 2 and 3, so that the
   rejection and GR_UNABLE paths get exercised too. */
static inline void
gr_ec_test_ring(gr_ctx_t ctx, flint_rand_t state)
{
    fmpz_t p;

    switch (n_randint(state, 7))
    {
        case 0:
            gr_ctx_init_fmpq(ctx);
            break;

        case 1:
            gr_ctx_init_fmpz(ctx);
            break;

        case 2:
        case 3:
            while (gr_ctx_init_nmod(ctx,
                        n_randprime(state, 4 + n_randint(state, 20), 1))
                    != GR_SUCCESS)
                ;
            break;

        case 4:
            fmpz_init(p);
            fmpz_set_ui(p, n_randprime(state, 4 + n_randint(state, 40), 1));
            gr_ctx_init_fmpz_mod(ctx, p);
            fmpz_clear(p);
            break;

        case 5:
            gr_ctx_init_fq_nmod(ctx, n_randprime(state, 3 + n_randint(state, 4), 1),
                    1 + n_randint(state, 3), "a");
            break;

        default:
            /* characteristic 2 and 3 */
            gr_ctx_init_fq_nmod(ctx, 2 + n_randint(state, 2),
                    1 + n_randint(state, 4), "a");
            break;
    }
}

/* An inexact ring, where zero tests are undecidable and essentially every
   operation is allowed to come back GR_UNABLE. */
static inline void
gr_ec_test_ring_inexact(gr_ctx_t ctx, flint_rand_t state)
{
    if (n_randint(state, 2))
        gr_ctx_init_real_arb(ctx, 64);
    else
        gr_ctx_init_complex_acb(ctx, 64);
}

/* A random curve. Both models are generated: short Weierstrass most of
   the time, and an explicit long Weierstrass equation often enough that
   the general formulas get exercised. Returns the gr status flag. */
static inline int
gr_ec_test_curve(gr_ec_ctx_t E, gr_ctx_t R, flint_rand_t state)
{
    slong which = n_randint(state, 4);
    gr_ptr a;
    int status = GR_UNABLE;
    slong iter;

    if (which == 0)
        return gr_ec_ctx_init_randtest(E, state, R);

    a = gr_heap_init_vec(5, R);

    for (iter = 0; iter < 10; iter++)
    {
        slong i;
        slong lo = (which == 1) ? 0 : 3;

        status = GR_SUCCESS;

        for (i = 0; i < 3; i++)
            status |= (i >= lo) ? gr_randtest(GR_ENTRY(a, i, R->sizeof_elem), state, R)
                                : gr_zero(GR_ENTRY(a, i, R->sizeof_elem), R);

        for (i = 3; i < 5; i++)
            status |= gr_randtest(GR_ENTRY(a, i, R->sizeof_elem), state, R);

        if (status == GR_SUCCESS)
            status = gr_ec_ctx_init(E, R,
                        GR_ENTRY(a, 0, R->sizeof_elem),
                        GR_ENTRY(a, 1, R->sizeof_elem),
                        GR_ENTRY(a, 2, R->sizeof_elem),
                        GR_ENTRY(a, 3, R->sizeof_elem),
                        GR_ENTRY(a, 4, R->sizeof_elem));

        if (status == GR_SUCCESS)
            break;
    }

    gr_heap_clear_vec(a, 5, R);

    return status;
}

/* Scalars have to stay tiny over infinite rings: without a normalization
   step the coordinates of 2^k P grow like 4^k, so even a 6-bit scalar is
   already expensive over ZZ or QQ. */
static inline slong
gr_ec_test_scalar_bits(gr_ctx_t R)
{
    return (gr_ctx_is_finite(R) == T_TRUE) ? 60 : 2;
}

#endif
