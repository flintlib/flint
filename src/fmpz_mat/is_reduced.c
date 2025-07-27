/*
    Copyright (C) 2014 Abhinav Baid
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_mat.h"
#include "fmpz_mat.h"

int
fmpz_mat_is_reduced(const fmpz_mat_t A, double delta, double eta)
{
    truth_t is_reduced = T_UNKNOWN;
    slong prec;
    gr_ctx_t ctx;
    gr_ptr Rdelta, Reta;
    gr_mat_t RA;
    slong exact_cutoff;

    /* To do: for very small matrices, consider doing a division-free version
       of the naive algorithm over Z instead of working over Q. */
    /* To do: this is not at all tuned. */
    exact_cutoff = fmpz_mat_max_bits(A);
    exact_cutoff = FLINT_ABS(exact_cutoff);
    exact_cutoff = (64 + exact_cutoff) * FLINT_MAX(A->r, A->c);

    for (prec = 64; ; prec *= 2)
    {
        /* flint_printf("fmpz_mat_is_reduced : prec %wd / %wd\n", prec, exact_cutoff); */
        if (prec >= exact_cutoff)
            gr_ctx_init_fmpq(ctx);
        else
            gr_ctx_init_real_arb(ctx, prec);

        gr_mat_init(RA, A->r, A->c, ctx);
        Rdelta = gr_heap_init(ctx);
        Reta = gr_heap_init(ctx);

        GR_MUST_SUCCEED(gr_mat_set_fmpz_mat(RA, A, ctx));
        GR_MUST_SUCCEED(gr_set_d(Rdelta, delta, ctx));
        GR_MUST_SUCCEED(gr_set_d(Reta, eta, ctx));

        is_reduced = gr_mat_is_row_lll_reduced_naive(RA, Rdelta, Reta, ctx);

        gr_mat_clear(RA, ctx);
        gr_heap_clear(Rdelta, ctx);
        gr_heap_clear(Reta, ctx);

        gr_ctx_clear(ctx);

        if (is_reduced != T_UNKNOWN)
            break;
    }

    return (is_reduced == T_TRUE);
}

