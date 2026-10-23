/*
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "factor_impl.h"

int
gr_poly_factor_kaltofen_shoup(gr_ptr c, gr_poly_vec_t fac, fmpz_vec_t exp,
    const gr_poly_t F, gr_ctx_t ctx)
{
    gr_poly_vec_t sqf, dd, ed;
    fmpz_vec_t sqf_exp, degs;
    gr_poly_t frob, frob2;
    gr_poly_preinv_t P;
    fmpz_t q;
    flint_rand_t state;
    slong i, j, k;
    int status = GR_SUCCESS;

    gr_poly_vec_set_length(fac, 0, ctx);
    fmpz_vec_set_length(exp, 0);

    status |= _gr_poly_factor_ff_info(NULL, NULL, NULL, ctx);
    if (status != GR_SUCCESS)
        return status;

    if (F->length == 0)
        return gr_zero(c, ctx);

    fmpz_init(q);
    gr_poly_vec_init(sqf, 0, ctx);
    gr_poly_vec_init(dd, 0, ctx);
    gr_poly_vec_init(ed, 0, ctx);
    fmpz_vec_init(sqf_exp, 0);
    fmpz_vec_init(degs, 0);
    gr_poly_preinv_init(P, ctx);
    gr_poly_init(frob, ctx);
    gr_poly_init(frob2, ctx);
    flint_rand_init(state);

    status |= gr_ctx_fq_order(q, ctx);
    status |= gr_poly_factor_squarefree(c, sqf, sqf_exp, F, ctx);

    for (i = 0; i < sqf->length && status == GR_SUCCESS; i++)
    {
        gr_poly_struct * f = sqf->entries + i;

        if (f->length == 2)
        {
            gr_poly_vec_append_swap(fac, f, ctx);
            fmpz_vec_append(exp, sqf_exp->entries + i);
            continue;
        }

        status |= gr_poly_preinv_set(P, f, ctx);
        status |= gr_poly_preinv_powmod_x_fmpz(frob, q, P, ctx);

        gr_poly_vec_set_length(dd, 0, ctx);
        fmpz_vec_set_length(degs, 0);
        if (f->length > gr_poly_factor_threaded_cutoff && flint_get_num_available_threads() > 1 && gr_ctx_is_threadsafe(ctx) == T_TRUE)
            status |= _gr_poly_factor_distinct_deg_with_frob_threaded(dd, degs, f, P, frob, q, ctx);
        else
            status |= _gr_poly_factor_distinct_deg_with_frob(dd, degs, f, P, frob, q, ctx);

        for (j = 0; j < dd->length && status == GR_SUCCESS; j++)
        {
            slong d = fmpz_get_si(degs->entries + j);

            gr_poly_vec_set_length(ed, 0, ctx);

            if (d > 1 && dd->entries[j].length - 1 > d)
                status |= gr_poly_rem(frob2, frob, dd->entries + j, ctx);

            status |= _gr_poly_factor_equal_deg_with_frob(ed, dd->entries + j, d, frob2, state, ctx);

            for (k = 0; k < ed->length; k++)
            {
                gr_poly_vec_append_swap(fac, ed->entries + k, ctx);
                fmpz_vec_append(exp, sqf_exp->entries + i);
            }
        }
    }

    fmpz_clear(q);
    gr_poly_vec_clear(sqf, ctx);
    gr_poly_vec_clear(dd, ctx);
    gr_poly_vec_clear(ed, ctx);
    fmpz_vec_clear(sqf_exp);
    fmpz_vec_clear(degs);
    gr_poly_preinv_clear(P, ctx);
    gr_poly_clear(frob, ctx);
    gr_poly_clear(frob2, ctx);
    flint_rand_clear(state);

    return status;
}
