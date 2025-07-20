/*
   Copyright (C) 2025 Marc Mezzarobba

   This file is part of FLINT.

   FLINT is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License (LGPL) as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_vec.h"
#include "perm.h"

static int
cmp(const void * i, const void * j, void * data)
{
    gr_vec_struct * shifts = data;
    gr_ctx_t ZZ;
    gr_ctx_init_fmpz(ZZ);
    int res = fmpz_cmp(gr_vec_entry_ptr(shifts, *(slong *) i, ZZ),
                       gr_vec_entry_ptr(shifts, *(slong *) j, ZZ));
    gr_ctx_clear(ZZ);
    return res;
}

int
_gr_poly_shiftless_decomposition_from_factors(
        gr_vec_t slfac, gr_vec_t slshifts, gr_vec_t slmult,
        const gr_vec_t fac, const gr_vec_t mult,
        gr_ctx_t ctx)
{
    gr_ctx_t pctx, ZZ, ZZvec;
    fmpz_t shift;
    gr_ptr delta;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_vector_gr_vec(ZZvec, ZZ);
    gr_ctx_init_gr_poly(pctx, ctx);
    fmpz_init(shift);
    GR_TMP_INIT(delta, ctx);

    int status = GR_SUCCESS;

    slfac->length = 0;
    slshifts->length = 0;
    slmult->length = 0;

    /* Detect shift-equivalent factors */

    for (slong i = 0; i < fac->length; i++)
    {
        const gr_poly_struct * f = gr_vec_entry_srcptr(fac, i, pctx);
        const fmpz * m = gr_vec_entry_srcptr(mult, i, ZZ);
        slong j = 0;

        for (; j < slfac->length; j++)
        {
            const gr_poly_struct * sl = gr_vec_entry_srcptr(slfac, j, pctx);

            int equiv = gr_poly_shift_equivalent(shift, sl, f, ctx);
            if (equiv == T_TRUE)
                break;
            if (equiv == T_UNKNOWN)
            {
                status = GR_UNABLE;
                goto cleanup;
            }
        }

        if (j == slfac->length)
        {
            status |= gr_vec_append(slfac, f, pctx);
            gr_vec_fit_length(slshifts, ++slshifts->length, ZZvec);
            gr_vec_fit_length(slmult, ++slmult->length, ZZvec);
        }

        gr_vec_append_swap(gr_vec_entry_ptr(slshifts, j, ZZvec), shift, ZZ);
        status |= gr_vec_append(gr_vec_entry_ptr(slmult, j, ZZvec), m, ZZ);
    }

    /* Sort the shifts of each factor and normalize the smallest to zero */

    for (slong j = 0; j < slfac->length; j++)
    {

        gr_poly_struct * f = gr_vec_entry_ptr(slfac, j, pctx);
        gr_vec_struct * shifts = gr_vec_entry_ptr(slshifts, j, ZZvec);
        gr_vec_struct * shmult = gr_vec_entry_ptr(slmult, j, ZZvec);

        slong * perm = _perm_init(shifts->length);

        flint_sort(perm, shifts->length, sizeof(slong), cmp, shifts);
        _perm_inv(perm, perm, shifts->length);

        status |= gr_vec_permute(shifts, shifts, perm, ZZ);
        status |= gr_vec_permute(shmult, shmult, perm, ZZ);

        status |= gr_set_fmpz(delta, shifts->entries, ctx);
        status |= gr_poly_taylor_shift(f, f, delta, ctx);
        status |= _gr_vec_sub_scalar((fmpz *) shifts->entries + 1,
                                     (fmpz *) shifts->entries + 1,
                                     shifts->length - 1,
                                     shifts->entries, ZZ);
        fmpz_zero(shifts->entries);

        _perm_clear(perm);
    }

cleanup:
    GR_TMP_CLEAR(delta, ctx);
    fmpz_clear(shift);
    gr_ctx_clear(pctx);
    gr_ctx_clear(ZZvec);
    gr_ctx_clear(ZZ);

    return status;
}

int
gr_poly_shiftless_decomposition_from_factors(
        gr_vec_t slfac, gr_vec_t slshifts, gr_vec_t slmult,
        const gr_vec_t fac, const gr_vec_t mult,
        gr_ctx_t ctx)
{
    gr_ctx_t pctx;
    gr_vec_t _slfac;

    gr_ctx_init_gr_poly(pctx, ctx);
    gr_vec_init(_slfac, 0, pctx);

    int status = _gr_poly_shiftless_decomposition_from_factors(
            _slfac, slshifts, slmult, fac, mult, ctx);
    FLINT_SWAP(gr_vec_struct, *slfac, *_slfac);

    gr_vec_clear(_slfac, pctx);
    gr_ctx_clear(pctx);

    return status;
}

int
gr_poly_shiftless_decomposition_factor(
        gr_ptr c, gr_vec_t slfac, gr_vec_t slshifts, gr_vec_t slmult,
        const gr_poly_t pol,
        gr_ctx_t ctx)
{
    gr_ctx_t pctx, ZZ;
    gr_vec_t mult;
    gr_poly_t polc;

    if (gr_poly_is_zero(pol, ctx) == T_TRUE)
        return GR_DOMAIN;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(pctx, ctx);

    gr_poly_init(polc, ctx);
    gr_vec_init(mult, 0, ZZ);

    int status = GR_SUCCESS;

    status |= gr_factor(polc, slfac, mult, pol, 0, pctx);
    if (status != GR_SUCCESS)
        goto cleanup;

    gr_poly_fit_length(polc, 1, ctx);
    gr_swap(c, polc->coeffs, ctx);
    _gr_poly_normalise(polc, ctx);
    status |= gr_poly_shiftless_decomposition_from_factors(
            slfac, slshifts, slmult, slfac, mult, ctx);

cleanup:
    gr_poly_clear(polc, ctx);
    gr_vec_clear(mult, ZZ);

    gr_ctx_clear(pctx);
    gr_ctx_clear(ZZ);

    return status;
}

int
gr_poly_shiftless_decomposition(
        gr_ptr c, gr_vec_t slfac, gr_vec_t slshifts, gr_vec_t slmult,
        const gr_poly_t pol,
        gr_ctx_t ctx)
{
    return gr_poly_shiftless_decomposition_factor(c, slfac, slshifts, slmult,
                                                  pol, ctx);
}
