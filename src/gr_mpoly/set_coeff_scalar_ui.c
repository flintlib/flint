/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr_mpoly.h"

int gr_mpoly_set_coeff_scalar_ui(gr_mpoly_t poly,
             gr_srcptr c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    slong i, nvars = mctx->nvars;
    fmpz * newexp;
    int status;
    TMP_INIT;

    TMP_START;
    newexp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (i = 0; i < nvars; i++)
        fmpz_init_set_ui(newexp + i, exp[i]);

    status = gr_mpoly_set_coeff_scalar_fmpz(poly, c, newexp, mctx, cctx);

    for (i = 0; i < nvars; i++)
        fmpz_clear(newexp + i);

    TMP_END;

    return status;
}

int gr_mpoly_set_coeff_ui_ui(
    gr_mpoly_t A,
    ulong c,
    const ulong * exp,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    int status;
    gr_ptr t;

    GR_TMP_INIT(t, cctx);
    status = gr_set_ui(t, c, cctx);
    status |= gr_mpoly_set_coeff_scalar_ui(A, t, exp, mctx, cctx);
    GR_TMP_CLEAR(t, cctx);

    return status;
}

int gr_mpoly_set_coeff_si_ui(
    gr_mpoly_t A,
    slong c,
    const ulong * exp,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    int status;
    gr_ptr t;

    GR_TMP_INIT(t, cctx);
    status = gr_set_si(t, c, cctx);
    status |= gr_mpoly_set_coeff_scalar_ui(A, t, exp, mctx, cctx);
    GR_TMP_CLEAR(t, cctx);

    return status;
}

int gr_mpoly_set_coeff_fmpz_ui(
    gr_mpoly_t A,
    const fmpz_t c,
    const ulong * exp,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    int status;
    gr_ptr t;

    GR_TMP_INIT(t, cctx);
    status = gr_set_fmpz(t, c, cctx);
    status |= gr_mpoly_set_coeff_scalar_ui(A, t, exp, mctx, cctx);
    GR_TMP_CLEAR(t, cctx);

    return status;
}

int gr_mpoly_set_coeff_fmpq_ui(
    gr_mpoly_t A,
    const fmpq_t c,
    const ulong * exp,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    int status;
    gr_ptr t;

    GR_TMP_INIT(t, cctx);
    status = gr_set_fmpq(t, c, cctx);
    if (status == GR_SUCCESS)
        status |= gr_mpoly_set_coeff_scalar_ui(A, t, exp, mctx, cctx);
    GR_TMP_CLEAR(t, cctx);

    return status;
}
