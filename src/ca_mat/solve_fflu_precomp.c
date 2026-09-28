/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"
#include "gr.h"
#include "gr_mat.h"

void
ca_mat_solve_fflu_precomp(ca_mat_t X, const slong * perm,
    const ca_mat_t A, const ca_t den, const ca_mat_t B, ca_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);

    /* The input is assumed (not checked) to be a valid factorization;
       like the former implementation, we do not abort on failed divisions. */
    GR_IGNORE(gr_mat_nonsingular_solve_fflu_precomp((gr_mat_struct *) X, perm,
        (const gr_mat_struct *) A, (const gr_mat_struct *) B, gr_ctx));

    ca_mat_div_ca(X, X, den, ctx);
}
