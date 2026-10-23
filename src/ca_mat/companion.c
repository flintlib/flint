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

int
ca_mat_companion(ca_mat_t A, const ca_poly_t poly, ca_ctx_t ctx)
{
    slong n = ca_mat_nrows(A);
    gr_ctx_t gr_ctx;

    if (n != poly->length - 1 || n != ca_mat_ncols(A))
        return 0;

    /* Special values are not valid elements of the gr field
       (ca_inv would map infinities to zero). */
    if (CA_IS_SPECIAL(poly->coeffs + n))
        return 0;

    /* The gr version does not invert the leading coefficient when n = 0.
       For a non-special x, 1/x is special iff x is not provably nonzero. */
    if (n == 0)
        return ca_check_is_zero(poly->coeffs, ctx) == T_FALSE;

    _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);
    return gr_mat_companion((gr_mat_struct *) A,
        (const gr_poly_struct *) poly, gr_ctx) == GR_SUCCESS;
}
