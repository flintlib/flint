/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"

/* todo: deal correctly with nested structures */

int gr_ctx_cmp_coercion(gr_ctx_t ctx1, gr_ctx_t ctx2)
{
    if (ctx1->which_ring < ctx2->which_ring)
        return -1;
    if (ctx1->which_ring > ctx2->which_ring)
        return 1;

    if (ctx1->which_ring == GR_CTX_GR_POLY)
    {
        return gr_ctx_cmp_coercion(POLYNOMIAL_ELEM_CTX(ctx1), POLYNOMIAL_ELEM_CTX(ctx2));
    }

    if (ctx1->which_ring == GR_CTX_GR_MAT)
    {
        return gr_ctx_cmp_coercion(MATRIX_CTX(ctx1)->base_ring, MATRIX_CTX(ctx2)->base_ring);
    }

    return 1;
}
