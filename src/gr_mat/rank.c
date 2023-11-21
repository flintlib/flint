/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "gr_mat.h"

/* todo: different algorithms */
int
gr_mat_rank(slong * rank, const gr_mat_t A, gr_ctx_t ctx)
{
    truth_t dom;

    dom = gr_ctx_is_integral_domain(ctx);

    if (dom == T_TRUE)
        return gr_mat_rank_fflu(rank, A, ctx);

    if (dom == T_FALSE)
        return GR_DOMAIN;

    return GR_UNABLE;
}
