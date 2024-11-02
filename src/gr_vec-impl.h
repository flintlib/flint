/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_VEC_IMPL_H
#define GR_VEC_IMPL_H

#include "gr_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

int _gr_vec_parallel_reduce(gr_ptr res, gr_method_vec_reduce_op basecase, gr_srcptr vec, slong n, gr_ctx_t ctx, int thread_limit, int flags);

#ifdef __cplusplus
}
#endif

#endif /* GR_VEC_IMPL_H */
