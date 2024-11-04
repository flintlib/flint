/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_POLY_IMPL_H
#define GR_POLY_IMPL_H

#include "gr_types.h"

#ifdef __cplusplus
extern "C" {
#endif

int
_gr_poly_divrem_divconquer_recursive(gr_ptr Q, gr_ptr BQ, gr_ptr W,
    gr_srcptr A, gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif /* GR_POLY_IMPL_H */
