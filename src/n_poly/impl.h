/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef N_POLY_IMPL_H
#define N_POLY_IMPL_H

#include "n_poly_types.h"
#include "fq_nmod_types.h"

void n_fq_poly_mullow_(n_fq_poly_t A, const n_fq_poly_t B, const n_fq_poly_t C, slong order, const fq_nmod_ctx_t ctx, n_poly_stack_t St);

#endif
