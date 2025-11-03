/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_MAT_IMPL_H
#define GR_MAT_IMPL_H

#include "gr_types.h"

int _gr_mat_write(gr_stream_t out, const gr_mat_t mat, int linebreaks, gr_ctx_t ctx);

#endif
