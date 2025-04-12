/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_ore_poly.h"

void
gr_ctx_init_gr_ore_poly(gr_ctx_t ctx, gr_ctx_t base_ring, slong base_var, const ore_poly_which_operator_t operator)
{
    gr_ore_poly_ctx_init(ctx, base_ring, base_var, operator);
}
