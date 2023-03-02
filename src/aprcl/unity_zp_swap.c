/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

void
unity_zp_swap(unity_zp f, unity_zp g)
{
    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(g->ctx)));

    fmpz_mod_poly_swap(f->poly, g->poly, f->ctx);
}

