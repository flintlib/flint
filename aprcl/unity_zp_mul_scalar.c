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
unity_zp_mul_scalar_fmpz(unity_zp_ptr f, unity_zp_srcptr g, const fmpz_t s)
{
    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(g->ctx)));

    fmpz_mod_poly_scalar_mul_fmpz(f->poly, g->poly, s, f->ctx);
}

void
unity_zp_mul_scalar_ui(unity_zp_ptr f, unity_zp_srcptr g, ulong s)
{
    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(g->ctx)));

    fmpz_mod_poly_scalar_mul_ui(f->poly, g->poly, s, f->ctx);
}
