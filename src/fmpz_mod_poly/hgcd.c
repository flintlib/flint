/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "gr_poly.h"

slong _fmpz_mod_poly_hgcd(fmpz **M, slong *lenM,
                     fmpz *A, slong *lenA, fmpz *B, slong *lenB,
                     const fmpz *a, slong lena, const fmpz *b, slong lenb,
                     const fmpz_mod_ctx_t ctx)
{
    slong sgnM;
    gr_ctx_t gr_ctx;

    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_hgcd(NULL, &sgnM, (gr_ptr *) M, lenM, A, lenA, B, lenB, a, lena, b, lenb, FMPZ_MOD_POLY_HGCD_CUTOFF, gr_ctx));

    return sgnM;
}
