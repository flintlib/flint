/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_rem(fmpz *R,
                        const fmpz *A, slong lenA, const fmpz *B, slong lenB,
                        const fmpz_t invB, const fmpz_mod_ctx_t ctx)
{
    fmpz * Q = _fmpz_vec_init(lenA - lenB + 1);
    _fmpz_mod_poly_divrem(Q, R, A, lenA, B, lenB, invB, ctx);
    _fmpz_vec_clear(Q, lenA - lenB + 1);
}

void fmpz_mod_poly_rem(fmpz_mod_poly_t R, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_t Q;
    fmpz_mod_poly_init(Q, ctx);
    fmpz_mod_poly_divrem(Q, R, A, B, ctx);
    fmpz_mod_poly_clear(Q, ctx);
}

void fmpz_mod_poly_rem_f(fmpz_t f, fmpz_mod_poly_t R, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_t Q;
    fmpz_mod_poly_init(Q, ctx);
    fmpz_mod_poly_divrem_f(f, Q, R, A, B, ctx);
    fmpz_mod_poly_clear(Q, ctx);
}
