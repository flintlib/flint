/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void radix_mulmid_naive(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong zlo, slong zhi, const radix_t radix)
{
    fmpz *Z, *A, *B;
    fmpz_t cy, n;
    slong i;
    nmod_t mod = radix->B;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(an >= bn);

    an = FLINT_MIN(an, zhi);
    bn = FLINT_MIN(bn, zhi);

    FLINT_ASSERT(zhi <= an + bn);
    FLINT_ASSERT(zhi > zlo);
    FLINT_ASSERT(zlo >= 0);

    A = _fmpz_vec_init(an);
    B = _fmpz_vec_init(bn);
    Z = _fmpz_vec_init(an + bn);
    fmpz_init(cy);
    fmpz_init(n);

    for (i = 0; i < an; i++)
        fmpz_set_ui(A + i, a[i]);
    for (i = 0; i < bn; i++)
        fmpz_set_ui(B + i, b[i]);

    if (an >= bn)
        _fmpz_poly_mul(Z, A, an, B, bn);
    else
        _fmpz_poly_mul(Z, B, bn, A, an);

    fmpz_set_ui(n, mod.n);
    for (i = zlo; i < zhi; i++)
    {
        fmpz_add(Z + i, Z + i, cy);
        fmpz_tdiv_qr(cy, Z + i, Z + i, n);
        z[i - zlo] = fmpz_get_ui(Z + i);
    }

    _fmpz_vec_clear(A, an);
    _fmpz_vec_clear(B, bn);
    _fmpz_vec_clear(Z, an + bn);
    fmpz_clear(cy);
    fmpz_clear(n);
}

