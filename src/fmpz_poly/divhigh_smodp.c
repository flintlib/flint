/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

void fmpz_poly_divhigh_smodp(fmpz * res, const fmpz_poly_t f,
                              const fmpz_poly_t g, const fmpz_t p, slong n)
{
    fmpz_mod_ctx_t ctx;
    fmpz *tA, *tB;
    slong glen, flen, Alen, Blen;

    glen = g->length;
    flen = f->length;
    Blen = FLINT_MIN(glen, n);
    Alen = FLINT_MIN(flen, n);

    if (Alen == 0)
    {
        _fmpz_vec_zero(res, n);
        return;
    }

    fmpz_mod_ctx_init(ctx, p);
    tA = _fmpz_vec_init(Alen);
    tB = _fmpz_vec_init(Blen);

    _fmpz_mod_vec_set_fmpz_vec(tA, f->coeffs + (flen - Alen), Alen, ctx);
    _fmpz_poly_reverse(tA, tA, Alen, Alen);
    _fmpz_mod_vec_set_fmpz_vec(tB, g->coeffs + (glen - Blen), Blen, ctx);
    _fmpz_poly_reverse(tB, tB, Blen, Blen);

    _fmpz_mod_poly_div_series(res, tA, Alen, tB, Blen, n, ctx);
    _fmpz_poly_reverse(res, res, n, n);
    _fmpz_mod_vec_get_fmpz_vec_smod(res, res, n, ctx);

    _fmpz_vec_clear(tA, Alen);
    _fmpz_vec_clear(tB, Blen);
    fmpz_mod_ctx_clear(ctx);
}

