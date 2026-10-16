/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"

void _fmpz_mod_vec_get_fmpz_vec_smod(fmpz * res, const fmpz * vec, slong len, const fmpz_mod_ctx_t ctx)
{
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);
    fmpz_t phalf;
    slong i;

    fmpz_init(phalf);
    fmpz_fdiv_q_2exp(phalf, p, 1);

    for (i = 0; i < len; i++)
    {
        if (fmpz_cmp(vec + i, phalf) > 0)
            fmpz_sub(res + i, vec + i, p);
        else
            fmpz_set(res + i, vec + i);
    }

    fmpz_clear(phalf);
}
