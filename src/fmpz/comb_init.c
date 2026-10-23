/*
    Copyright (C) 2008, 2009, William Hart
    Copyright (C) 2010, 2026 Fredrik Johansson
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"

void fmpz_comb_init2(fmpz_comb_t C, nn_srcptr m, slong len, int flags)
{
    if (len < 1)
        flint_throw(FLINT_ERROR, "fmpz_comb_init: len should be positive");

    C->crt = flint_malloc(sizeof(flint_mpn_crt_struct));
    flint_mpn_crt_init2(C->crt, m, len, flags);
    C->num_primes = len;
    C->primes = C->crt->primes;
}

void fmpz_comb_init(fmpz_comb_t C, nn_srcptr m, slong len)
{
    fmpz_comb_init2(C, m, len, FMPZ_COMB_MOD | FMPZ_COMB_CRT);
}

void fmpz_comb_temp_init(fmpz_comb_temp_t CT, const fmpz_comb_t C)
{
    CT->tmp = flint_malloc(sizeof(ulong) * C->crt->tmp_limbs);
    CT->out = (C->crt->flags & FLINT_MPN_CRT_CRT) ?
                flint_malloc(sizeof(ulong) * C->crt->prod_len) : NULL;
}
