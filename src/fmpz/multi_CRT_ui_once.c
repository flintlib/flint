/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz/impl.h"

void
fmpz_multi_CRT_ui_once(fmpz_t output, fmpz_t prod, nn_srcptr residues,
                       nn_srcptr primes, slong num_primes, int sign)
{
    nn_ptr out;
    slong outn;
    int negative;

    if (num_primes == 0)
    {
        fmpz_zero(output);
        if (prod != NULL)
            fmpz_one(prod);
        return;
    }

    out = flint_malloc(num_primes * sizeof(ulong) * (prod != NULL ? 2 : 1));
    negative = flint_mpn_multi_crt_once(out, &outn, (prod != NULL) ? out + num_primes : NULL,
                                        residues, primes, num_primes, sign);
    _fmpz_set_ui_array_neg(output, out, outn, negative);
    if (prod != NULL)
        fmpz_set_ui_array(prod, out + num_primes, outn);
    flint_free(out);
}

void
fmpz_multi_mod_ui_once(ulong * out, const fmpz_t in, nn_srcptr primes, slong num_primes)
{
    fmpz_comb_t comb;
    fmpz_comb_temp_t temp;

    fmpz_comb_init2(comb, primes, num_primes, FMPZ_COMB_MOD);
    fmpz_comb_temp_init(temp, comb);
    fmpz_multi_mod_ui(out, in, comb, temp);
    fmpz_comb_temp_clear(temp);
    fmpz_comb_clear(comb);
}
