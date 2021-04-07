/*
    Copyright (C) 2008, 2009, William Hart 
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "nmod_vec.h"


void fmpz_multi_mod_ui(
    mp_limb_t * out,
    const fmpz_t input,
    const fmpz_comb_t C,
    fmpz_comb_temp_t CT)
{
    slong i, k, l;
    slong stride = 1;
    fmpz * A = CT->A;
    mod_lut_entry * lu;
    slong * offsets;
    slong klen = C->mod_klen;
    fmpz_t ttt;

    /* high level split */
    if (klen == 1)
    {
        *ttt = A[0];
        A[0] = *input;
    }
    else
    {
        _fmpz_multi_mod_precomp(A, C->mod_P, input, -1, CT->T);
    }

    offsets = C->mod_offsets;
    lu = C->mod_lu;

    for (k = 0, i = 0, l = 0; k < klen; k++)
    {
        slong j = offsets[k];

        for ( ; i < j; i++)
        {
            /* mid level split: depends on FMPZ_MOD_UI_CUTOFF */
            mp_limb_t t = fmpz_get_nmod(A + k, lu[i].mod);

            /* low level split: 1, 2, or 3 small primes */
            if (lu[i].mod2.n != 0)
            {
                FLINT_ASSERT(l + 3 <= C->num_primes);
                NMOD_RED(out[l*stride], t, lu[i].mod0); l++;
                NMOD_RED(out[l*stride], t, lu[i].mod1); l++;
                NMOD_RED(out[l*stride], t, lu[i].mod2); l++;
            }
            else if (lu[i].mod1.n != 0)
            {
                FLINT_ASSERT(l + 2 <= C->num_primes);
                NMOD_RED(out[l*stride], t, lu[i].mod0); l++;
                NMOD_RED(out[l*stride], t, lu[i].mod1); l++;
            }
            else
            {
                FLINT_ASSERT(l + 1 <= C->num_primes);
                out[l*stride] = t; l++;
            }
        }
    }

    FLINT_ASSERT(l == C->num_primes);

    if (klen == 1)
        A[0] = *ttt;
}

