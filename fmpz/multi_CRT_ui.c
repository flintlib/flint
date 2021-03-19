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
#include "mpn_extras.h"
#include "fmpz.h"
#include "nmod_vec.h"

#define MAC(h, l, a, b)                 \
do {                                    \
    mp_limb_t p1, p0;                   \
    umul_ppmm(p1, p0, a, b);            \
    add_ssaaaa(h, l, h, l, p1, p0);     \
} while (0)

/*
    TODO: a slight speedup (< 10%) can be achived in the case of few primes
          where the switching in the low level combination is a bottleneck.
*/
void fmpz_multi_CRT_ui(
    fmpz_t b,
    mp_srcptr in,
    const fmpz_comb_t C,
    fmpz_comb_temp_t CT,
    int sign)
{
    slong i, k, l;
    slong klen = C->crt_klen;
    slong * step = C->step;
    crt_lut_entry * lu = C->crt_lu;
    fmpz * T = CT->T;
    fmpz * A = CT->A;
    slong * offsets = C->crt_offsets;
    const mp_limb_t * md = C->packed_multipliers;

    for (k = 0, i = 0, l = 0; k < klen; k++)
    {
        slong s = step[k];
        slong j = offsets[k];
        mpz_ptr az = _fmpz_promote(A + k);
        mp_limb_t * ad = FLINT_MPZ_REALLOC(az, s + 2);
        mp_limb_t t, hi, lo;

        flint_mpn_zero(ad, s + 2);

        for ( ; i < j; md += s, i++)
        {
            /* low level combination: 1, 2, or 3 small primes */
            FLINT_ASSERT(l + 1 <= C->num_primes);
            umul_ppmm(hi, lo, in[l*1], lu[i].i0); l++;

            if (lu[i].i2 != 0)
            {
                FLINT_ASSERT(l + 2 <= C->num_primes);
                MAC(hi, lo, in[l*1], lu[i].i1); l++;
                MAC(hi, lo, in[l*1], lu[i].i2); l++;
            }
            else if (lu[i].i1 != 0)
            {
                FLINT_ASSERT(l + 1 <= C->num_primes);
                MAC(hi, lo, in[l*1], lu[i].i1); l++;
            }

            FLINT_ASSERT(hi < lu[i].mod.n);
            NMOD_RED2(t, hi, lo, lu[i].mod);

            /* mid level combination: depends on FMPZ_CRT_UI_CUTOFF */
            hi = mpn_addmul_1(ad, md, s, t);
            add_ssaaaa(ad[s + 1], ad[s], ad[s + 1], ad[s], UWORD(0), hi);
        }

        s += 2;
        MPN_NORM(ad, s);
        az->_mp_size = s;
        _fmpz_demote_val(A + k);

        _fmpz_mods(A + k, A + k, C->crt_P->moduli + k, sign, T + 0);
    }

    FLINT_ASSERT(l == C->num_primes);

    /* high level combination */
    fmpz_swap(T + 0, b);
    _fmpz_multi_CRT_run(T, C->crt_P, A, sign);
    fmpz_swap(T + 0, b);
}

