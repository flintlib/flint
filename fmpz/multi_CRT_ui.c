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


void fmpz_multi_CRT_ui(
    fmpz_t b,
    mp_srcptr in,
    const fmpz_comb_t C,
    fmpz_comb_temp_t CT,
    int sign)
{
    slong i, j, k, l, s;
    slong klen = C->crt_klen;
    slong * step = C->step;
    crt_lut_entry * lu = C->crt_lu;
    fmpz * T = CT->T;
    fmpz * A = CT->A;
    slong * offsets = C->crt_offsets;
    const mp_limb_t * md = C->packed_multipliers;
    mpz_ptr az;
    mp_limb_t * ad;
    mp_limb_t hi, lo, t;

    for (k = 0, i = 0, l = 0; k < klen; k++)
    {
        s = step[k];
        j = offsets[k];
        az = _fmpz_promote(A + k);

        if (s < 0)
        {
            /*
                every low level combination in this chunk has 1 prime
                and md already has lu[i].i0 pre-multiplied in.
            */
            s = -s - 1;

            ad = FLINT_MPZ_REALLOC(az, s + 2);

            flint_mpn_zero(ad, s + 2);
            hi = lo = 0;

            for ( ; i < j; md += s, l++, i++)
            {
                FLINT_ASSERT(lu[i].i0 != 0);
                FLINT_ASSERT(lu[i].i1 == 0);
                FLINT_ASSERT(lu[i].i2 == 0);

                t = mpn_addmul_1(ad, md, s, in[l*1]);
                add_ssaaaa(hi, lo, hi, lo, UWORD(0), t);
            }

            ad[s] = lo;
            ad[s + 1] = hi;
        }
        else
        {
            ad = FLINT_MPZ_REALLOC(az, s + 2);

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
                /*
                    We have lu[i].mod.n = p0*p1*p2, and each lu[i].i{0|1|2} is
                    strictly less than p1*p2*p3, and the inputs are reduced mod pi.
                    Therefore, the sum is at most (p0*p1*p2-1)*(p0-1+p1-1+p2-1).
                    Since p0*p1*p2 fits into a word, the sum fits into two words
                    and the hi word is less than p0*p1*p2.
                */
                }
                else if (lu[i].i1 != 0)
                {
                    FLINT_ASSERT(l + 1 <= C->num_primes);
                    MAC(hi, lo, in[l*1], lu[i].i1); l++;
                    /* Ditto for two */
                }

                FLINT_ASSERT(hi < lu[i].mod.n);
                NMOD_RED2(t, hi, lo, lu[i].mod);

                /* mid level combination: depends on FMPZ_CRT_UI_CUTOFF */
                hi = mpn_addmul_1(ad, md, s, t);
                add_ssaaaa(ad[s + 1], ad[s], ad[s + 1], ad[s], UWORD(0), hi);
            }
        }

        s += 2;

        MPN_NORM(ad, s);
        az->_mp_size = s;
        _fmpz_demote_val(A + k);

        _fmpz_smod(A + k, A + k, C->crt_P->moduli + k, sign, T + 0);
    }

    FLINT_ASSERT(l == C->num_primes);

    /* high level combination */
    fmpz_swap(T + 0, b);
    _fmpz_multi_CRT_precomp(T, C->crt_P, A, sign);
    fmpz_swap(T + 0, b);
}

