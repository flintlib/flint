/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_DIV_SMALL_H
#define FMPZ_DIV_SMALL_H

#include <gmp.h>
#include "fmpz.h"

/*
    Division of a 1- or 2-limb mpz-sized g by a nonzero divisor of
    magnitude ah (< 2^FLINT_BITS) and sign hneg, without touching the mpz
    layer (one or two hardware divisions). rnd is 0 for truncation, -1 for
    floor and 1 for ceiling; q and r may each be NULL, and may alias g.
    If rem != NULL the nonnegative remainder magnitude is stored there.
    Returns 1 if the case was handled (|g| has at most two limbs), 0
    otherwise.
*/
static inline int
_fmpz_div_qr_small_divisor_ui(fmpz_t q, fmpz_t r, ulong * rem,
    const fmpz_t g, ulong ah, int hneg, int rnd)
{
    mpz_srcptr mg = COEFF_TO_PTR(*g);
    slong gn = mg->_mp_size;
    ulong g0, g1, q0, q1, r0;
    int gneg = gn < 0, qneg = gneg ^ hneg, rneg = gneg;

    if (FLINT_ABS(gn) > 2)
        return 0;

    g0 = mg->_mp_d[0];
    if (gn == 1 || gn == -1)
    {
        q1 = 0;
        q0 = g0 / ah;
        r0 = g0 - q0 * ah;
    }
    else
    {
        g1 = mg->_mp_d[1];
        q1 = g1 / ah;
        r0 = g1 - q1 * ah;
        udiv_qrnnd(q0, r0, r0, g0, ah);
    }

    /* floor rounds down when the truncated remainder has the wrong sign
       (that of g rather than h), ceiling when it has the sign of h */
    if (r0 != 0 && ((rnd < 0 && gneg != hneg) || (rnd > 0 && gneg == hneg)))
    {
        q0++;
        q1 += (q0 == 0);
        r0 = ah - r0;
        rneg = (rnd < 0) ? hneg : !hneg;
    }

    if (q != NULL)
    {
        if (qneg)
            fmpz_neg_uiui(q, q1, q0);
        else
            fmpz_set_uiui(q, q1, q0);
    }

    if (r != NULL)
    {
        if (rneg)
            fmpz_neg_ui(r, r0);
        else
            fmpz_set_ui(r, r0);
    }

    if (rem != NULL)
        *rem = r0;

    return 1;
}

/* the same with a signed small divisor; the magnitude is formed with an
   unsigned negation so that h = WORD_MIN (magnitude 2^(FLINT_BITS-1)) is
   handled without signed overflow (FLINT_ABS(WORD_MIN) is undefined
   behaviour, which optimising compilers exploit) */
static inline int
_fmpz_div_qr_small_divisor(fmpz_t q, fmpz_t r, const fmpz_t g, slong h, int rnd)
{
    ulong ah = (h < 0) ? -(ulong) h : (ulong) h;
    return _fmpz_div_qr_small_divisor_ui(q, r, NULL, g, ah, h < 0, rnd);
}

#endif
