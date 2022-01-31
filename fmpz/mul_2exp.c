/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include <string.h>

void
fmpz_mul_2exp(fmpz_t f, const fmpz_t g, ulong exp)
{
    slong c1 = *g;
    ulong c1abs, c1bits;
    __mpz_struct * mf;

    if (c1 == 0)
    {
        fmpz_zero(f);
        return;
    }

    c1abs = FLINT_ABS(c1);
    c1bits = FLINT_BIT_COUNT(c1abs);

    if (c1bits + exp <= SMALL_FMPZ_BITCOUNT_MAX)  /* Result fits inside a small fmpz */
    {
        if (COEFF_IS_MPZ(*f))
            _fmpz_clear_mpz(*f);

        *f = c1 << exp;
    }
    else if (c1bits <= SMALL_FMPZ_BITCOUNT_MAX)   /* g is small */
    {
        ulong expred = exp % GMP_LIMB_BITS;
        int alloc = 1 + (FLINT_BITS != GMP_LIMB_BITS) + exp / GMP_LIMB_BITS
            + ((c1bits + expred) > GMP_LIMB_BITS);
        mp_limb_t * limbs;
        slong tmp;

        /* Ensure enough limbs are allocated for f */
        if (!COEFF_IS_MPZ(*f))
        {
            /* TODO: Initialize the new mpz with alloc limbs instead of
             * reallocating them. */
            mf = _fmpz_new_mpz();
            *f = PTR_TO_COEFF(mf);
            _mpz_realloc(mf, alloc);
        }
        else
        {
            mf = COEFF_TO_PTR(*f);
            if (mf->_mp_alloc < alloc)
                _mpz_realloc(mf, alloc);
        }
        limbs = mf->_mp_d;
        mf->_mp_size = (c1 > 0) ? alloc : -alloc;
        memset(limbs, 0, sizeof(mp_limb_t) * alloc);

        tmp = (slong) (c1bits + expred) - GMP_LIMB_BITS;
        if (tmp <= 0)
        {
            limbs[alloc - 1] = c1abs << expred;
#if FLINT_BITS != GMP_LIMB_BITS
            return;
#endif
        }
#if FLINT_BITS == GMP_LIMB_BITS
        else
        {
            limbs[alloc - 1] = c1abs >> (FLINT_BITS - expred);
            limbs[alloc - 2] = c1abs << expred;
        }
#else   /* FLINT_BITS == 2 * GMP_LIMB_BITS */
        tmp -= GMP_LIMB_BITS;
        if (tmp <= 0)
        {
            ulong * haha = (ulong *) (limbs + alloc - 2);
            *haha = c1abs << expred;
        }
        else
        {
            ulong * haha = (ulong *) (limbs + alloc - 3);
            *haha = c1abs << expred;
            limbs[alloc - 1] = c1abs >> tmp;
        }
#endif
    }
    else                                /* g is large */
    {
        __mpz_struct * mg = COEFF_TO_PTR(c1);

        if (!COEFF_IS_MPZ(*f))
        {
            /* TODO: Initialize the new mpz with alloc limbs instead of
             * reallocating them. */
            mf = _fmpz_new_mpz();
            *f = PTR_TO_COEFF(mf);
            _mpz_realloc(mf, FLINT_ABS(mg->_mp_size) + exp / GMP_LIMB_BITS + 1);
        }
        else
            mf = COEFF_TO_PTR(*f);

        mpz_mul_2exp(mf, mg, exp);
    }
}
