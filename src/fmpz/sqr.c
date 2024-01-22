/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"

__GMP_DECLSPEC extern void * (*__gmp_allocate_func)(size_t);
__GMP_DECLSPEC extern void   (*__gmp_free_func)(void *, size_t);

#define GMP_ALLOCATE_LIMBS(nlimbs) __gmp_allocate_func(sizeof(mp_limb_t) * (nlimbs))
#define GMP_FREE_LIMBS(ptr, nlimbs) __gmp_free_func(ptr, sizeof(mp_limb_t) * (nlimbs))

void
fmpz_sqr(fmpz_t res, const fmpz_t xp)
{
    slong xs = *xp;
    mpz_ptr mres;
    mp_ptr resd;
    mp_srcptr xd;
    mp_size_t xsz;
    mp_limb_t upperlimb;

    if (!COEFF_IS_MPZ(xs))
    {
        ulong t[2];
        smul_ppmm(t[1], t[0], xs, xs);
        fmpz_set_uiui(res, t[1], t[0]);
        return;
    }

    {
        mpz_srcptr mxp = COEFF_TO_PTR(xs);
        xd = mxp->_mp_d;
        xsz = FLINT_ABS(mxp->_mp_size);
    }

    mres = _fmpz_promote(res);
    resd = mres->_mp_d;

    /* ADX version work with aliasing up to n = 3. mpn_sqr never works with
     * aliasing. */
#if FLINT_HAVE_ADX
# define N 3
#else
# define N 0
#endif
    if (mres->_mp_alloc < 2 * xsz || (resd == xd && xsz > N))
        resd = GMP_ALLOCATE_LIMBS(2 * xsz + 1); /* Allocate one extra limb */
#undef N

    upperlimb = flint_mpn_sqr(resd, xd, xsz);
    mres->_mp_size = 2 * xsz - (upperlimb == 0);

    if (resd != mres->_mp_d)
    {
        /* FIXME: Is this check needed? */
        if (mres->_mp_alloc)
            GMP_FREE_LIMBS(mres->_mp_d, mres->_mp_alloc);
        mres->_mp_d = resd;
        mres->_mp_alloc = 2 * xsz + 1;
    }
}
