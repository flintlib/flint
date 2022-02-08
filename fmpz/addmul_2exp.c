/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <gmp.h>
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_addmul_2exp(fmpz_t f, const fmpz_t g, ulong exp)
{
    slong f1 = *f;
    slong g1 = *g;
    mp_size_t limbshift;
    __mpz_struct * mf;

    /* if (g1 == 0) */
    /*     return; */
    /* else if (exp == 0) */
    /* { */
    /*     fmpz_add(f, f, g); */
    /*     return; */
    /* } */

    limbshift = exp / FLINT_BITS;
    exp %= FLINT_BITS;

    if (!COEFF_IS_MPZ(f1))
    {
        ulong absx = FLINT_ABS(f1);
        if (!COEFF_IS_MPZ(g1))
        {
            ulong x1, x2;
            x2 = FLINT_ABS(g1);
            x1 = x2 << exp;
            x2 >>= FLINT_BITS - exp;

            if ((f1 ^ g1) >= 0) /* Same sign */
            {
                if (limbshift == 0)
                {
                    add_ssaaaa(x2, x1, x2, x1, 0, absx);
                    if (x2 == 0 && x1 <= COEFF_MAX)
                    {
                        *f = (f1 >= 0) ? x1 : -x1;
                    }
                    else
                    {
smallmpz:
                        mf = _fmpz_new_mpz();
                        *f = PTR_TO_COEFF(mf);
                        mf->_mp_d[0] = x1;
                        mf->_mp_d[1] = x2;
                        mf->_mp_size = 1 + (x2 != 0);
                        if (g1 < 0)
                            mf->_mp_size = -mf->_mp_size;
                    }
                }
                else
                {
                    mf = _fmpz_new_mpz();
                    *f = PTR_TO_COEFF(mf);
                    if (mf->_mp_alloc < (limbshift + 2))
                    {
                        mf->_mp_d = flint_realloc(mf->_mp_d,
                                        sizeof(mp_limb_t) * (limbshift + 2));
                        mf->_mp_alloc = limbshift + 2;
                    }
                    memset(mf->_mp_d, 0, sizeof(mp_limb_t) * limbshift);
                    mf->_mp_d[0] = absx;
                    mf->_mp_d[limbshift] = x1;
                    mf->_mp_d[limbshift + 1] = x2;
                    mf->_mp_size = limbshift + 1 + (x2 != 0);
                    if (f1 < 0)
                        mf->_mp_size = -mf->_mp_size;
                }
            }
            else /* Different signs */
            {
                if (limbshift == 0)
                {
                    if (x2 == 0)
                    {
                        if (x1 >= absx)
                            x1 -= absx;
                        else
                        {
                            g1 = f1;
                            x1 = absx - x1;
                        }
                    }
                    else
                        sub_ddmmss(x2, x1, x2, x1, 0, absx);

                    if (x2 == 0 && x1 <= COEFF_MAX)
                        *f = (g1 >= 0) ? x1 : -x1;
                    else
                        goto smallmpz;
                }
                else
                {
                    mf = _fmpz_new_mpz();
                    sub_ddmmss(x2, x1, x2, x1, 0, 1);
                    *f = PTR_TO_COEFF(mf);
                    if (mf->_mp_alloc < (limbshift + 2))
                    {
                        mf->_mp_d = flint_realloc(mf->_mp_d,
                                        sizeof(mp_limb_t) * (limbshift + 2));
                        mf->_mp_alloc = limbshift + 2;
                    }
                    memset(mf->_mp_d, UINT_MAX, sizeof(mp_limb_t) * limbshift);
                    mf->_mp_d[0] = g1;
                    mf->_mp_d[limbshift] = x1;
                    mf->_mp_d[limbshift + 1] = x2;
                    mf->_mp_size = limbshift + (x2 != 0) + (x2 != 0 && x1 != 0);
                    if (g1 < 0)
                        mf->_mp_size = -mf->_mp_size;
                    else
                        mf->_mp_d[0] = -g1;
                }
            }
        }
        else /* f is small, g is large */
        {
            __mpz_struct * mres = _fmpz_new_mpz(); /* mres represents f */
            mp_limb_t * mreslimbs, * mflimbs;
            mp_size_t abssz;
            *f = PTR_TO_COEFF(mres);
            mf = COEFF_TO_PTR(g1); /* mf represents g */
            abssz = FLINT_ABS(mf->_mp_size);
            mflimbs = mf->_mp_d;
            if (mres->_mp_alloc < (limbshift + abssz + 1))
            {
                mres->_mp_d = flint_realloc(mres->_mp_d,
                                sizeof(mp_limb_t) * (limbshift + abssz + 1));
                mres->_mp_alloc = limbshift + abssz + 1;
            }
            mreslimbs = mf->_mp_d + limbshift;

            /* Shift g */
            g1 = mpn_lshift(mreslimbs, mflimbs, abssz, exp);
            mreslimbs[abssz] = g1;

            /* And add original f to the shifted g */
            mreslimbs -= limbshift;
            abssz += limbshift + 1;
            if ((mf->_mp_size > 0) != (f1 > 0)) /* Same sign */
                mpn_add_1(mreslimbs, mreslimbs, abssz, absx); /* No carry */
            else
                mpn_sub_1(mreslimbs, mreslimbs, abssz, absx); /* No borrow */

            /* Set size */
            mres->_mp_size
                = abssz - (mreslimbs[abssz - 1] == 0)
                - (mreslimbs[abssz - 1] == 0 && mreslimbs[abssz - 2] == 0);
            if (mf->_mp_size < 0)
                mres->_mp_size = -mres->_mp_size;
        }
    }
    else
    {
        if (!COEFF_IS_MPZ(g1)) /* f is large, g is small */
        {
        }
        else /* both large */
        {
        }
    }
}
