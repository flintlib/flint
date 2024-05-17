/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "longlong.h"
#include "fmpz.h"

void
fmpz_mul2_uiui(fmpz_t f, const fmpz_t g, ulong h0, ulong h1)
{
    fmpz gs = *g;

    if (h0 == 0 || h1 == 0 || gs == 0)
    {
        fmpz_zero(f);
        return;
    }

    /* Now |f| >= |g| */

    /* NOTE: Order important for x86 */
    umul_ppmm(h1, h0, h0, h1);

    /* Treat g first to reduce latency penalty of checking h1 */
    if (!COEFF_IS_MPZ(gs))
    {
        int sign = gs < 0;
        gs = FLINT_ABS(gs);

        if (h1 == UWORD(0))
        {
            mpz_ptr mf;

            umul_ppmm(h1, h0, h0, gs);

            if (!COEFF_IS_MPZ(*f))
            {
                if (h1 == UWORD(0) && h0 <= COEFF_MAX)
                {
                    *f = (sign) ? -h0 : h0;
                    return;
                }
                else
                {
                    mf = _fmpz_new_mpz();
semi_small:         mf->_mp_d[0] = h0;
                    mf->_mp_d[1] = h1;
                    mf->_mp_size = 1 + (h1 != 0);
                    if (sign)
                        mf->_mp_size = -mf->_mp_size;

                    *f = PTR_TO_COEFF(mf);
                    return;
                }
            }
            else
            {
                if (h1 == UWORD(0) && h0 <= COEFF_MAX)
                {
                    _fmpz_clear_mpz(*f);
                    *f = (sign) ? -h0 : h0;
                    return;
                }
                else
                {
                    mf = COEFF_TO_PTR(*f);
                    goto semi_small;
                }
            }
        }
        else
        {
            /* h1 != 0, but gs small */
            ulong f0, f1, f2;
            mpz_ptr mf;

            umul_ppmm(f1, f0, gs, h0);
            umul_ppmm(f2, gs, gs, h1);
            mf = _fmpz_promote(f);
            add_ssaaaa(f2, f1, f2, f1, 0, gs);

            if (f2 == UWORD(0))
            {
                mf->_mp_d[0] = f0;
                mf->_mp_d[1] = f1;
                mf->_mp_size = (sign) ? -2 : 2;
            }
            else
            {
                mp_ptr fp = FLINT_MPZ_REALLOC(mf, 3);
                fp[0] = f0;
                fp[1] = f1;
                fp[2] = f2;
                mf->_mp_size = (sign) ? -3 : 3;
            }
        }
    }
    else
    {
        mpz_srcptr mg = COEFF_TO_PTR(gs);
        mpz_ptr mf = _fmpz_promote(f);
        mpz_t mu;
        mp_limb_t u[2] = {h0, h1};

        mu->_mp_d = u;
        mu->_mp_size = 1 + (h1 != UWORD(0));

        mpz_mul(mf, mg, mu);
    }
}
