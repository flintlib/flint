/*
    Copyright (C) 2009 William Hart

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

/* defined in addmul.c */
void
_flint_mpz_addmul_large(mpz_ptr z, mpz_srcptr x, mpz_srcptr y, int negate);

void fmpz_submul(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1, c2, c3;
    __mpz_struct * mf;

    c1 = *g;
	c2 = *h;
    c3 = *f;

    /* todo: are the zero checks worth it for small input? */

    if (c1 == 0 || c2 == 0)
        return;

    if (c3 == 0)
    {
        fmpz_mul(f, g, h);
        fmpz_neg(f, f);
        return;
    }

    if (!COEFF_IS_MPZ(c1))  /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* both are small */
        {
            ulong p1, p0;
            smul_ppmm(p1, p0, -c1, c2);

            if (!COEFF_IS_MPZ(c3))
            {
                ulong F1 = FLINT_SIGN_EXT(c3);
                add_ssaaaa(p1, p0, p1, p0, F1, c3);

                fmpz_set_signed_uiui(f, p1, p0);
            }
            else
            {
                mpz_ptr pF = COEFF_TO_PTR(c3);
                flint_mpz_add_signed_uiui(pF, pF, p1, p0);
                _fmpz_demote_val(f);  /* cancellation may have occurred	*/
            }
        }
        else
        {
            fmpz_addmul_si(f, h, -c1);
        }
    }
    else if (!COEFF_IS_MPZ(c2))  /* h is small */
	{
		fmpz_addmul_si(f, g, -c2);
	}
    else
    {
        mf = _fmpz_promote_val(f);
        _flint_mpz_addmul_large(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2), 1);
        _fmpz_demote_val(f);  /* cancellation may have occurred	*/
    }
}
