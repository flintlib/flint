/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"

void fmpz_sqrt(fmpz_t f, const fmpz_t g)
{
    if (fmpz_sgn(g) < 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_sqrt). g is negative.\n");
    }

    if (!COEFF_IS_MPZ(*g))
        fmpz_set_ui(f, n_sqrt(*g));
    else
    {
        mpz_srcptr mg = COEFF_TO_PTR(*g);
        mpz_ptr mf;

        if (mg->_mp_size <= 2)
        {
            /* one- or two-limb input: the root fits in a limb */
            mp_limb_t sd;
            _flint_mpn_sqrtrem(&sd, NULL, mg->_mp_d, mg->_mp_size);
            fmpz_set_ui(f, sd);
            return;
        }

        mf = _fmpz_promote(f);
        flint_mpz_sqrt(mf, mg);
        _fmpz_demote_val(f);
    }
}
