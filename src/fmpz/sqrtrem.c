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

void fmpz_sqrtrem(fmpz_t f, fmpz_t r, const fmpz_t g)
{
    if (fmpz_sgn(g) < 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_sqrtrem). g is negative.\n");
    }

    if (!COEFF_IS_MPZ(*g))
    {
        if (COEFF_IS_MPZ(*r))
            _fmpz_clear_mpz(*r);
        fmpz_set_ui(f, n_sqrtrem((ulong *) r, *g));
    }
    else
    {
        mpz_srcptr g_mpz_ptr = COEFF_TO_PTR(*g);
        mpz_ptr r_mpz_ptr, f_mpz_ptr;

        if (g_mpz_ptr->_mp_size <= 2)
        {
            /* one- or two-limb input: the root fits in a limb and the
               remainder in two, so no mpz outputs are needed */
            mp_limb_t sd, rd[2];
            _flint_mpn_sqrtrem(&sd, rd, g_mpz_ptr->_mp_d, g_mpz_ptr->_mp_size);
            fmpz_set_ui(f, sd);
            fmpz_set_uiui(r, rd[1], rd[0]);
            return;
        }

        _fmpz_promote(f); /* must not hang on to pointer whilst promoting */
        r_mpz_ptr = _fmpz_promote(r);
        f_mpz_ptr = COEFF_TO_PTR(*f);

        flint_mpz_sqrtrem(f_mpz_ptr, r_mpz_ptr, g_mpz_ptr);
        _fmpz_demote_val(f);
        _fmpz_demote_val(r);
    }
}
