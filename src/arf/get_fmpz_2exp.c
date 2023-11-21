/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

void
arf_get_fmpz_2exp(fmpz_t man, fmpz_t exp, const arf_t x)
{
    if (arf_is_special(x))
    {
        fmpz_zero(man);
        fmpz_zero(exp);
    }
    else
    {
        mp_srcptr xptr;
        mp_size_t xn;
        int shift;

        ARF_GET_MPN_READONLY(xptr, xn, x);

        shift = flint_ctz(xptr[0]);

        fmpz_sub_ui(exp, ARF_EXPREF(x), xn * FLINT_BITS - shift);

        if (xn == 1)
        {
            if (ARF_SGNBIT(x))
                fmpz_neg_ui(man, xptr[0] >> shift);
            else
                fmpz_set_ui(man, xptr[0] >> shift);
        }
        else
        {
            __mpz_struct * z = _fmpz_promote(man);

            if (z->_mp_alloc < xn)
                mpz_realloc(z, xn);

            if (shift == 0)
                flint_mpn_copyi(z->_mp_d, xptr, xn);
            else
                mpn_rshift(z->_mp_d, xptr, xn, shift);

            /* top limb cannot be zero */
            z->_mp_size = ARF_SGNBIT(x) ? -xn : xn;
        }
    }
}

