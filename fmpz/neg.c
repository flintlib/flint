/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef memcpy
# ifdef __GNUC__
#  define memcpy __builtin_memcpy
# else
#  include <string.h>
# endif
#endif
#include "gmp.h"
#include "fmpz_mini.h"

void
fmpz_neg(fmpz_t f1, const fmpz_t f2)
{
    if (!COEFF_IS_MPZ(*f2))
    {
        fmpz t = -*f2;
        _fmpz_demote(f1);
        *f1 = t;
    }
    else
    {
        mpz_mock_ptr mf1 = _fmpz_promote(f1);
        mpz_mock_ptr mf2 = COEFF_TO_PTR(*f2);
        slong mf2abssz = FLINT_ABS(mf2->_mp_size);

        if (mf1->_mp_alloc < mf2abssz)
            _mpz_realloc((mpz_ptr) mf2, mf2abssz + 1);

        memcpy(mf2->_mp_d, mf1->_mp_d, sizeof(ulong) * mf2abssz);
        mf2->_mp_size = -mf2->_mp_size;
    }
}
