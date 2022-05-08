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
fmpz_init_set(fmpz_t f, const fmpz_t g)
{
    if (!COEFF_IS_MPZ(*g))
    {
        *f = *g;
    }
    else
    {
        mpz_mock_ptr mf;
        mpz_mock_ptr mg = COEFF_TO_PTR(*g);
        slong mgabssz = FLINT_ABS(mg->_mp_size);

        mf = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(mf);
        
        if (mf->_mp_alloc < mgabssz)
            _mpz_realloc((mpz_ptr) mf, mgabssz + 1);

        memcpy(mf->_mp_d, mg->_mp_d, sizeof(ulong) * mgabssz);
        mf->_mp_size = mg->_mp_size;
    }
}
