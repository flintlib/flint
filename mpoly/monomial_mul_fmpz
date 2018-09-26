/*
    Copyright (C) 2018 Daniel Schultz
    This file is part of FLINT.
    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/
#include "mpoly.h"


void mpoly_monomial_mul_fmpz(ulong * exp2, const ulong * exp3,
                                                            slong N, fmpz_t c)
{
    FLINT_ASSERT(fmpz_sgn(c) >= 0);

    if (!COEFF_IS_MPZ(*c))
    {
        mpn_mul_1(exp2, exp3, N, *c);
    }
    else
    {
        slong i, cn;
        ulong * cp;

        cn = COEFF_TO_PTR(*c)->_mp_size;
        cp = COEFF_TO_PTR(*c)->_mp_d;
        FLINT_ASSERT(0 < cn);
        FLINT_ASSERT(cn <= N);

        if (exp2 != exp3)
        {
            mpn_mul_1(exp2, exp3, N, cp[0]);
            for (i = 1; i < cn; i++)
                mpn_addmul_1(exp2 + i, exp3, N - i, cp[i]);
        }
        else
        {
            ulong * t;
            TMP_INIT;
            TMP_START;
            t = TMP_ALLOC(N*sizeof(ulong));

            mpn_mul_1(t, exp3, N, cp[0]);
            for (i = 1; i < cn; i++)
                mpn_addmul_1(t + i, exp3, N - i, cp[i]);

            for (i = 0; i < N; i++)
                exp2[i] = t[i];

            TMP_END;            
        }
    }
}
