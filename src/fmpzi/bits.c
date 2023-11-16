/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"

slong fmpzi_bits(const fmpzi_t x)
{
    ulong v;
    fmpz a, b;

    a = *fmpzi_realref(x);
    b = *fmpzi_imagref(x);

    if (!COEFF_IS_MPZ(a))
    {
        if (!COEFF_IS_MPZ(b))
        {
            v = FLINT_ABS(a) | FLINT_ABS(b);
            return FLINT_BIT_COUNT(v);
        }
        else
        {
            return mpz_sizeinbase(COEFF_TO_PTR(b), 2);
        }
    }
    else
    {
        if (!COEFF_IS_MPZ(b))
        {
            return mpz_sizeinbase(COEFF_TO_PTR(a), 2);
        }
        else
        {
            __mpz_struct *z1, *z2;
            slong s1, s2;

            z1 = COEFF_TO_PTR(a);
            z2 = COEFF_TO_PTR(b);

            s1 = FLINT_ABS(z1->_mp_size);
            s2 = FLINT_ABS(z2->_mp_size);

            if (s1 == s2)
            {
                v = z1->_mp_d[s1 - 1] | z2->_mp_d[s1 - 1];
            }
            else if (s1 > s2)
            {
                v = z1->_mp_d[s1 - 1];
            }
            else
            {
                v = z2->_mp_d[s2 - 1];
                s1 = s2;
            }

            return (s1 - 1) * FLINT_BITS + FLINT_BIT_COUNT(v);
        }
    }
}
