/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"

slong
_fmpz_vec_max_bits(const fmpz * vec, slong len)
{
    slong i, sign, max_limbs;
    ulong max_limb;
    slong limbs;

    sign = 1;
    max_limb = 0;

    for (i = 0; i < len; i++)
    {
        fmpz c = vec[i];

        if (c >= 0)
        {
            if (c > COEFF_MAX)
                goto bignum;
            max_limb |= c;
        }
        else
        {
            if (c < COEFF_MIN)
                goto bignum;
            max_limb |= -c;
            sign = -1;
        }
    }
    return sign * FLINT_BIT_COUNT(max_limb);

bignum:
    max_limbs = 1;

    for ( ; i < len; i++)
    {
        fmpz c = vec[i];

        if (COEFF_IS_MPZ(c))
        {
            mpz_ptr z = COEFF_TO_PTR(c);
            limbs = z->_mp_size;

            if (limbs < 0)
            {
                sign = -1;
                limbs = -limbs;
            }

            if (limbs == max_limbs)
                max_limb |= z->_mp_d[limbs - 1];
            else if (limbs > max_limbs)
            {
                max_limb = z->_mp_d[limbs - 1];
                max_limbs = limbs;
            }
        }
        else if (c < 0)
            sign = -1;
    }
    return sign * ((max_limbs - 1) * FLINT_BITS + FLINT_BIT_COUNT(max_limb));
}

slong
_fmpz_vec_max_bits_ref(const fmpz * vec, slong len)
{
    slong i, bits, max_bits = 0, sign = 1;

    for (i = 0; i < len; i++)
    {
        bits = fmpz_bits(vec + i);
        if (bits > max_bits)
            max_bits = bits;
        if (fmpz_sgn(vec + i) < 0)
            sign = WORD(-1);
    }

    return max_bits * sign;
}
