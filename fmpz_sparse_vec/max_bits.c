/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_sparse_vec.h"


slong fmpz_sparse_vec_max_bits(const fmpz_sparse_vec_t v)
{
    slong i, sign, max_limbs;
    mp_limb_t max_limb;
    mp_size_t limbs;

    sign = 1;
    max_limb = 0;

    for (i = 0; i < v->nnz; i++)
    {
        fmpz c = *(v->entries[i].val);

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

    for ( ; i < v->nnz; i++)
    {
        fmpz c = *(v->entries[i].val);

        if (COEFF_IS_MPZ(c))
        {
            __mpz_struct * z = COEFF_TO_PTR(c);
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