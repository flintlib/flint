/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2023, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "gmpcompat.h"
#include "fmpz.h"
#include "fmpz/impl.h"
#include "mpn_extras.h"

char * fmpz_get_str(char * str, int b, const fmpz_t f)
{
    FLINT_ASSERT(b >= 2 && b <= 62);

    if (!COEFF_IS_MPZ(*f))
    {
        fmpz c;
        ulong d;
        c = *f;

        d = FLINT_ABS(c);

        /* Need a special case for zero, which may as well handle
           single digits. */
        if ((slong) d < FLINT_MIN(b, 10))
        {
            if (str == NULL)
                str = flint_malloc(3);
            str[0] = '-';
            str[0 + (c < 0)] = d + '0';
            str[1 + (c < 0)] = '\0';
        }
        else if (b == 10)
        {
            unsigned char tmp[FLINT_BITS + 3];
            slong i, len;
            /* The compiler might generate faster code for 32-bit divisions */
            unsigned int dl;

            len = 0;
#if FLINT_BITS == 64
            while (d >= (UWORD(1) << 32))
            {
                tmp[len] = d % 10;
                d /= 10;
                len++;
            }
#endif
            dl = d;
            while (dl != 0)
            {
                tmp[len] = dl % 10;
                dl /= 10;
                len++;
            }

            if (str == NULL)
                str = flint_malloc(len + 2);

            str[0] = '-';
            for (i = 0; i < len; i++)
                str[i + (c < 0)] = tmp[len - 1 - i] + '0';
            str[len + (c < 0)] = '\0';

            return str;
        }
        else
        {
            mpz_t z;

            z->_mp_d = &d;
            z->_mp_alloc = 1;
            z->_mp_size = (c > 1) ? 1 : -1;

            if (str == NULL)
                str = flint_malloc(mpz_sizeinbase(z, b) + 2);
            str = mpz_get_str(str, b, z);
        }
    }
    else
    {
        __mpz_struct * z;
        z = COEFF_TO_PTR(*f);
        str = flint_mpn_get_str(str, b, z->_mp_d, FLINT_ABS(z->_mp_size), z->_mp_size < 0);
    }

    return str;
}

