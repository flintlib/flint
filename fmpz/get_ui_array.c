/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

/*
    Assuming that "in" is non negative and has a limb count <= out_len,
    write the limbs to "out" and zero extend to "out_len" limbs.
*/

void fmpz_get_ui_array(ulong * out, slong out_len, const fmpz_t in)
{
    slong size = 0;
    FLINT_ASSERT(out_len > 0);

    /* copy limbs */
    if (fmpz_abs_fits_ui(in))
    {
        *out++ = fmpz_get_ui(in);
        size++;
    } else
    {
        __mpz_struct * mpz = COEFF_TO_PTR(*in);
        FLINT_ASSERT(mpz->_mp_size <= out_len);
        while (size < mpz->_mp_size)
            *out++ = mpz->_mp_d[size++];
    }

    /* zero extend to out_len */
    while (size++ < out_len)
        *out++ = UWORD(0);
}
