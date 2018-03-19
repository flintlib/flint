/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "assert.h"

void fmpz_get_ui_array(ulong * out, slong out_len, const fmpz_t in) {
    slong size = 0;
    assert(out_len > 0);

    /* copy limbs */
    if (fmpz_abs_fits_ui(in)) {
        *out++ = fmpz_get_ui(in);
        size++;
    } else {
        __mpz_struct * mpz = COEFF_TO_PTR(*in);
        assert(mpz->_mp_size <= out_len);
        while (size < mpz->_mp_size)
            *out++ = mpz->_mp_d[size++];
    }

    /* zero extend to out_len */
    while (size++ < out_len)
        *out++ = 0;
}
