/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "assert.h"

void fmpz_set_ui_array(fmpz_t out, const ulong * in, slong in_len) {

    slong size = in_len;

    /* find end of zero extension */
    while (size > 1 && in[size - 1] == 0)
        size--;

    /* copy limbs */
    if (size == 1) {
        fmpz_set_ui(out, in[0]);
    } else {
        __mpz_struct * mpz = _fmpz_promote(out);
        mpz_realloc2(mpz, FLINT_BITS*size);
        mpn_copyi(mpz->_mp_d, in, size);
        mpz->_mp_size = size;
    }
}
