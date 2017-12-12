/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/*
    ensure that bits is a valid (or good) number of bits
*/
slong mpoly_fix_bits(slong bits, const mpoly_ctx_t mctx) {
    slong nfields = mctx->nfields;
    if (bits <= FLINT_BITS)
    {
        bits = FLINT_MAX(bits, MPOLY_MIN_BITS);

        while (bits < FLINT_BITS &&   (nfields - 1)/(FLINT_BITS/(bits    ))
                                   == (nfields - 1)/(FLINT_BITS/(bits + 1))
              )
        {
            bits++;
        }
    } else {

        bits = (bits + FLINT_BITS - 1)/FLINT_BITS*FLINT_BITS;
    }

    return bits;
}
