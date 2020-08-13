/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* !!! this file DOES need to change with new orderings */

/* bits required to represent the generator of index var */
flint_bitcnt_t mpoly_gen_bits_required(slong var, const mpoly_ctx_t mctx)
{
    return MPOLY_MIN_BITS; /* 2 bits suffice in every case */
}
