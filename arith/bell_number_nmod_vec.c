/*
    Copyright (C) 2011, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

void
arith_bell_number_nmod_vec(mp_ptr b, slong len, nmod_t mod)
{
    if (len < 300)
    {
        arith_bell_number_nmod_vec_recursive(b, len, mod);
    }
    else
    {
        if (mod.n >= len && arith_bell_number_nmod_vec_series(b, len, mod))
            return;

        if (len < 500 + NMOD_BITS(mod) * NMOD_BITS(mod))
            arith_bell_number_nmod_vec_recursive(b, len, mod);
        else
            arith_bell_number_nmod_vec_ogf(b, len, mod);
    }
}
