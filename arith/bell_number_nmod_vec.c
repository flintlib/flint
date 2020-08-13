/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

void
arith_bell_number_nmod_vec(mp_ptr b, slong n, nmod_t mod)
{
    if (n < 2000 || mod.n <= n)
        arith_bell_number_nmod_vec_recursive(b, n, mod);
    else
        arith_bell_number_nmod_vec_series(b, n, mod);
}
