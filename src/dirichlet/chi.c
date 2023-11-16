/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

ulong
dirichlet_chi(const dirichlet_group_t G, const dirichlet_char_t chi, ulong n)
{
    if (n_gcd(G->q, n) > 1)
    {
        return DIRICHLET_CHI_NULL;
    }
    else
    {
        ulong v;
        dirichlet_char_t x;
        dirichlet_char_init(x, G);
        dirichlet_char_log(x, G, n);

        v = dirichlet_pairing_char(G, chi, x);

        dirichlet_char_clear(x);
        return v;
    }
}
