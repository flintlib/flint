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
dirichlet_pairing(const dirichlet_group_t G, ulong m, ulong n)
{
    ulong x;
    dirichlet_char_t a, b;

    if (n_gcd(G->q, m) > 1 || n_gcd(G->q, n) > 1)
        return DIRICHLET_CHI_NULL;

    dirichlet_char_init(a, G);
    dirichlet_char_init(b, G);
    dirichlet_char_log(a, G, m);
    dirichlet_char_log(b, G, n);

    x = dirichlet_pairing_char(G, a, b);

    dirichlet_char_clear(a);
    dirichlet_char_clear(b);

    return x;
}
