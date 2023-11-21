/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

void
dirichlet_chi_vec(ulong *v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv)
{
    if (2 * nv > G->phi_q)
        dirichlet_chi_vec_loop(v, G, chi, nv);
    else
        dirichlet_chi_vec_primeloop(v, G, chi, nv);
}

void
dirichlet_chi_vec_order(ulong *v, const dirichlet_group_t G, const dirichlet_char_t chi, ulong order, slong nv)
{
    if (2 * nv > G->phi_q)
        dirichlet_chi_vec_loop_order(v, G, chi, order, nv);
    else
        dirichlet_chi_vec_primeloop_order(v, G, chi, order, nv);
}
