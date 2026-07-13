/*
    Copyright (C) 2016 Pascal Molin
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "gr_dft.h"

/* Thin wrappers around the gr_dft module (see gr_dft/dirichlet.c),
   which computes the product DFT over the Conrey components in
   fixed-point arithmetic with rigorous error bounds when possible. */

/* dft, lexicographic conrey indexing, array size G->phi_q */
void
acb_dirichlet_dft_index(acb_ptr w, acb_srcptr v, const dirichlet_group_t G,
        slong prec)
{
    gr_dft_acb_dirichlet_index(w, v, G, prec);
}

/* dft, number indexing, array size G->q */
void
acb_dirichlet_dft(acb_ptr w, acb_srcptr v, const dirichlet_group_t G,
        slong prec)
{
    gr_dft_acb_dirichlet(w, v, G, prec);
}
