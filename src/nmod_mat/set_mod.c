/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_mat.h"

void nmod_mat_set_mod(nmod_mat_t mat, ulong n)
{
    mat->mod.n = n;
    mat->mod.norm = flint_clz(n);
    mat->mod.ninv = n_preinvert_limb_prenorm(n << mat->mod.norm);
}
