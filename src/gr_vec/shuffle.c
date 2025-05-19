/*
   Copyright (C) 2025 Marc Mezzarobba

   This file is part of FLINT.

   FLINT is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License (LGPL) as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
   */

#include "gr_vec.h"
#include "perm.h"

void
_gr_vec_shuffle(gr_ptr vec, flint_rand_t state, slong len, gr_ctx_t ctx)
{
    slong * perm = _perm_init(len);
    _perm_randtest(perm, len, state);

    _gr_vec_permute(vec, perm, len, ctx);

    _perm_clear(perm);
}
