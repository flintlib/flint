/*
    Copyright (C) 2020 Fredrik Johansson
    Copyright (C) 2025 Andrii Yanovets

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

int
fmpz_mod_mpoly_buchberger_with_limits(fmpz_mod_mpoly_vec_t G, const fmpz_mod_mpoly_vec_t F,
    slong ideal_len_limit, slong poly_len_limit, int (*solver)(fmpz_mod_mpoly_vec_t, const fmpz_mod_mpoly_vec_t,
    slong, slong), const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;

    success = solver(G, F, ideal_len_limit, poly_len_limit);

    return success;
}
