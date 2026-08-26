/*
    Copyright (C) 2008, 2009, William Hart
    Copyright (C) 2010, 2026 Fredrik Johansson
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"

void fmpz_comb_temp_clear(fmpz_comb_temp_t CT)
{
    flint_free(CT->tmp);
    flint_free(CT->out);
}

void fmpz_comb_clear(fmpz_comb_t C)
{
    flint_mpn_crt_clear(C->crt);
    flint_free(C->crt);
}
