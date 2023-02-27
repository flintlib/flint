/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#define ulong ulongxx /* interferes with system includes */
#include <math.h>
#undef ulong
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

/* TODO: improve the implementation of this function */
ulong n_CRT(ulong A1, ulong M1, ulong A2, ulong M2)
{
    ulong R;
    fmpz_t r, a1, m1, a2, m2;

    fmpz_init(r);
    fmpz_init_set_ui(a1, A1);
    fmpz_init_set_ui(m1, M1);
    fmpz_init_set_ui(a2, A2);
    fmpz_init_set_ui(m2, M2);

    fmpz_CRT(r, a1, m1, a2, m2, 0);
    FLINT_ASSERT(fmpz_abs_fits_ui(r));
    R = fmpz_get_ui(r);   

    fmpz_clear(r);
    fmpz_clear(a1);
    fmpz_clear(m1);
    fmpz_clear(a2);
    fmpz_clear(m2);

    return R;
}
