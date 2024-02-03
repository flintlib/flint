/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"

void
fmpz_mod_poly_frobenius_powers_clear(fmpz_mod_poly_frobenius_powers_t pow,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong i;

    for (i = 0; i <= pow->len; i++)
       fmpz_mod_poly_clear(pow->pow + i, ctx);

    flint_free(pow->pow);
}
