/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_IMPL_H
#define FMPZ_IMPL_H

#include <gmp.h>
#include "fmpz_types.h"

void _flint_mpz_addmul_large(mpz_ptr z, mpz_srcptr x, mpz_srcptr y, int negate);
void _fmpz_CRT(fmpz_t out, const fmpz_t r1, const fmpz_t m1, const fmpz_t r2, const fmpz_t m2, const fmpz_t m1m2, fmpz_t c, int sign);
char * fmpz_get_str_bsplit_threaded(char * s, const fmpz_t f);

#endif
