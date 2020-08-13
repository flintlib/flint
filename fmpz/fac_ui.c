/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

#if FLINT64
#define FLINT_NUM_TINY_FACTORIALS 21
#else
#define FLINT_NUM_TINY_FACTORIALS 13
#endif

const mp_limb_t flint_tiny_factorials[] =
{
  UWORD(1), UWORD(1), UWORD(2), UWORD(6), UWORD(24), UWORD(120), UWORD(720), UWORD(5040), UWORD(40320), UWORD(362880),
  UWORD(3628800), UWORD(39916800), UWORD(479001600),
#if FLINT64
  UWORD(6227020800), UWORD(87178291200), UWORD(1307674368000), UWORD(20922789888000),
  UWORD(355687428096000), UWORD(6402373705728000), UWORD(121645100408832000),
  UWORD(2432902008176640000),
#endif
};

void fmpz_fac_ui(fmpz_t f, ulong n)
{
    if (n < FLINT_NUM_TINY_FACTORIALS)
        fmpz_set_ui(f, flint_tiny_factorials[n]);
    else
        flint_mpz_fac_ui(_fmpz_promote(f), n);
}
