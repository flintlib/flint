/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

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
  1UL, 1UL, 2UL, 6UL, 24UL, 120UL, 720UL, 5040UL, 40320UL, 362880UL,
  3628800UL, 39916800UL, 479001600UL,
#if FLINT64
  6227020800UL, 87178291200UL, 1307674368000UL, 20922789888000UL,
  355687428096000UL, 6402373705728000UL, 121645100408832000UL,
  2432902008176640000UL,
#endif
};

void fmpz_fac_ui(fmpz_t f, ulong n)
{
    if (n < FLINT_NUM_TINY_FACTORIALS)
        fmpz_set_ui(f, flint_tiny_factorials[n]);
    else
        mpz_fac_ui(_fmpz_promote(f), n);
}
