/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"

/* Calculating the GCD is a very costly operation.
  
   Hence, our attempt is to find a single-limbed integer of which we can start
   the process with, because calculating the GCD of {u, n} and {v, 1} is much
   more effective than calculating if of two n-limbed integers.
*/

/* TODO
   * Do we want to utilize mpn_gcd_22 coupled with mpn_divrem_2? Both functions
     should be available in all new versions of GMP. */

void
_fmpz_vec_content_chained(fmpz_t res, const fmpz * vec, slong len, const fmpz_t in)
{
    flint_bitcnt_t exp2;
    mp_size_t limb2;
    slong ix;

    if (in != NULL && !COEFF_IS_MPZ(*in))
    {
        fmpz_set(res, in);
        goto gcd1;
    }

gcd1: /* We will assume that `in` is taken care of in `res`. */
}
