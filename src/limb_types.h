/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef LIMB_TYPES_H
#define LIMB_TYPES_H

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FLINT_MAX_FACTORS_IN_LIMB 15

typedef struct
{
   int num;
   int exp[FLINT_MAX_FACTORS_IN_LIMB];
   ulong p[FLINT_MAX_FACTORS_IN_LIMB];
}
n_factor_t;

typedef struct
{
    slong small_i;
    slong small_num;
    unsigned int * small_primes;

    ulong sieve_a;
    ulong sieve_b;
    slong sieve_i;
    slong sieve_num;
    char * sieve;
}
n_primes_struct;

typedef n_primes_struct n_primes_t[1];

#ifdef __cplusplus
}
#endif

#endif /* LIMB_TYPES_H */
