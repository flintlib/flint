/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Here we mock a few GMP types so we do not have to include gmp.h in every
 * source file. */

#ifndef MOCK_GMP_TYPES_H
#define MOCK_GMP_TYPES_H

#include "flint-config.h"

#ifdef __cplusplus
extern "C" {
#endif

#define GMP_LIMB_BITS FLINT_GMP_LIMB_BITS

typedef FLINT_MP_LIMB_T mp_limb_t;
typedef FLINT_MP_LIMB_SIGNED_T mp_limb_signed_t;
typedef long int mp_size_t;

typedef mp_limb_t *	mp_ptr;
typedef const mp_limb_t * mp_srcptr;

typedef struct
{
    int _mp_alloc;
    int _mp_size;
    mp_limb_t *_mp_d;
}
__mpz_struct;

typedef __mpz_struct mpz_t[1];
typedef __mpz_struct * mpz_ptr;
typedef const __mpz_struct * mpz_srcptr;

typedef enum
{
  GMP_RAND_ALG_DEFAULT = 0,
  GMP_RAND_ALG_LC = GMP_RAND_ALG_DEFAULT
}
gmp_randalg_t;

typedef struct
{
    mpz_t _mp_seed;
    gmp_randalg_t _mp_alg;
    union { void *_mp_lc; } _mp_algdata;
}
__gmp_randstate_struct;

typedef __gmp_randstate_struct gmp_randstate_t[1];

#ifdef __cplusplus
}
#endif

#endif /* MOCK_GMP_TYPES_H */
