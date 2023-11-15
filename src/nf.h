/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NF_H
#define NF_H

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
   fmpq_poly_t pol;  /* defining polynomial */
   union {
      /* insert any precomputed inverse for zz case here */
      fmpz_preinvn_t qq; /* precomputed inverse for leading coeff of num(pol), QQ case */
   } pinv;
   union { /* powers of the generator mod pol */
      fmpq_poly_powers_precomp_t qq;
      fmpz_poly_powers_precomp_t zz;
   } powers;
   fmpq_poly_t traces; /* S_k = sum_i \theta_i^k for k = 0, 1, 2, ..., (n-1) */
   ulong flag;       /* 1 = pol monic over ZZ, 2, = linear, 4 = quadratic field */
} nf_struct;

typedef nf_struct nf_t[1];

#define NF_POWERS_CUTOFF 30 /* maximum length of pol where we precompute powers */

#define NF_GENERIC 0
#define NF_MONIC 1
#define NF_LINEAR 2
#define NF_QUADRATIC 4
#define NF_GAUSSIAN 8

/******************************************************************************

    Initialisation

******************************************************************************/

void nf_init(nf_t nf, const fmpq_poly_t pol);

void nf_init_randtest(nf_t nf, flint_rand_t state, slong len,  mp_bitcnt_t bits_in);

void nf_clear(nf_t nf);

void nf_print(const nf_t nf);

#ifdef __cplusplus
}
#endif

#endif
