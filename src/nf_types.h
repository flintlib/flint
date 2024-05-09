/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NF_TYPES_H
#define NF_TYPES_H

#include "fmpq_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* macros ********************************************************************/

#define NF_GENERIC 0
#define NF_MONIC 1
#define NF_LINEAR 2
#define NF_QUADRATIC 4
#define NF_GAUSSIAN 8

/* number field **************************************************************/

typedef struct
{
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
}
nf_struct;

typedef nf_struct nf_t[1];

/* number field element ******************************************************/

typedef struct /* element of a linear number field */
{
   fmpz_t num;
   fmpz_t den;
}
lnf_elem_struct;

typedef lnf_elem_struct lnf_elem_t[1];

typedef struct /* element of a quadratic number field */
{
   fmpz num[3]; /* extra coeff for delayed reduction */
   fmpz_t den;
}
qnf_elem_struct;

typedef qnf_elem_struct qnf_elem_t[1];

typedef union /* element in a number field (specified by an nf_t) */
{
   fmpq_poly_t elem; /* general case */
   lnf_elem_t lelem; /* linear number field */
   qnf_elem_t qelem; /* quadratic number field */
}
nf_elem_struct;

typedef nf_elem_struct nf_elem_t[1];

#ifdef __cplusplus
}
#endif

#endif /* NF_TYPES_H */
