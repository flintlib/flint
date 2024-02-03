/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_TYPES_H
#define FMPZ_TYPES_H

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    int sign;
    fmpz * p;
    ulong * exp;
    slong alloc;
    slong num;
}
fmpz_factor_struct;

typedef fmpz_factor_struct fmpz_factor_t[1];

typedef struct
{
   mp_ptr dinv;
   slong n;
   flint_bitcnt_t norm;
} fmpz_preinvn_struct;

typedef fmpz_preinvn_struct fmpz_preinvn_t[1];

typedef struct
{
    fmpz * coeffs;
    slong alloc;
    slong length;
}
fmpz_poly_struct;

typedef fmpz_poly_struct fmpz_poly_t[1];

typedef struct
{
    fmpz c;
    fmpz_poly_struct *p;
    slong *exp;
    slong num;
    slong alloc;
}
fmpz_poly_factor_struct;

typedef fmpz_poly_factor_struct fmpz_poly_factor_t[1];

typedef struct
{
    fmpz * entries;
    slong r;
    slong c;
    fmpz ** rows;
}
fmpz_mat_struct;

typedef fmpz_mat_struct fmpz_mat_t[1];

typedef struct
{
    fmpz_poly_struct * entries;
    slong r;
    slong c;
    fmpz_poly_struct ** rows;
}
fmpz_poly_mat_struct;

typedef fmpz_poly_mat_struct fmpz_poly_mat_t[1];

typedef struct
{
   fmpz * coeffs; /* alloc fmpzs */
   ulong * exps;
   slong alloc;
   slong length;
   flint_bitcnt_t bits;     /* number of bits per exponent */
}
fmpz_mpoly_struct;

typedef fmpz_mpoly_struct fmpz_mpoly_t[1];

typedef struct
{
    fmpz_t constant;
    fmpz_t constant_den;        /* should be one after normal operations */
    fmpz_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
}
fmpz_mpoly_factor_struct;

typedef fmpz_mpoly_factor_struct fmpz_mpoly_factor_t[1];

typedef struct
{
    fmpz_poly_struct *num;
    fmpz_poly_struct *den;
}
fmpz_poly_q_struct;

typedef fmpz_poly_q_struct fmpz_poly_q_t[1];

typedef struct
{
    fmpz_mpoly_struct num;
    fmpz_mpoly_struct den;
}
fmpz_mpoly_q_struct;

typedef fmpz_mpoly_q_struct fmpz_mpoly_q_t[1];

typedef struct
{
    fmpz a;
    fmpz b;
}
fmpzi_struct;

typedef fmpzi_struct fmpzi_t[1];

#ifdef __cplusplus
}
#endif

#endif /* FMPZ_TYPES_H */
