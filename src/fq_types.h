/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_TYPES_H
#define FQ_TYPES_H

#include "fmpz_mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef fmpz_poly_t fq_t;
typedef fmpz_poly_struct fq_struct;

typedef struct
{
    fmpz_mod_ctx_t ctxp;

    int sparse_modulus;
    int is_conway; /* whether field was initialized with the Flint Conway tables  (assures primitivity) */

    fmpz * a;
    slong * j;
    slong len;

    fmpz_mod_poly_t modulus;
    fmpz_mod_poly_t inv;

    char * var;
}
fq_ctx_struct;

typedef fq_ctx_struct fq_ctx_t[1];

typedef struct
{
    fq_struct * entries;
    slong r;
    slong c;
    fq_struct ** rows;
}
fq_mat_struct;

typedef fq_mat_struct fq_mat_t[1];

typedef struct
{
    fq_struct * coeffs;
    slong alloc;
    slong length;
}
fq_poly_struct;

typedef fq_poly_struct fq_poly_t[1];

typedef struct
{
    fq_poly_struct * poly;
    slong * exp;
    slong num;
    slong alloc;
}
fq_poly_factor_struct;

typedef fq_poly_factor_struct fq_poly_factor_t[1];

#ifdef __cplusplus
}
#endif

#endif /* FQ_TYPES_H */
