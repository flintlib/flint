/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef N_MOD_TYPES_H
#define N_MOD_TYPES_H

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
   ulong nu;
   ulong nn;
   ulong ninv;
   unsigned int norm;
}
n_mod_ctx_struct;

typedef n_mod_ctx_struct n_mod_ctx_t[1];
typedef n_mod_ctx_struct * n_mod_ctx_ptr;
typedef const n_mod_ctx_struct * n_mod_ctx_srcptr;

typedef struct
{
    nn_ptr entries;
    slong r;
    slong c;
    nn_ptr * rows;
}
n_mod_mat_struct;

typedef n_mod_mat_struct n_mod_mat_t[1];
typedef n_mod_mat_struct * n_mod_mat_ptr;
typedef const n_mod_mat_struct * n_mod_mat_srcptr;

typedef struct
{
    nn_ptr coeffs;
    slong alloc;
    slong length;
}
n_mod_poly_struct;

typedef n_mod_poly_struct n_mod_poly_t[1];
typedef n_mod_poly_struct * n_mod_poly_ptr;
typedef const n_mod_poly_struct * n_mod_poly_srcptr;

#ifdef __cplusplus
}
#endif

#endif /* N_MOD_TYPES_H */
