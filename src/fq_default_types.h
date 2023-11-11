/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_DEFAULT_TYPES_H
#define FQ_DEFAULT_TYPES_H

#include "fq_types.h"
#include "fq_nmod_types.h"
#include "fq_zech_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FQ_DEFAULT_FQ_ZECH  1
#define FQ_DEFAULT_FQ_NMOD  2
#define FQ_DEFAULT_FQ       3
#define FQ_DEFAULT_NMOD     4
#define FQ_DEFAULT_FMPZ_MOD 5

typedef union fq_default_struct
{
    fq_t fq;
    fq_nmod_t fq_nmod;
    fq_zech_t fq_zech;
    ulong nmod;
    fmpz_t fmpz_mod;
} fq_default_struct;

typedef fq_default_struct fq_default_t[1];

typedef struct
{
    int type;
    union ctx
    {
        fq_ctx_t fq;
        fq_nmod_ctx_t fq_nmod;
        fq_zech_ctx_t fq_zech;
        struct {
            nmod_t mod;
            mp_limb_t a;    /* minpoly is x - a */
        } nmod;
        struct {
            fmpz_mod_ctx_t mod;
            fmpz_t a;       /* minpoly is x - a */
        } fmpz_mod;
    } ctx;
} fq_default_ctx_struct;

typedef fq_default_ctx_struct fq_default_ctx_t[1];

typedef union fq_default_mat_struct
{
    fq_mat_t fq;
    fq_nmod_mat_t fq_nmod;
    fq_zech_mat_t fq_zech;
    nmod_mat_t nmod;
    fmpz_mod_mat_t fmpz_mod;
} fq_default_mat_struct;

typedef fq_default_mat_struct fq_default_mat_t[1];

typedef union fq_default_poly_struct
{
    fq_poly_t fq;
    fq_nmod_poly_t fq_nmod;
    fq_zech_poly_t fq_zech;
    nmod_poly_t nmod;
    fmpz_mod_poly_t fmpz_mod;
}
fq_default_poly_struct;

typedef fq_default_poly_struct fq_default_poly_t[1];

#ifdef __cplusplus
}
#endif

#endif /* FQ_DEFAULT_TYPES_H */
