/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_ZECH_TYPES_H
#define FQ_ZECH_TYPES_H

#include "fq_nmod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    mp_limb_t value;
}
fq_zech_struct;

typedef fq_zech_struct fq_zech_t[1];

typedef struct
{
    ulong qm1;     /* q - 1 */
    ulong qm1o2;   /* (q - 1) / 2 or 1 when p == 2 */
    ulong qm1opm1; /* (q - 1) / (p - 1) */
    ulong p;
    double ppre;
    ulong prime_root;       /* primitive root for prime subfield */
    ulong * zech_log_table;
    ulong * prime_field_table;
    ulong * eval_table;

    fq_nmod_ctx_struct * fq_nmod_ctx;
    int owns_fq_nmod_ctx;
    int is_conway; /* whether field was generated using Flint Conway tables (assures primitivity) */
}
fq_zech_ctx_struct;

typedef fq_zech_ctx_struct fq_zech_ctx_t[1];

typedef struct
{
    fq_zech_struct * entries;
    slong r;
    slong c;
    fq_zech_struct ** rows;
}
fq_zech_mat_struct;

typedef fq_zech_mat_struct fq_zech_mat_t[1];

typedef struct
{
    fq_zech_struct * coeffs;
    slong alloc;
    slong length;
}
fq_zech_poly_struct;

typedef fq_zech_poly_struct fq_zech_poly_t[1];

typedef struct
{
    fq_zech_poly_struct * poly;
    slong * exp;
    slong num;
    slong alloc;
}
fq_zech_poly_factor_struct;

typedef fq_zech_poly_factor_struct fq_zech_poly_factor_t[1];

#ifdef __cplusplus
}
#endif

#endif /* FQ_ZECH_TYPES_H */
