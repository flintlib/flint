/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_TYPES_H
#define ACB_TYPES_H

#include "arb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    arb_struct real;
    arb_struct imag;
}
acb_struct;

typedef acb_struct acb_t[1];
typedef acb_struct * acb_ptr;
typedef const acb_struct * acb_srcptr;

typedef struct
{
    acb_ptr entries;
    slong r;
    slong c;
    acb_ptr * rows;
}
acb_mat_struct;

typedef acb_mat_struct acb_mat_t[1];

typedef struct
{
    acb_ptr coeffs;
    slong alloc;
    slong length;
}
acb_poly_struct;

typedef acb_poly_struct acb_poly_t[1];

#ifdef __cplusplus
}
#endif

#endif /* ACB_TYPES_H */
