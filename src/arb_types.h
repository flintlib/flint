/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_TYPES_H
#define ARB_TYPES_H

#include "arf_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    fmpz exp;
    mp_limb_t man;
}
mag_struct;

typedef mag_struct mag_t[1];
typedef mag_struct * mag_ptr;
typedef const mag_struct * mag_srcptr;

typedef struct
{
    arf_struct mid;
    mag_struct rad;
}
arb_struct;

typedef arb_struct arb_t[1];
typedef arb_struct * arb_ptr;
typedef const arb_struct * arb_srcptr;

typedef struct
{
    arb_ptr entries;
    slong r;
    slong c;
    arb_ptr * rows;
}
arb_mat_struct;

typedef arb_mat_struct arb_mat_t[1];

typedef struct
{
    arb_ptr coeffs;
    slong alloc;
    slong length;
}
arb_poly_struct;

typedef arb_poly_struct arb_poly_t[1];

#ifdef __cplusplus
}
#endif

#endif /* ARB_TYPES_H */
