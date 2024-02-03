/*
    Copyright (C) 2011, 2012, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef PADIC_TYPES_H
#define PADIC_TYPES_H

#include "fmpz_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    fmpz u;
    slong v;
    slong N;
} padic_struct;

typedef padic_struct padic_t[1];

enum padic_print_mode
{
    PADIC_TERSE,
    PADIC_SERIES,
    PADIC_VAL_UNIT
};

typedef struct
{
    fmpz_t p;
    double pinv;
    fmpz * pow;
    slong min;
    slong max;
    enum padic_print_mode mode;
} padic_ctx_struct;

typedef padic_ctx_struct padic_ctx_t[1];

typedef struct
{
    slong n;
    fmpz * pow;
} padic_inv_struct;

typedef padic_inv_struct padic_inv_t[1];

typedef struct
{
    fmpz_mat_struct mat;
    slong val;
    slong N;
} padic_mat_struct;

typedef padic_mat_struct padic_mat_t[1];

typedef struct
{
    fmpz * coeffs;
    slong alloc;
    slong length;
    slong val;
    slong N;
} padic_poly_struct;

typedef padic_poly_struct padic_poly_t[1];

#ifdef __cplusplus
}
#endif

#endif /* PADIC_TYPES_H */
