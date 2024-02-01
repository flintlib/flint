/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPQ_TYPES_H
#define FMPQ_TYPES_H

#include "fmpz_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    fmpq * entries;
    slong r;
    slong c;
    fmpq ** rows;
}
fmpq_mat_struct;

typedef fmpq_mat_struct fmpq_mat_t[1];

typedef struct
{
    fmpz * coeffs;
    slong alloc;
    slong length;
    fmpz_t den;
}
fmpq_poly_struct;

typedef fmpq_poly_struct fmpq_poly_t[1];

/*
    A polynomial f is represented as
        content * zpoly,
    where zpoly should have positive leading coefficient and trivial content.
    If f is zero, then the representation should have
        content = 0 and zpoly = 0
*/

typedef struct
{                       /* non zero case:                   |  zero case: */
    fmpq_t content;     /* positive or negative content     |  zero       */
    fmpz_mpoly_t zpoly; /* contentless poly, lc is positive |  zero       */
}
fmpq_mpoly_struct;

typedef fmpq_mpoly_struct fmpq_mpoly_t[1];

typedef struct
{
    fmpq_t constant;
    fmpq_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
}
fmpq_mpoly_factor_struct;

typedef fmpq_mpoly_factor_struct fmpq_mpoly_factor_t[1];

#ifdef __cplusplus
}
#endif

#endif /* FMPQ_TYPES_H */
