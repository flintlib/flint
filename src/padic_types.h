/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef PADIC_TYPES_H
#define PADIC_TYPES_H

#include "flint.h"

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

#define padic_val(x)   ((x)->v)
#define padic_prec(x)  ((x)->N)

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

#ifdef __cplusplus
}
#endif

#endif /* PADIC_TYPES_H */
