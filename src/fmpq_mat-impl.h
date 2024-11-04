/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPQ_MAT_IMPL_H
#define FMPQ_MAT_IMPL_H

#include "fmpq_types.h"

#ifdef __cplusplus
extern "C" {
#endif

int
_fmpq_mat_check_solution_fmpz_mat(const fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B);

#ifdef __cplusplus
}
#endif

#endif /* FMPQ_MAT_IMPL_H */
