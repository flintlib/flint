/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPQ_MAT_IMPL_H
#define FMPQ_MAT_IMPL_H

#include <stdio.h>
#include "fmpq_vec.h"
#include "fmpq_mat.h"

/* defined in solve_multi_mod.c */
int
_fmpq_mat_check_solution_fmpz_mat(const fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B);

#endif
