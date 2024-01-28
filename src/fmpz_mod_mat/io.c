/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz_mod_mat.h"

int fmpz_mod_mat_fprint(FILE * file, const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx) { return fmpz_mat_fprint(file, mat); }
int fmpz_mod_mat_fprint_pretty(FILE * file, const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx) { return fmpz_mat_fprint_pretty(file, mat); }
int fmpz_mod_mat_print(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx) { return fmpz_mat_print(mat); }
void fmpz_mod_mat_print_pretty(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx) { fmpz_mat_print_pretty(mat); }
