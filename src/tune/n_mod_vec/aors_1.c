/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_mod_vec.h"

#undef N_MOD_VEC_ADD_METHOD
#undef N_MOD_VEC_SUB_METHOD
#define TUNE_PROGRAM 1

#define N_MOD_VEC_ADD_METHOD 1
#define N_MOD_VEC_SUB_METHOD 1

#define _n_mod_vec_add _n_mod_vec_add_1
#define _n_mod_vec_sub _n_mod_vec_sub_1

#include "n_mod_vec/aors.c"
