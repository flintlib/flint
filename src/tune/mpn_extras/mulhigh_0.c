/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint-mparam.h"
#include "mpn_extras.h"
#include "tune.h"

#undef FLINT_MPN_MULHIGH_K_TAB_SIZE
#undef FLINT_MPN_MULHIGH_K_TAB
#define TUNE_PROGRAM 1

#define FLINT_MPN_MULHIGH_K_TAB_SIZE FLINT_MPN_MULHIGH_K_TAB_MAX_SIZE
#define flint_mpn_mulhigh_k_tab flint_mpn_mulhigh_k_tab_0
#define _flint_mpn_mulhigh_n_mulders_recursive _flint_mpn_mulhigh_n_mulders_recursive_0
#define _flint_mpn_mulhigh_n_mulders _flint_mpn_mulhigh_n_mulders_0

#include "mpn_extras/mulhigh.c"
