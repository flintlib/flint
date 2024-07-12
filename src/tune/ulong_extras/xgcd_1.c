/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint-mparam.h"
#include "ulong_extras.h"

#undef N_GCDEXT_METHOD
#define TUNE_PROGRAM 1

#define N_GCDEXT_METHOD 1

#define n_xgcd n_xgcd_1

#include "ulong_extras/xgcd.c"
