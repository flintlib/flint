/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ULONG_EXTRAS_IMPL_H
#define ULONG_EXTRAS_IMPL_H

#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <pthread.h>
#include "mpfr.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "fmpz.h"

/* NOTE: These lie in rootrem.c */

/* A table of precomputed inverses of values from 1 to 64  
   inv_table[n-1] = 1/n for all n in range[1, 64] */
extern const double inv_table[64];

/* This table has the max possible base for a given root. 
   max_base[n-1] = UWORD_MAX^(1/n) for n in range [1, FLINT_BITS] */
extern const mp_limb_t max_base[FLINT_BITS];

#endif
