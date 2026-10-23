/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MP_REAL_DENSE_H
#define MP_REAL_DENSE_H

/* What the dense-table per-size kernels (the *_opt_<n>.c files emitted
   from dev/gen_mp_real_dense_bodies.py) need besides impl.h: their
   tables (the first two levels from the 21-limb tables they share with
   the run-time kernels) and the small reciprocal of mpn_extras/recip_small.h.
   Included from impl.h. */

#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "dense_tables.h"
#include "dense_rt_tables.h"
#include "mpn_extras/recip_small.h"

#endif
