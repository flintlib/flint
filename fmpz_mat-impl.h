/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MAT_IMPL_H
#define FMPZ_MAT_IMPL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "thread_support.h"
#include "perm.h"
#include "fmpz_mat.h"
#include "fmpq_vec.h"
#include "fmpq_mat.h"
#include "fq_nmod.h"
#include "fmpz_factor.h"
#include "fft.h"

#define __FLINT_LOG2E  1.44269504088896340736  /* log2(e) */

static __inline__ long double _log2(const long double x)
{
    return log(x) * __FLINT_LOG2E;
}

#endif
