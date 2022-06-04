/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_IMPL_H
#define FMPZ_IMPL_H

#include <stdio.h>

#ifdef __unix__
#include <unistd.h> /* sysconf */
#endif

#if defined(_WIN32) || defined(WIN32)
#include <windows.h> /* GetSytemInfo */
#endif

#if defined(_MSC_VER) && FLINT_USES_PTHREAD
#include <atomic.h>
#endif

#include "flint.h"

#if FLINT_USES_PTHREAD
#include <pthread.h>
#endif

#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "mpfr.h"
#include "long_extras.h"
#include "mpn_extras.h"
#include "fmpz_factor.h"
#include "aprcl.h"

#if FLINT64   /* 2^53 */
#define DOUBLE_MAX WORD(9007199254740992)
#define DOUBLE_MIN (-DOUBLE_MAX)
#else
#define DOUBLE_MAX COEFF_MAX
#define DOUBLE_MIN COEFF_MIN
#endif

#endif
