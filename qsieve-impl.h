/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef QSIEVE_IMPL_H
#define QSIEVE_IMPL_H

#if (!defined (__WIN32) || defined(__CYGWIN__)) && !defined(_MSC_VER)
#include <unistd.h>
#endif
#if defined (__WIN32) && !defined(__CYGWIN__)
#include <windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "fmpz_vec.h"
#include "qsieve.h"

#endif
