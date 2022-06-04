/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MPN_EXTRAS_IMPL_H
#define MPN_EXTRAS_IMPL_H

#include "flint.h"

#if FLINT_REENTRANT && !FLINT_USES_TLS
#include <pthread.h>
#endif

#include "ulong_extras.h"
#include "mpn_extras.h"

#endif
