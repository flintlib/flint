/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2025 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdarg.h>
#include "flint.h"

char flint_version[] = FLINT_VERSION;

void flint_va_end(va_list vlist) { va_end(vlist); }
