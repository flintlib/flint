/*
    Copyright (C) 2025 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"
#include "fq_mat.h"

#ifdef T
#undef T
#endif

#define T fq
#define CAP_T FQ
#include "fq_mat_templates/test/t-transpose.c"
#undef CAP_T
#undef T
