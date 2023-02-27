/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_mat.h"

#ifdef T
#undef T
#endif

#define T fq
#define CAP_T FQ
#include "fq_mat_templates/set_nmod_mat.c"
#undef CAP_T
#undef T
