/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_sparse_mat.h"

#ifdef T
#undef T
#endif

#define T fq_nmod
#define CAP_T FQ_NMOD
#include "fq_sparse_mat_templates/from_entries.c"
#undef CAP_T
#undef T
