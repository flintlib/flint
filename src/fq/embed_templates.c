/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"
#include "fq.h"
#include "fq_embed.h"

#ifdef T
#undef T
#endif
#ifdef B
#undef B
#endif

#define T fq
#define CAP_T FQ
#define B fmpz_mod

#include "fq_embed_templates/composition_matrix.c"

#undef B
#undef CAP_T
#undef T
