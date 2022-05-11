/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_vec.h"
#include "fq_zech_poly_mini.h"

#ifdef T
#undef T
#endif

#define T fq_zech
#define CAP_T FQ_ZECH
#include "fq_poly_templates/rem.c"
#undef CAP_T
#undef T
