/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_poly.h"



#ifdef T
#undef T
#endif

#define T fq
#define CAP_T FQ
#include "fq_poly_templates/add_si.c"
#undef CAP_T
#undef T
