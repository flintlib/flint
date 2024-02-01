/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2018 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_factor.h"
#include "fq_zech.h"

#ifdef T
#undef T
#endif

#define T fq_zech
#define CAP_T FQ_ZECH

#include "fq_templates/div.c"
#include "fq_templates/is_invertible.c"
#include "fq_templates/is_invertible_f.c"
#include "fq_templates/multiplicative_order.c"

#undef CAP_T
#undef T
