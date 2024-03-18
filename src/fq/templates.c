/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2018 Luca De Feo
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_factor.h"
#include "fq.h"

#ifdef T
#undef T
#endif

#define T fq
#define CAP_T FQ

#include "fq_templates/div.c"
#include "fq_templates/is_invertible.c"
#include "fq_templates/is_invertible_f.c"
#include "fq_templates/is_square.c"
#include "fq_templates/multiplicative_order.c"
#include "fq_templates/sqrt.c"

#undef CAP_T
#undef T
