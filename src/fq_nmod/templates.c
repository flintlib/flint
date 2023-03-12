/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2020 William Hart
    Copyright (C) 2018 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "fq_nmod.h"

#ifdef T
#undef T
#endif

#define T fq_nmod
#define CAP_T FQ_NMOD

#include "fq_templates/div.c"
#include "fq_templates/is_invertible.c"
#include "fq_templates/is_invertible_f.c"
#include "fq_templates/is_square.c"
#include "fq_templates/multiplicative_order.c"
#include "fq_templates/sqrt.c"

#undef CAP_T
#undef T
