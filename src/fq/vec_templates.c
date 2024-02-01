/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz_poly.h"
#include "fq.h"
#include "fq_vec.h"

#ifdef T
#undef T
#endif

#define T fq
#define CAP_T FQ

#include "fq_vec_templates/add.c"
#include "fq_vec_templates/clear.c"
#include "fq_vec_templates/equal.c"
#include "fq_vec_templates/init.c"
#include "fq_vec_templates/io.c"
#include "fq_vec_templates/is_zero.c"
#include "fq_vec_templates/neg.c"
#include "fq_vec_templates/randtest.c"
#include "fq_vec_templates/scalar_addmul_fq.c"
#include "fq_vec_templates/scalar_mul_fq.c"
#include "fq_vec_templates/scalar_submul_fq.c"
#include "fq_vec_templates/set.c"
#include "fq_vec_templates/sub.c"
#include "fq_vec_templates/swap.c"
#include "fq_vec_templates/zero.c"

#undef CAP_T
#undef T
