/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2015 Elena Sergeicheva
    Copyright (C) 2018 Tommy Hofmann
    Copyright (C) 2020, 2021 William Hart
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz.h"
#include "fmpz_mod_mat.h"
#include "fq.h"
#include "fq_vec.h"
#include "fq_mat.h"
#include "fq_poly.h"

#ifdef T
#undef T
#endif

#define T fq
#define CAP_T FQ

#include "fq_mat_templates/add.c"
#include "fq_mat_templates/can_solve.c"
#include "fq_mat_templates/charpoly.c"
#include "fq_mat_templates/clear.c"
#include "fq_mat_templates/concat_horizontal.c"
#include "fq_mat_templates/concat_vertical.c"
#include "fq_mat_templates/equal.c"
#include "fq_mat_templates/init.c"
#include "fq_mat_templates/init_set.c"
#include "fq_mat_templates/inv.c"
#include "fq_mat_templates/io.c"
#include "fq_mat_templates/is_one.c"
#include "fq_mat_templates/is_zero.c"
#include "fq_mat_templates/lu.c"
#include "fq_mat_templates/lu_classical.c"
#include "fq_mat_templates/lu_recursive.c"
#include "fq_mat_templates/mat_entry_set.c"
#include "fq_mat_templates/mat_invert_cols.c"
#include "fq_mat_templates/mat_swap_cols.c"
#include "fq_mat_templates/mat_swap_entrywise.c"
#include "fq_mat_templates/minpoly.c"
#include "fq_mat_templates/mul.c"
#include "fq_mat_templates/mul_classical.c"
#include "fq_mat_templates/mul_KS.c"
#include "fq_mat_templates/mul_vec.c"
#include "fq_mat_templates/neg.c"
#include "fq_mat_templates/nullspace.c"
#include "fq_mat_templates/one.c"
#include "fq_mat_templates/randops.c"
#include "fq_mat_templates/randpermdiag.c"
#include "fq_mat_templates/randrank.c"
#include "fq_mat_templates/randtest.c"
#include "fq_mat_templates/randtril.c"
#include "fq_mat_templates/randtriu.c"
#include "fq_mat_templates/rank.c"
#include "fq_mat_templates/rref.c"
#include "fq_mat_templates/set.c"
#include "fq_mat_templates/set_fmpz_mod_mat.c"
#include "fq_mat_templates/set_nmod_mat.c"
#include "fq_mat_templates/similarity.c"
#include "fq_mat_templates/solve.c"
#include "fq_mat_templates/solve_tril.c"
#include "fq_mat_templates/solve_tril_classical.c"
#include "fq_mat_templates/solve_tril_recursive.c"
#include "fq_mat_templates/solve_triu.c"
#include "fq_mat_templates/solve_triu_classical.c"
#include "fq_mat_templates/solve_triu_recursive.c"
#include "fq_mat_templates/sub.c"
#include "fq_mat_templates/submul.c"
#include "fq_mat_templates/swap.c"
#include "fq_mat_templates/vec_mul.c"
#include "fq_mat_templates/window_clear.c"
#include "fq_mat_templates/window_init.c"
#include "fq_mat_templates/zero.c"

#undef CAP_T
#undef T
