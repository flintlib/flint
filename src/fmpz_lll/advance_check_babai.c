/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_lll.h"

#ifdef FUNC_HEAD
#undef FUNC_HEAD
#endif

#ifdef LIMIT
#undef LIMIT
#endif

#ifdef COMPUTE
#undef COMPUTE
#endif

#ifdef TYPE
#undef TYPE
#endif

#define FUNC_HEAD int fmpz_lll_advance_check_babai(int cur_kappa, int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double *s, \
       d_mat_t appB, int *expo, fmpz_gram_t A, \
       int a, int zeros, int kappamax, int n, const fmpz_lll_t fl)
#define LIMIT cur_kappa
#define COMPUTE(G, I, J, C)                                 \
do {                                                        \
    d_mat_entry(G, I, J) =                                  \
            _d_vec_dot(appB->rows[I], appB->rows[J], C);    \
} while (0)
#define TYPE 2
#include "babai.c"
#undef FUNC_HEAD
#undef LIMIT
#undef COMPUTE
#undef TYPE
