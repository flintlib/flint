/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

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

#define FUNC_HEAD int fmpz_lll_check_babai(int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double *s, \
       d_mat_t appB, int *expo, fmpz_gram_t A, \
       int a, int zeros, int kappamax, int n, const fmpz_lll_t fl)
#define LIMIT kappa
#define COMPUTE(G, I, J, C)                                         \
do {                                                                \
    if (I != J)                                                     \
        d_mat_entry(G, I, J) =                                      \
                    _d_vec_dot(appB->rows[I], appB->rows[J], C);    \
    else                                                            \
        d_mat_entry(G, I, J) = _d_vec_norm(appB->rows[I], C);       \
} while (0)
#define TYPE 1
#include "babai.c"
#undef FUNC_HEAD
#undef LIMIT
#undef COMPUTE
#undef TYPE
