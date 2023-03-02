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

#ifdef CALL_BABAI
#undef CALL_BABAI
#endif

#ifdef TYPE
#undef TYPE
#endif

#define FUNC_HEAD int fmpz_lll_d_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
#define CALL_BABAI(NFF, BO, HF)                                        \
do {                                                                   \
    if (NFF < 50)                                                      \
    {                                                                  \
        BO =                                                           \
            fmpz_lll_check_babai(kappa, B, U, mu, r, s, appB, expo, A, \
                                 alpha[kappa], zeros, kappamax,        \
                                 FLINT_MIN(kappamax + 1 + shift,       \
                                           n), fl);                    \
    }                                                                  \
    else                                                               \
    {                                                                  \
        BO = -1;                                                       \
    }                                                                  \
    if (BO == -1)                                                      \
    {                                                                  \
        NFF++;                                                         \
        HF =                                                           \
            fmpz_lll_check_babai_heuristic_d(kappa, B, U, mu, r, s,    \
                                             appB, expo, A,            \
                                             alpha[kappa], zeros,      \
                                             kappamax,                 \
                                             FLINT_MIN(kappamax +      \
                                                       1 + shift,      \
                                                       n), fl);        \
    }                                                                  \
} while (0)
#define TYPE 1                  /* indicates removals are desired */
#include "d_lll.c"
#undef FUNC_HEAD
#undef CALL_BABAI
#undef TYPE
