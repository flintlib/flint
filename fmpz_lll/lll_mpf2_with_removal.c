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

#ifdef USE_NEWD
#undef USE_NEWD
#endif

#define FUNC_HEAD int fmpz_lll_mpf2_with_removal(fmpz_mat_t B, mp_bitcnt_t prec, const fmpz_t gs_B, const fmpz_lll_t fl)
#define USE_NEWD(ND, FLAG, GSN)                                        \
do {                                                                   \
    ND = B->r;                                                         \
    fmpz_init(GSN);                                                    \
    for (i = d - 1; (i >= 0) && (FLAG > 0); i--)                       \
    {                                                                  \
        fmpz_set_mpf(GSN, mpf_mat_entry(r, i, i));                     \
        if ((FLAG = fmpz_cmp(GSN, gs_B)) > 0)                          \
        {                                                              \
            ND--;                                                      \
        }                                                              \
    }                                                                  \
    fmpz_clear(GSN);                                                   \
} while (0)
#include "mpf2_lll.c"
#undef FUNC_HEAD
#undef USE_NEWD
