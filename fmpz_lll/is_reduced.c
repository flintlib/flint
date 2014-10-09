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

int
fmpz_lll_is_reduced(const fmpz_mat_t B, const fmpz_lll_t fl, mp_bitcnt_t prec)
{
    return ((fmpz_lll_is_reduced_d(B, fl)
             || fmpz_lll_is_reduced_mpfr(B, fl, prec))
            || ((fl->rt == Z_BASIS) ?
                fmpz_mat_is_reduced(B, fl->delta,
                                    fl->eta) : fmpz_mat_is_reduced_gram(B,
                                                                        fl->delta,
                                                                        fl->eta)));
}
