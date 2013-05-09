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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_factor.h"
#include "ulong_extras.h"

void
_fmpz_factor_extend_factor_ui(fmpz_factor_t factor, mp_limb_t n)
{
    len_t i, len;
    n_factor_t nfac;

    if (n == 0)
    {
        _fmpz_factor_set_length(factor, 0);
        factor->sign = 0;
        return;
    }

    n_factor_init(&nfac);
    n_factor(&nfac, n, 0);

    len = factor->num;

    _fmpz_factor_fit_length(factor, len + nfac.num);
    _fmpz_factor_set_length(factor, len + nfac.num);

    for (i = 0; i < nfac.num; i++)
    {
        fmpz_set_ui(factor->p + len + i, nfac.p[i]);
        factor->exp[len + i] = nfac.exp[i];
    }
}
