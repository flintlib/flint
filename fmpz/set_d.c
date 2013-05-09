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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

#if FLINT64   /* 2^53 */
#define DOUBLE_MAX 9007199254740992.0
#define DOUBLE_MIN -9007199254740992.0
#else
#define DOUBLE_MAX COEFF_MAX
#define DOUBLE_MIN COEFF_MIN
#endif

void
fmpz_set_d(fmpz_t f, double c)
{
    if (c >= DOUBLE_MIN && c <= DOUBLE_MAX)
    {
        _fmpz_demote(f);
        /* guaranteed to fit, since c gets truncated */
        *f = (len_t) c;
    }
    else
    {
        __mpz_struct * z = _fmpz_promote(f);
        mpz_set_d(z, c);
        _fmpz_demote_val(f);
    }
}
