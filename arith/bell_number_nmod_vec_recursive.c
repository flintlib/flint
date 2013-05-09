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
#include "arith.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

void
arith_bell_number_nmod_vec_recursive(mp_ptr b, len_t n, nmod_t mod)
{
    len_t i, k;
    mp_ptr t;

    if (n < BELL_NUMBER_TAB_SIZE)
    {
        for (i = 0; i < n; i++)
            b[i] = n_mod2_preinv(bell_number_tab[i], mod.n, mod.ninv);
        return;
    }

    n -= 1;
    t = _nmod_vec_init(n);

    t[0] = b[0] = b[1] = 1;

    for (i = 1; i < n; i++)
    {
        t[i] = t[0];
        for (k = i; k > 0; k--)
            t[k - 1] = n_addmod(t[k - 1], t[k], mod.n);

        b[i + 1] = t[0];
    }

    _nmod_vec_clear(t);
}
