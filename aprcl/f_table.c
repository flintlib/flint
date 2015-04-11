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

    Copyright (C) 2015 Vladimir Glazachev
   
******************************************************************************/

#include "aprcl.h"
#include "ulong_extras.h"

mp_ptr f_table(const ulong q)
{
    int i;
    ulong g, g_pow, g_comp, align;
    mp_ptr table;

    g = n_primitive_root_prime(q);
    table = _nmod_vec_init(q - 2);
  
    g_pow = g;
    for (i = 0; i < q - 2; i++)
    {
        g_comp = 1 - g_pow;
        align = UWORD_MAX - g_comp;
        align /= q;
        align += 1;
        g_comp += align * q;
        table[i] = n_discrete_log_bsgs(g_comp, g, q);
        g_pow *= g;
    }
    return table;
}

