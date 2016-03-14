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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

void _nmod_vec_add(mp_ptr res, mp_srcptr vec1, 
                   mp_srcptr vec2, slong len, nmod_t mod)
{
    slong i;

    if (mod.norm)
    {
        for (i = 0 ; i < len; i++)
        res[i] = _nmod_add(vec1[i], vec2[i], mod);
    } else
    {
        for (i = 0 ; i < len; i++)
        res[i] = nmod_add(vec1[i], vec2[i], mod);
    }
}
