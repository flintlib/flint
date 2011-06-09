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

#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

mp_limb_t
_nmod_vec_dot(mp_srcptr vec1, mp_srcptr vec2, long len, nmod_t mod, int nlimbs)
{
    mp_limb_t t0, t1;
    mp_limb_t s0, s1, s2;
    long i;

    s0 = s1 = s2 = 0UL;

    switch (nlimbs)
    {
        case 1:
            for (i = 0; i < len; i++)
            {
                s0 += vec1[i] * vec2[i];
            }
            NMOD_RED(s0, s0, mod);
            break;

        case 2:
            for (i = 0; i < len; i++)
            {
                umul_ppmm(t1, t0, vec1[i], vec2[i]);
                add_ssaaaa(s1, s0, s1, s0, t1, t0);
            }
            NMOD2_RED2(s0, s1, s0, mod);
            break;

        default:
            for (i = 0; i < len; i++)
            {
                umul_ppmm(t1, t0, vec1[i], vec2[i]);
                add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0);
            }
            NMOD_RED(s2, s2, mod);
            NMOD_RED3(s0, s2, s1, s0, mod);
            break;
    }

    return s0;
}
