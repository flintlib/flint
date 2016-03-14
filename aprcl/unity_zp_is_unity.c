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

slong
unity_zp_is_unity(unity_zp f)
{
    ulong result;
    ulong i, p_pow;
    unity_zp unity;

    p_pow = n_pow(f->p, f->exp);
    unity_zp_init(unity, f->p, f->exp, f->n);

    /* if the power was not found returns -1 */
    result = -1;
    for (i = 0; i < p_pow; i++)
    {
        /* set unity = \zeta_{p^k}^i */
        unity_zp_set_zero(unity);
        unity_zp_coeff_set_ui(unity, i, 1);

        /* check if f = zeta_{p^k}^i */
        if (unity_zp_equal(unity, f) == 1)
        {
            /* if so, returns \zeta_{p^k} power */
            result = i;
            break;
        }
    }

    unity_zp_clear(unity);
    return result;
}

