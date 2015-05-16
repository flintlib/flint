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

void unity_automorphism_inv(unity_root_mod a, ulong x, unity_root_mod b)
{
    int i, j;

    ulong p_pow, p_pow_kk, m, p_pow_inv;
    p_pow = a->power;
    p_pow_kk = p_pow / a->p;
    m = (a->p - 1) * p_pow_kk;
    p_pow_inv = n_preinvert_limb(p_pow);

    for (i = 0; i < m; i++)
    {
        slong a_ind;
        fmpz_t a_value;

        a_ind = n_mulmod2_preinv(x, i, p_pow, p_pow_inv);
        fmpz_mod_poly_get_coeff_fmpz(a_value, a->poly, a_ind);
        fmpz_mod_poly_set_coeff_fmpz(b->poly, i, a_value);
    }

    for (i = m; i < p_pow; i++)
    {
        slong a_ind = n_mulmod2_preinv(x, i, p_pow, p_pow_inv);
        for (j = 1; j < a->p; j++)
        {
            slong b_ind;
            fmpz_t a_value, b_value;

            b_ind = i - j * p_pow_kk;
            fmpz_mod_poly_get_coeff_fmpz(a_value, a->poly, a_ind);
            fmpz_mod_poly_get_coeff_fmpz(b_value, b->poly, b_ind);
            fmpz_sub(b_value, b_value, a_value);
            fmpz_mod_poly_set_coeff_fmpz(b->poly, b_ind, b_value);
        }
    }
}

