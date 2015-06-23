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

void 
_jacobi_pq_general(unity_zp f, const mp_ptr table
        , ulong p, ulong q, ulong k, ulong a, ulong b)
{
    int i, j;
    ulong size, pow, pow_dec;

    unity_zp_set_zero(f);

    pow_dec = n_pow(p, k - 1);
    size = (p - 1) * pow_dec;
    pow = pow_dec * p;

    for (i = 1; i < q - 1; i++)
    {
        ulong l = (a * (i) + b * table[i]) % pow;
        if (l < size) 
            unity_zp_coeff_inc(f, l); 
        else
        {
            for (j = 0; j < p - 1; j++)
            {
                l -= pow_dec;
                unity_zp_coeff_dec(f, l);
            }
        }
    }
}

void jacobi_pq(unity_zp f, ulong q, ulong p)
{
    ulong k;
    mp_ptr table;

    table = f_table(q);
    k = p_power_in_q(q - 1, p);

    _jacobi_pq_general(f, table, p, q, k, 1, 1);

    _nmod_vec_clear(table);
}

void jacobi_2q_one(unity_zp f, ulong q)
{
    ulong k;
    mp_ptr table;

    table = f_table(q);
    k = p_power_in_q(q - 1, 2);

    _jacobi_pq_general(f, table, 2, q, k, 2, 1);

    _nmod_vec_clear(table);
}

void jacobi_2q_two(unity_zp f, ulong q)
{
    ulong a, b, k;
    mp_ptr table;

    table = f_table(q);
    k = p_power_in_q(q - 1, 2);
    b = n_pow(2, k - 3);
    a = 3 * b;

    _jacobi_pq_general(f, table, 2, q, k, a, b);

    _nmod_vec_clear(table);
}

