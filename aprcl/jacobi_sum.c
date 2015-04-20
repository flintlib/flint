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

void jacobi_pq(unity_root result, ulong q, ulong p)
{
    int i;
    mp_ptr table;

    unity_init(result, n_pow(p, (q - 1) / p));
    table = f_table(q);

    for (i = 0; i < q - 2; i++)
    {
        unity_root temp;
        unity_init(temp, result->power);
        unity_nth_root(temp, i + table[i] + 1);
        unity_roots_add(result, result, temp);
        unity_clear(temp);
    }
}

void jacobi_pq_general(unity_root result, const mp_ptr table, ulong p, ulong q, ulong k, ulong a, ulong b)
{
    int i, j;
    ulong size, pow, pow_dec;

    pow_dec = n_pow(p, k - 1);
    size = (p - 1) * pow_dec;
    pow = pow_dec * p;

    unity_init(result, size);

    for (i = 0; i < q - 2; i++)
    {
        ulong l = (a * (i + 1) + b * table[i]) % pow;
        if (l < size) 
        {
            unity_nth_root_inc(result, l); 
        } else 
        {
            for (j = 0; j < p - 1; j++)
            {
                l -= pow_dec;
                unity_nth_root_dec(result, l);
            }
        }
    }
}

void jacobi_pq_not2(unity_root result, ulong q, ulong p)
{
    ulong k;
    mp_ptr table;

    table = f_table(q);
    k = q / p;

    jacobi_pq_general(result, table, p, q, k, 1, 1);
}

void jacobi_2q_one(unity_root result, ulong q)
{
    ulong k;
    mp_ptr table;

    table = f_table(q);
    k = q / 2;

    jacobi_pq_general(result, table, 2, q, k, 2, 1);
}

void jacobi_2q_two(unity_root result, ulong q)
{
    ulong a, b, k;
    mp_ptr table;

    table = f_table(q);
    k = q / 2;
    b = n_pow(2, k - 3);
    a = 3 * b;

    jacobi_pq_general(result, table, 2, q, k, a, b);
}

