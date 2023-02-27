/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"
#include "ulong_extras.h"

/*
    Computes sum \zeta_{p^k}^{a * x + b * f(x)} for x = 1, 2, ..., q - 2.
*/
void 
_unity_zp_jacobi_sum_pq_general(unity_zp f, const mp_ptr table,
        ulong p, ulong q, ulong k, ulong a, ulong b)
{
    int i, j;
    ulong size, pow, pow_dec;

    unity_zp_set_zero(f);

    pow_dec = n_pow(p, k - 1);
    size = (p - 1) * pow_dec;
    pow = pow_dec * p;

    for (i = 1; i < q - 1; i++)
    {
        /* l = a * i + b * f[i] */ 
        ulong l = (a * (i) + b * table[i]) % pow;

        /* if l < (p - 1)p^{k - 1} increase f[l] by one */
        if (l < size)
        {
            unity_zp_coeff_inc(f, l);
        }
        else /* else decrease f[l - jp^{k - 1}] by one for j in [0; p - 2] */
        {
            for (j = 0; j < p - 1; j++)
            {
                l -= pow_dec;
                unity_zp_coeff_dec(f, l);
            }
        }
    }
}

/*
    Computes sum \zeta_{p^k}^{x + f(x)} for x = 1, 2, ..., q - 2.
*/
void
unity_zp_jacobi_sum_pq(unity_zp f, ulong q, ulong p)
{
    ulong k;
    mp_ptr table;

    table = aprcl_f_table(q);
    k = aprcl_p_power_in_q(q - 1, p);

    _unity_zp_jacobi_sum_pq_general(f, table, p, q, k, 1, 1);

    _nmod_vec_clear(table);
}

/*
    Computes sum \zeta_{p^k}^{2 * x + b * f(x)} for x = 1, 2, ..., q - 2.
*/
void
unity_zp_jacobi_sum_2q_one(unity_zp f, ulong q)
{
    ulong k;
    mp_ptr table;

    table = aprcl_f_table(q);
    k = aprcl_p_power_in_q(q - 1, 2);

    _unity_zp_jacobi_sum_pq_general(f, table, 2, q, k, 2, 1);

    _nmod_vec_clear(table);
}

/*
    Computes sum \zeta_{p^k}^{3 * 2^{k - 3} * x + 2^{k - 3} * f(x)}
    for x = 1, 2, ..., q - 2.
*/
void
unity_zp_jacobi_sum_2q_two(unity_zp f, ulong q)
{
    ulong a, b, k;
    mp_ptr table;

    table = aprcl_f_table(q);
    k = aprcl_p_power_in_q(q - 1, 2);
    b = n_pow(2, k - 3);    /* b = 2^{k - 3}        */
    a = 3 * b;              /* a = 3 * 2^{k - 3}    */

    _unity_zp_jacobi_sum_pq_general(f, table, 2, q, k, a, b);

    _nmod_vec_clear(table);
}

