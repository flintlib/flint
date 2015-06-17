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
_is_prime_jacobi_check_pk(const unity_zp j, const fmpz_t u, ulong v)
{
    slong h;
    ulong i, p_pow;
    unity_zp j1, j2, j_pow, temp, aut;

    p_pow = n_pow(j->p, j->exp);
    unity_zp_init(j_pow, j->p, j->exp, j->n);
    unity_zp_init(j1, j->p, j->exp, j->n);
    unity_zp_init(j2, j->p, j->exp, j->n);
    unity_zp_init(temp, j->p, j->exp, j->n);
    unity_zp_init(aut, j->p, j->exp, j->n);

    unity_zp_copy(j_pow, j);
    unity_zp_coeff_set_ui(j1, 0, 1);
    unity_zp_coeff_set_ui(j2, 0, 1);

    for (i = 1; i <= p_pow; i++)
    {
        if (i % j->p == 0) continue;

        unity_zp_pow_ui(j_pow, j, i);
        _unity_zp_reduce_cyclotomic(j_pow);
        unity_zp_aut_inv(aut, j_pow, i);
        unity_zp_mul(temp, j1, aut);
        unity_zp_swap(temp, j1);
    
        unity_zp_pow_ui(j_pow, j, (v * i) / p_pow);
        _unity_zp_reduce_cyclotomic(j_pow);
        unity_zp_aut_inv(aut, j_pow, i);
        unity_zp_mul(temp, j2, aut);
        unity_zp_swap(temp, j2);
    }

    unity_zp_pow_fmpz(j_pow, j1, u);
    unity_zp_mul(j1, j2, j_pow);
    
    h = unity_zp_is_unity(j1);

    unity_zp_clear(aut);
    unity_zp_clear(j_pow);
    unity_zp_clear(j1);
    unity_zp_clear(j2);
    unity_zp_clear(temp);

    return h;
}

int
_is_prime_jacobi(const fmpz_t n, const aprcl_config config)
{
    int *lambdas;
    ulong i, j, k, nmod4;
    primality_test_status result;
    fmpz_t temp;
    fmpz_t p2;

    fmpz_init(temp);
    fmpz_init(p2);
    
    result = PROBABPRIME;
    lambdas = (int*) malloc(sizeof(int) * config->rs.num);
    for (i = 0; i < config->rs.num; i++)
    {
        ulong p = config->rs.p[i];
        if (p > 2)
        {
            fmpz_set_ui(p2, p * p);
            fmpz_powm_ui(temp, n, p - 1, p2);
            if (fmpz_equal_ui(temp, 1) == 0)
                lambdas[i] = 1;
            else
                lambdas[i] = 0;
        }
        else
            lambdas[i] = 0;
    }

    
    /* 1) l_p = 1 if p >= 3 and n^{p-1} != 1 mod p^2; l_p = 0 otherwise. */


    free(lambdas);
    fmpz_clear(temp);
    return 0;
}

