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

/* if returns >= 0 and h % p != 0 lambda_p = 1 */
slong
_is_prime_jacobi_check_pk(const unity_zp j, const fmpz_t u, ulong v)
{
    slong h;
    ulong i, p_pow;
    unity_zp j1, j2, temp, aut;

    p_pow = n_pow(j->p, j->exp);
    unity_zp_init(j1, j->p, j->exp, j->n);
    unity_zp_init(j2, j->p, j->exp, j->n);
    unity_zp_init(temp, j->p, j->exp, j->n);
    unity_zp_init(aut, j->p, j->exp, j->n);

    unity_zp_coeff_set_ui(j1, 0, 1);
    unity_zp_coeff_set_ui(j2, 0, 1);

    for (i = 1; i <= p_pow; i++)
    {
        if (i % j->p == 0) continue;

        unity_zp_pow_ui(temp, j, i);
        _unity_zp_reduce_cyclotomic(temp);
        unity_zp_aut_inv(aut, temp, i);
        unity_zp_mul(temp, j1, aut);
        unity_zp_swap(temp, j1);
    
        unity_zp_pow_ui(temp, j, (v * i) / p_pow);
        _unity_zp_reduce_cyclotomic(temp);
        unity_zp_aut_inv(aut, temp, i);
        unity_zp_mul(temp, j2, aut);
        unity_zp_swap(temp, j2);
    }

    unity_zp_pow_fmpz(temp, j1, u);
    unity_zp_mul(j1, j2, temp);
    
    h = unity_zp_is_unity(j1);

    unity_zp_clear(aut);
    unity_zp_clear(j1);
    unity_zp_clear(j2);
    unity_zp_clear(temp);

    return h;
}

/* if returns 1 and n = 1 mod 3 then lambda_2 = 1*/
int
_is_prime_jacobi_check_21(ulong q, const fmpz_t n)
{
    int result;
    fmpz_t qpow, ncmp, temp;

    fmpz_init(temp);
    fmpz_init_set_ui(qpow, q);
    fmpz_init_set(ncmp, n);

    fmpz_sub_ui(ncmp, ncmp, 1);
    fmpz_fdiv_q_2exp(temp, ncmp, 1);
    fmpz_powm(qpow, qpow, temp, n);

    result = 0;
    if (fmpz_equal_ui(qpow, 1) || fmpz_equal(qpow, ncmp))
        result = 1;

    fmpz_clear(temp);
    fmpz_clear(qpow);
    fmpz_clear(ncmp);

    return result;
}

slong
_is_prime_jacobi_check_22(const unity_zp j, const fmpz_t u, ulong v, ulong q)
{
    slong h;
    unity_zp j1, j2, j_pow;

    unity_zp_init(j_pow, 2, 2, j->n);
    unity_zp_init(j1, 2, 2, j->n);
    unity_zp_init(j2, 2, 2, j->n);

    unity_zp_mul(j_pow, j, j);
    unity_zp_mul_scalar_ui(j1, j_pow, q);

    if (v == 1)
        unity_zp_coeff_set_ui(j2, 0, 1);
    else if (v == 3)
        unity_zp_swap(j2, j_pow);

    unity_zp_pow_fmpz(j_pow, j1, u);
    unity_zp_mul(j1, j2, j_pow);

    h = unity_zp_is_unity(j1);

    unity_zp_clear(j_pow);
    unity_zp_clear(j1);
    unity_zp_clear(j2);

    return h;
}

slong
_is_prime_jacobi_check_2k(const unity_zp j, const unity_zp j2_1
        , const unity_zp j2_2, const fmpz_t u, ulong v)
{
    slong h;
    ulong i, p_pow;
    unity_zp j_j2, j1, j2, temp, aut;

    p_pow = n_pow(j->p, j->exp);
    unity_zp_init(temp, 2, j->exp, j->n);
    unity_zp_init(j_j2, 2, j->exp, j->n);
    unity_zp_init(aut, 2, j->exp, j->n);
    unity_zp_init(j1, 2, j->exp, j->n);
    unity_zp_init(j2, 2, j->exp, j->n);

    unity_zp_coeff_set_ui(j1, 0, 1);
    unity_zp_coeff_set_ui(j2, 0, 1);

    unity_zp_mul(j_j2, j, j2_1);

    for (i = 1; i < p_pow;)
    {
        unity_zp_pow_ui(temp, j_j2, i);
        _unity_zp_reduce_cyclotomic(temp);
        unity_zp_aut_inv(aut, temp, i);
        unity_zp_mul(temp, j1, aut);
        unity_zp_swap(temp, j1);

        unity_zp_pow_ui(temp, j_j2, (v * i) / p_pow);
        _unity_zp_reduce_cyclotomic(temp);
        unity_zp_aut_inv(aut, temp, i);
        unity_zp_mul(temp, j2, aut);
        unity_zp_swap(temp, j2);

        i += 2;

        unity_zp_pow_ui(temp, j_j2, i);
        _unity_zp_reduce_cyclotomic(temp);
        unity_zp_aut_inv(aut, temp, i);
        unity_zp_mul(temp, j1, aut);
        unity_zp_swap(temp, j1);

        unity_zp_pow_ui(temp, j_j2, (v * i) / p_pow);
        _unity_zp_reduce_cyclotomic(temp);
        unity_zp_aut_inv(aut, temp, i);
        unity_zp_mul(temp, j2, aut);
        unity_zp_swap(temp, j2);

        i += 6;
    }

    if (v % 8 != 1 && v % 8 != 3)
    {
        unity_zp_mul(temp, j2_2, j2_2);
        unity_zp_mul(j_j2, j2, temp);
        unity_zp_swap(j_j2, j2);
    }

    unity_zp_pow_fmpz(temp, j1, u);
    unity_zp_mul(j1, j2, temp);
    
    h = unity_zp_is_unity(j1);

    unity_zp_clear(aut);
    unity_zp_clear(j1);
    unity_zp_clear(j2);
    unity_zp_clear(j_j2);
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

