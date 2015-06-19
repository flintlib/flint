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

int
_is_prime_jacobi_additional_test(const fmpz_t n, ulong p)
{
    int result, p_counter, counter;
    ulong q, k;
    fmpz_t npow, qmod;

    result = 0;
    p_counter = 50;
    counter = 9;

    fmpz_init(npow);
    fmpz_init(qmod);

    while (p_counter > 0)
    {
        q = 2 * counter * p + 1;
        if (n_is_prime(q))
        {
            fmpz_set_ui(qmod, q);
            fmpz_powm_ui(npow, n, (q - 1) / p, qmod);
            if (fmpz_equal_ui(npow, 1) == 0)
                break;
            p_counter--;
        }
        counter += 2;
    }

    k = p_power_in_q(q - 1, p);
    if (p_counter != 0)
    {
        ulong v, h;
        fmpz_t u;
        unity_zp jacobi_sum;

        fmpz_init(u);
        unity_zp_init(jacobi_sum, p, k, n);

        jacobi_pq(jacobi_sum, q, p);
        fmpz_tdiv_q_ui(u, n, n_pow(p, k));
        v = fmpz_tdiv_ui(n, n_pow(p, k));

        if (p == 2)
        {
            h = _is_prime_jacobi_check_22(jacobi_sum, u, v, q);
            if (h < 0 || h % 2 == 0)
                result = 2;
            else
            {
                fmpz_t ndec, ndecdiv, qpow;
                fmpz_init_set(ndec, n);
                fmpz_init(ndecdiv);
                fmpz_init_set_ui(qpow, q);
                fmpz_sub_ui(ndec, ndec, 1);
                fmpz_fdiv_q_2exp(ndecdiv, ndec, 1);
                fmpz_powm(qpow, qpow, ndecdiv, n);

                if (fmpz_equal(qpow, ndec) == 1)
                    result = 1;
                else
                    result = 2;

                fmpz_clear(ndec);
                fmpz_clear(ndecdiv);
                fmpz_clear(qpow);
            }
        }
        else
        {
            h = _is_prime_jacobi_check_pk(jacobi_sum, u, v);
            if (h < 0 || h % p == 0)
                result = 2;
            else
                result = 1;
        }

        fmpz_clear(u);
        unity_zp_clear(jacobi_sum);
    }

    fmpz_clear(npow);
    fmpz_clear(qmod);
    return result;
}

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
slong
_is_prime_jacobi_check_21(ulong q, const fmpz_t n)
{
    slong h;
    fmpz_t qpow, ndec, temp;

    fmpz_init(temp);
    fmpz_init_set_ui(qpow, q);
    fmpz_init_set(ndec, n);

    fmpz_sub_ui(ndec, ndec, 1);
    fmpz_fdiv_q_2exp(temp, ndec, 1);
    fmpz_powm(qpow, qpow, temp, n);

    h = -1;
    if (fmpz_equal_ui(qpow, 1))
        h = 0;
    if (fmpz_equal(qpow, ndec))
        h = 1;

    fmpz_clear(temp);
    fmpz_clear(qpow);
    fmpz_clear(ndec);

    return h;
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
    fmpz_t ndec, ndecdiv;

    fmpz_init(temp);
    fmpz_init(p2);
    fmpz_init(ndecdiv);
    fmpz_init_set(ndec, n);
    fmpz_sub_ui(ndec, ndec, 1);
    fmpz_fdiv_q_2exp(ndecdiv, ndec, 1);
    
    result = PROBABPRIME;
    lambdas = (int*) malloc(sizeof(int) * config->rs.num);

    /* nmod4 = n % 4 */
    nmod4 = fmpz_tdiv_ui(n, 4);

    /* 1) l_p = 1 if p >= 3 and n^{p-1} != 1 mod p^2; l_p = 0 otherwise. */
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

    /* check that s*R and n are coprime */
    if (is_mul_coprime_ui_fmpz(config->R, config->s, n) == 0)
        result = COMPOSITE;

    for (i = 0; i < config->qs->num; i++)
    {
        n_factor_t q_factors;
        ulong q;
        if (result == COMPOSITE)
            break;

        q = fmpz_get_ui(config->qs->p + i);

        if (fmpz_equal_ui(n, q))
        {
            result = PRIME;
            break;
        }

        /* find prime factors of q - 1 */
        n_factor_init(&q_factors);
        n_factor(&q_factors, q - 1, 1);

        for (j = 0; j < q_factors.num; j++)
        {
            int state, pind;
            ulong v, p, exp, p_pow, k;
            unity_zp jacobi_sum, jacobi_sum2_1, jacobi_sum2_2;
            fmpz_t u, q_pow;
            if (result == COMPOSITE) break;

            p = q_factors.p[j];
            exp = q_factors.exp[j];
            p_pow = n_pow(p, exp);
            pind = _p_ind(config, p);
            state = lambdas[pind];

            fmpz_init_set_ui(q_pow, q);
            if (state == 0)
            {
                fmpz_powm(q_pow, q_pow, ndec, n);
            }

            fmpz_init(u);
            fmpz_tdiv_q_ui(u, n, p_pow);
            v = fmpz_tdiv_ui(n, p_pow);

            unity_zp_init(jacobi_sum, p, exp, n);
            unity_zp_init(jacobi_sum2_1, p, exp, n);
            unity_zp_init(jacobi_sum2_2, p, exp, n);

            jacobi_pq(jacobi_sum, q, p);
            if (p == 2 && exp >= 3)
            {
                jacobi_2q_one(jacobi_sum2_1, q);
                jacobi_2q_two(jacobi_sum2_2, q);
            }

            k = p_power_in_q(q - 1, p);

            slong h;
            if (p == 2 && k == 1)
            {
                h = _is_prime_jacobi_check_21(q, n);

                if (h < 0)
                {
                    result = COMPOSITE;
                }

                if (state == 0 && h == 1 && nmod4 == 1)
                {
                    state = 1;
                    lambdas[pind] = state;
                }
                continue;
            }

            if (p == 2 && k == 2)
            {
                h = _is_prime_jacobi_check_22(jacobi_sum, u, v, q);

                if (h < 0)
                {
                    result = COMPOSITE;
                }
                if (h % 2 != 0 && state == 0 && fmpz_equal(q_pow, ndec))
                {
                    state = 1;
                    lambdas[pind] = state;
                }

                continue;
            }

            if (p == 2 && k >= 3)
            {
                h = _is_prime_jacobi_check_2k(jacobi_sum
                    , jacobi_sum2_1, jacobi_sum2_2, u, v);

                if (h < 0)
                {
                    result = COMPOSITE;
                }

                if (h % 2 != 0 && state == 0 && fmpz_equal(q_pow, ndec))
                {
                    state = 1;
                    lambdas[pind] = state;
                }

                continue;
            }

            if (p != 2)
            {
                h = _is_prime_jacobi_check_pk(jacobi_sum, u, v);

                if (h < 0)
                {
                    result = COMPOSITE;
                }

                if (h % p != 0 && state == 0)
                {
                    state = 1;
                    lambdas[pind] = state;
                }
            }

            fmpz_clear(u);
            fmpz_clear(q_pow);
            unity_zp_clear(jacobi_sum);
            unity_zp_clear(jacobi_sum2_1);
            unity_zp_clear(jacobi_sum2_2);
        }

    }

    if (result == PROBABPRIME)
    {
        for (i = 0; i < config->rs.num; i++)
        {
            if (lambdas[i] == 0)
            {
                int r = _is_prime_jacobi_additional_test(n, config->rs.p[i]);
                if (r == 2)
                    result = COMPOSITE;
                else if (r == 1)
                    lambdas[i] = 1;
                else
                    result = UNKNOWN;
            }
        }
    }

    if (result == PROBABPRIME)
        if (is_prime_final_division(n, config->s, config->R) == 1)
            result = PRIME;

    free(lambdas);
    fmpz_clear(p2);
    fmpz_clear(ndec);
    fmpz_clear(ndecdiv);
    fmpz_clear(temp);
    return result;
}

int
is_prime_jacobi(const fmpz_t n)
{
    primality_test_status result;
    aprcl_config config;
    aprcl_config_init_min_R(config, n, 2);

    result = _is_prime_jacobi(n, config);
    aprcl_config_clear(config);

    if (result == PRIME)
        return 1;
    return 0;
}

