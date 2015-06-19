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

/*
    Below is the implementation of primality test using Jacibi sums.
    Many steps are well described in:
        [1] "A Course in Computational Algebraic Number Theory" by H. Cohen
        [2] "Implementation of a New Primality Test" by H. Cohen and A. K. Lenstra

    The algorithm consist of N steps:
        (1.) Precomutation;
        (2.) Pseudoprime tests with Jacobi sums;
        (3.) Additionl tests;
        (4.) Final trial division and primality proving.

    The file contains the implementation of steps (2.) and (3.).
    It also contains the Jacobi sum primality test.
    A small part of implementation of step (1.) is here and most are in
    aprcl_config file.
    (4.) implemented in function is_prime_final_division().

    Standard variables:
        n - number to check for primality;
        
        R - configuration parameter;
        s - configuration parameter depends on the R; s^2 > n;

        q - prime number such that (q | s);
        p - prime number such that (p | q - 1);
        r - prime power p^k such that (r | q - 1) and not (r * p | q - 1);
        k - power of p in r;

        v - n % r;
        u - n / r;

    Standard notation:
        \zeta_p - p-th root of unity;
        \sigma_x - automorphism of Z[\zeta_p] for which
                \sigma_x(\zeta_p) = \zeta_p^x;
        \sigma_x^{-1} - inverse of \sigma_x;

        \chi_{p, q} - character defined by 
                \chi_{p, q}(g^x) = \zeta_{p^k}^x, 
                there g is a primitive root modulo q;

        J(p, q) - jacobi sum j(\chi_{p, q}, \chi_{p, q});
        J_2(q)  - jacobi sum 
                ( j(\chi_{2, q}^{2^{k - 3}}, \chi{2, q}^{3 * 2^{k - 3}}) )^2;
        J_3(q)  - jacobi sum
                j(\chi_{2, q}, \chi_{2, q}, \chi_{2, q}) =
                J(2, q) * j(\chi{2, q}^2, \chi{2, q});
*/

/*----------------------------------------------------------------------------*/

/*
    (2.a)
    Checks the case p != 2.
    Computes j0 = j_{0, p, q}, jv = j_{v, p, q} and checks that
    j0^u * jv is unity root.

    Parameters:
        j = J(p, q);
        u, v from standard variables;

    Returns:
        If there exist h such that j0^u * jv = \zeta_{p^k}^h returns h;
        otherwise returns -1.

    For details about j0 and jv see (i1a) in [2] or
    algorithm (9.1.28) step 4.a in [1].
*/
slong
_is_prime_jacobi_check_pk(const unity_zp j, const fmpz_t u, ulong v)
{
    slong h;
    ulong i, r;
    unity_zp j0, jv, temp, aut;

    /* initialization */
    r = n_pow(j->p, j->exp);            /* r = p^k */
    unity_zp_init(j0, j->p, j->exp, j->n);
    unity_zp_init(jv, j->p, j->exp, j->n);
    unity_zp_init(temp, j->p, j->exp, j->n);
    unity_zp_init(aut, j->p, j->exp, j->n);

    unity_zp_coeff_set_ui(j0, 0, 1);    /* j0 = 1 */
    unity_zp_coeff_set_ui(jv, 0, 1);    /* jv = 1 */

    /* for i in 1..p^k */
    for (i = 1; i <= r; i++)
    {
        /* only for i that coprime to p */
        if (i % j->p == 0) continue;

        /* update j0 = \prod{\sigma_i^{-1}(j^i)} */
        unity_zp_pow_ui(temp, j, i);
        _unity_zp_reduce_cyclotomic(temp);
        /* aut = \sigma_i^{-1}(temp) */
        unity_zp_aut_inv(aut, temp, i);

        /* j0 *= aut */
        unity_zp_mul(temp, j0, aut);
        unity_zp_swap(temp, j0);

        /* update jv = \prod{\sigma_i^{-1}(j^{(v * i) / r})} */
        unity_zp_pow_ui(temp, j, (v * i) / r);
        _unity_zp_reduce_cyclotomic(temp);
        /* aut = \sigma_i^{-1}(temp) */
        unity_zp_aut_inv(aut, temp, i);

        /* jv *= aut */
        unity_zp_mul(temp, jv, aut);
        unity_zp_swap(temp, jv);
    }

    /* temp = j0^u */
    unity_zp_pow_fmpz(temp, j0, u);
    /* j0 = j0^u * jv */
    unity_zp_mul(j0, jv, temp);
    
    /* try to find h */
    h = unity_zp_is_unity(j0);

    /* clear */
    unity_zp_clear(aut);
    unity_zp_clear(j0);
    unity_zp_clear(jv);
    unity_zp_clear(temp);

    return h;
}

/*
    (2.b)
    Check the case p = 2 and k = 1.

    Computes j0 = j_{0, 2, q}, jv = j_{v, 2, q} and checks that
    j0^u * jv is unity root. j^0^u * jv = (-q)^{(n - 1) / 2}.

    Parameters:
        q, n from standard variables;

    Returns:
        If j0^u * jv = 1 returns 0;
        If j0^u * jv = -1 returns 1;
        otherwise returns -1.

    For details see (i1b) in [2] or
    algorithm (9.1.28) step 4.d in [1].
*/
slong
_is_prime_jacobi_check_21(ulong q, const fmpz_t n)
{
    slong h;
    fmpz_t qpow, ndec, temp;

    /* initialization */
    fmpz_init(temp);
    fmpz_init_set_ui(qpow, q);
    fmpz_init_set(ndec, n);

    /* qpow = -q mod n */
    fmpz_sub(qpow, n, qpow);
    fmpz_sub_ui(ndec, ndec, 1);
    /* temp = (n - 1) / 2 */
    fmpz_fdiv_q_2exp(temp, ndec, 1);
    /* qpow = (-q)^{(n - 1) / 2} */
    fmpz_powm(qpow, qpow, temp, n);

    h = -1;
    /* check if qpow == +-1 */
    if (fmpz_equal_ui(qpow, 1))
        h = 0;
    if (fmpz_equal(qpow, ndec))
        h = 1;

    /* clear */
    fmpz_clear(temp);
    fmpz_clear(qpow);
    fmpz_clear(ndec);

    return h;
}

/*
    (2.c)
    Check the case p = 2 and k = 2.

    Computes j0 = j_{0, 2, q}, jv = j_{v, 2, q} and checks that
    j0^u * jv is unity root.

    Parameters:
        j = J(2, q);
        u, v, q from standard variables;

    Returns:
        If there exist h such that j0^u * jv = \zeta_{4}^h returns h
        \zeta_4^h \in (1, i, -1, -i);
        otherwise returns -1.

    For details see (i1c) in [2] or
    algorithm (9.1.28) step 4.c in [1].
*/
slong
_is_prime_jacobi_check_22(const unity_zp j, const fmpz_t u, ulong v, ulong q)
{
    slong h;
    unity_zp j0, jv, j_pow;

    /* initialization */
    unity_zp_init(j_pow, 2, 2, j->n);
    unity_zp_init(j0, 2, 2, j->n);
    unity_zp_init(jv, 2, 2, j->n);

    /* set j0 = q * j^2 */
    unity_zp_mul(j_pow, j, j);
    unity_zp_mul_scalar_ui(j0, j_pow, q);

    /*
        if v == 1 jv = 1
        if v == 3 jv = j^2
    */
    if (v == 1)
        unity_zp_coeff_set_ui(jv, 0, 1);
    else if (v == 3)
        unity_zp_swap(jv, j_pow);

    /* j0 = j0^u * jv */
    unity_zp_pow_fmpz(j_pow, j0, u);
    unity_zp_mul(j0, jv, j_pow);

    /* try to find h */
    h = unity_zp_is_unity(j0);

    /* clear */
    unity_zp_clear(j_pow);
    unity_zp_clear(j0);
    unity_zp_clear(jv);

    return h;
}

slong
_is_prime_jacobi_check_2k(const unity_zp j, const unity_zp jv_1
        , const unity_zp jv_2, const fmpz_t u, ulong v)
{
    slong h;
    ulong i, r;
    unity_zp j_jv, j0, jv, temp, aut;

    r = n_pow(j->p, j->exp);
    unity_zp_init(temp, 2, j->exp, j->n);
    unity_zp_init(j_jv, 2, j->exp, j->n);
    unity_zp_init(aut, 2, j->exp, j->n);
    unity_zp_init(j0, 2, j->exp, j->n);
    unity_zp_init(jv, 2, j->exp, j->n);

    unity_zp_coeff_set_ui(j0, 0, 1);
    unity_zp_coeff_set_ui(jv, 0, 1);

    unity_zp_mul(j_jv, j, jv_1);

    for (i = 1; i < r;)
    {
        unity_zp_pow_ui(temp, j_jv, i);
        _unity_zp_reduce_cyclotomic(temp);
        unity_zp_aut_inv(aut, temp, i);
        unity_zp_mul(temp, j0, aut);
        unity_zp_swap(temp, j0);

        unity_zp_pow_ui(temp, j_jv, (v * i) / r);
        _unity_zp_reduce_cyclotomic(temp);
        unity_zp_aut_inv(aut, temp, i);
        unity_zp_mul(temp, jv, aut);
        unity_zp_swap(temp, jv);

        i += 2;

        unity_zp_pow_ui(temp, j_jv, i);
        _unity_zp_reduce_cyclotomic(temp);
        unity_zp_aut_inv(aut, temp, i);
        unity_zp_mul(temp, j0, aut);
        unity_zp_swap(temp, j0);

        unity_zp_pow_ui(temp, j_jv, (v * i) / r);
        _unity_zp_reduce_cyclotomic(temp);
        unity_zp_aut_inv(aut, temp, i);
        unity_zp_mul(temp, jv, aut);
        unity_zp_swap(temp, jv);

        i += 6;
    }

    if (v % 8 != 1 && v % 8 != 3)
    {
        unity_zp_mul(temp, jv_2, jv_2);
        unity_zp_mul(j_jv, jv, temp);
        unity_zp_swap(j_jv, jv);
    }

    unity_zp_pow_fmpz(temp, j0, u);
    unity_zp_mul(j0, jv, temp);
    
    h = unity_zp_is_unity(j0);

    unity_zp_clear(aut);
    unity_zp_clear(j0);
    unity_zp_clear(jv);
    unity_zp_clear(j_jv);
    unity_zp_clear(temp);

    return h;
}

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

int
_is_prime_jacobi(const fmpz_t n, const aprcl_config config)
{
    int *lambdas;
    ulong i, j, nmod4;
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
            ulong v, p, exp, r, k;
            unity_zp jacobi_sum, jacobi_sum2_1, jacobi_sum2_2;
            fmpz_t u, q_pow;
            if (result == COMPOSITE) break;

            p = q_factors.p[j];
            exp = q_factors.exp[j];
            r = n_pow(p, exp);
            pind = _p_ind(config, p);
            state = lambdas[pind];

            fmpz_init_set_ui(q_pow, q);
            if (state == 0)
            {
                fmpz_powm(q_pow, q_pow, ndec, n);
            }

            fmpz_init(u);
            fmpz_tdiv_q_ui(u, n, r);
            v = fmpz_tdiv_ui(n, r);

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

