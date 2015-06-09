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
    Returns 1 if gcd(q * r, n) == 1 or gcd(q * r, n) == n
*/
int _is_coprime(ulong q, ulong r, const fmpz_t n)
{
    int result;
    fmpz_t qr, gcd;

    fmpz_init(gcd);
    fmpz_init_set_ui(qr, q);

    fmpz_mul_ui(qr, qr, r);
    fmpz_gcd(gcd, n, qr);

    result = 1;
    if (fmpz_equal_ui(gcd, 1) == 0 && fmpz_equal(gcd, n) == 0)
        result = 0;

    fmpz_clear(qr);
    fmpz_clear(gcd);

    return result;
}

/*
    Returns the index of divisor p on R factors list.
*/
int _p_ind(const aprcl_config conf, ulong p)
{
    int i;
    for (i = 0; i < conf->rs.num; i++)
        if (p == conf->rs.p[i])
            return i;
    return -1;
}

/*
    Returns 1 if \tau^{\sigma_n-n}(\chi)=-1; otherwise returns 0.
    It is equal to check: 
        (\chi(-1) * q)^((n - 1) / 2) congruent -1 mod n
*/
int _is_gausspower_2q_equal_first(ulong q, const fmpz_t n)
{
    int result;
    fmpz_t npow, nval, ncmp;

    fmpz_init_set(npow, n);
    fmpz_init_set_ui(nval, q);
    fmpz_init_set(ncmp, n);

    /* ncmp = -1 mod n */
    fmpz_sub_ui(ncmp, ncmp, 1); 

    /* nval = (\chi(-1) * q) = ((-1)^((q - 1) / 2) * q) */
    if ((q - 1) % 2 == 1)       
    {
        fmpz_neg(nval, nval);
        fmpz_add(nval, nval, n);
    }

    /* npow = (n - 1) / 2 */
    fmpz_sub_ui(npow, npow, 1);
    fmpz_fdiv_q_2exp(npow, npow, 1); 

    /* nval = (\chi(-1) * q)^((n - 1) / 2) mod n */
    fmpz_powm(nval, nval, npow, n);
    
    result = 0;
    if (fmpz_equal(nval, ncmp))
        result = 1;

    fmpz_clear(npow);
    fmpz_clear(nval);
    fmpz_clear(ncmp);

    return result;
}

/*
    Returns 1 if \tau^{\sigma_n-n}(\chi^{p / 2}) = -1; otherwise returns 0.
    It is equal to check:
        q^((n - 1) / 2) congruent -1 mod n
*/
int _is_gausspower_2q_equal_second(ulong q, const fmpz_t n)
{
    int result;
    fmpz_t npow, nval, ncmp;

    fmpz_init_set(npow, n);
    fmpz_init_set_ui(nval, q);
    fmpz_init_set(ncmp, n);

    /* ncmp = -1 mod n */
    fmpz_sub_ui(ncmp, ncmp, 1);

    /* nval = q^((n - 1) / 2) mod n */
    fmpz_sub_ui(npow, npow, 1);
    fmpz_fdiv_q_2exp(npow, npow, 1);
    fmpz_powm(nval, nval, npow, n);    
    
    result = 0;
    if (fmpz_equal(nval, ncmp))
        result = 1;

    fmpz_clear(npow);
    fmpz_clear(nval);
    fmpz_clear(ncmp);

    return result;
}

/*  
    Returns non-negative value if \tau(\chi)^(\sigma_n - n) from <\zeta_p> 
    cyclic group. It is equal to check that for some i from 0 to p - 1 
    \tau(\chi^n) == \zeta_p^i * \tau^n(\chi). If such i exists returns i.
    Otherwise returns -1.
*/
int _is_gausspower_from_unity_p(ulong q, ulong r, const fmpz_t n)
{
    int result;
    ulong i;
    unity_zpq temp, gauss, gausspow, gausssigma;

    unity_zpq_init(gauss, q, r, n);
    unity_zpq_init(gausssigma, q, r, n);
    unity_zpq_init(gausspow, q, r, n);
    unity_zpq_init(temp, q, r, n);

    /* gauss = \tau(\chi) */
    unity_zpq_gauss_sum(gauss, q, r); 
    /* gausssigma = \tau(\chi^n) */
    unity_zpq_gauss_sum_sigma_pow(gausssigma, q, r);
    /* gausspow = \tau^n(\chi) */
    unity_zpq_pow(gausspow, gauss, n);

    result = -1;
    for (i = 0; i < r; i++)
    {
        /* temp = \zeta_p^i * \tau^n(\chi) */
        unity_zpq_mul_unity_p_pow(temp, gausspow, i);
        if (unity_zpq_equal(gausssigma, temp))
        {
            result = i;
            break;
        }
    }

    unity_zpq_clear(gauss);
    unity_zpq_clear(gausssigma);
    unity_zpq_clear(gausspow);
    unity_zpq_clear(temp);

    return result;
}

/*  
    Returns 0 if for some a = n^k mod s, where k from 1 to r - 1 we have that 
    a | n; otherwise returns 1.
*/
int _is_prime_final_division(const fmpz_t n, const fmpz_t s, ulong r)
{
    int result;
    ulong i;
    fmpz_t npow, nmul, rem;

    fmpz_init(rem);
    fmpz_init_set(npow, n);
    fmpz_mod(npow, npow, s); /* npow = n mod s */
    fmpz_init_set(nmul, npow);

    result = 1;
    for (i = 1; i < r; i++)
    {
        fmpz_mod(rem, n, npow);
        /* if npow | n */
        if (fmpz_is_zero(rem))
        {
            /* if npow != n and npow != 1 */
            if ((fmpz_equal(n, npow) == 0) && (fmpz_equal_ui(npow, 1) == 0))
            {
                /* npow | n, so n is composite */
                result = 0;
                break;
            }
        }

        /* npow = n^i mod s */
        fmpz_mul(npow, npow, nmul);
        fmpz_mod(npow, npow, s);
    }

    fmpz_clear(npow);
    fmpz_clear(nmul);
    fmpz_clear(rem);

    return result;
}

int is_prime_gauss(const fmpz_t n)
{
    int result;
    ulong i, j, k;
    ulong nmod4;
    aprcl_config conf;
    aprcl_config_init(conf, n);

    int *lambdas = (int *)malloc(sizeof(int) * conf->rs.num);
    for (i = 0; i < conf->rs.num; i++)
        lambdas[i] = 0;

    result = 1;

    /* nmod4 = N % 4 */
    nmod4 = fmpz_tdiv_ui(n, 4);

    for (i = 0; i < conf->qs->num; i++)
    {
        if (result == 0) break;

        n_factor_t q_factors;
        ulong q = fmpz_get_ui(conf->qs->p + i);

        if (fmpz_equal_ui(n, q))
        {
            result = 2;
            break;
        }

        n_factor_init(&q_factors);
        n_factor(&q_factors, q - 1, 1);

        for (j = 0; j < q_factors.num; j++)
        {
            if (result == 0) break;

            int state, pind;
            ulong p = q_factors.p[j];

            pind = _p_ind(conf, p);
            state = lambdas[pind];

            if (p == 2 && state == 0 && nmod4 == 1)
            {
                if (_is_gausspower_2q_equal_first(q, n) == 1)
                {
                    state = 3;
                    lambdas[pind] = state;
                }
            }

            if (p == 2 && (state == 0 || state == 2) && nmod4 == 3)
            {
                if (_is_gausspower_2q_equal_second(q, n) == 1)
                {
                    if (state == 2)
                        state = 3;
                    else
                        state = 1;
                    lambdas[pind] = state;
                }
            }

            for (k = 1; k <= q_factors.exp[j]; k++)
            {
                int unity_power;
                ulong r;

                r = n_pow(p, k);

                if (_is_coprime(q, r, n) == 0)
                {
                    result = 0; break;
                }

                unity_power = _is_gausspower_from_unity_p(q, r, n);
                if (unity_power < 0)
                {
                    result = 0; break;
                }

                if (p > 2 && state == 0 && unity_power > 0)
                {
                    ulong upow = unity_power;
                    if (n_gcd(r, upow) == 1)
                    {
                        state = 3;
                        lambdas[pind] = state;
                    }
                }

                if (p == 2 && unity_power > 0 && (state == 0 || state == 1) && nmod4 == 3)
                {
                    ulong upow = unity_power;
                    if (n_gcd(r, upow) == 1)
                    {
                        if (state == 0) { state = 2; lambdas[pind] = state; }
                        if (state == 1) { state = 3; lambdas[pind] = state; }
                    }
                }
            }
        }
    }

    if (result != 2)
        for (i = 0; i < conf->rs.num; i++)
            if (lambdas[i] != 3) result = 0;

    free(lambdas);

    if (result == 1)
        result = _is_prime_final_division(n, conf->s, conf->R);
    aprcl_config_clear(conf);
    if (result == 2)
        return 1;

    return result;
}

