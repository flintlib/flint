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
*/
int _is_gausspower_2q_equal_first(ulong q, const fmpz_t n)
{
    int result;

    fmpz_t npow, nval, ncmp;
    fmpz_init_set_ui(nval, q);
    fmpz_init_set(ncmp, n);
    fmpz_sub_ui(ncmp, ncmp, 1);
    if ((q - 1) % 2 == 1)
    {
        fmpz_neg(nval, nval);
        fmpz_add(nval, nval, n);
    }

    fmpz_init_set(npow, n);
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
    Returns 1 if \tau^{\sigma_n-n}(\chi^{p / 2}) = -1; otherwise returns 0.
*/
int _is_gausspower_2q_equal_second(ulong q, const fmpz_t n)
{
    int result;

    fmpz_t npow, nval, ncmp;
    fmpz_init_set_ui(nval, q);
    fmpz_init_set(ncmp, n);
    fmpz_sub_ui(ncmp, ncmp, 1);

    fmpz_init_set(npow, n);
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

int is_condition_one()
{
    return 0;
}

/*  
    if for some i from 0 to p-1 \tau(\chi^n) == \tau^n(\chi)*\zeta_p^i return 1; 
    otherwise return 0 
*/
int is_condition_two(const unity_zpq f, const unity_zpq g)
{
    int result;
    ulong i;
    unity_zpq temp;
    unity_zpq_init(temp, f->q, f->p, f->n);

    result = 0;
    for (i = 0; i < f->p; i++)
    {
        unity_zpq_mul_unity_p_pow(temp, g, i);
        if (unity_zpq_equal(f, g))
        { 
            result = 1;
            break;
        }
    }

    unity_zpq_clear(temp);
    return result;
}

/*  
    if for some a = n^k mod s, where k from 1 to r - 1 we have that 
    a | n returns 0; otherwise returns 1.
*/
int is_prime_gauss_final_division(const fmpz_t n, const fmpz_t s, ulong r)
{
    int result = 0;
    ulong i;
    fmpz_t npow, nmul, rem;
    fmpz_init_set(npow, n);
    fmpz_mod(npow, npow, s);
    fmpz_init_set(nmul, npow);
    fmpz_init(rem);

    result = 1;
    for (i = 1; i < r; i++)
    {
        fmpz_mod(rem, n, npow);
        if (fmpz_is_zero(rem))
        {
            result = 0;
            break;
        }
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

    flint_printf("\nR = %wu\n", conf->R);
    flint_printf("\nNMOD4 = %we\n", nmod4);

    for (i = 0; i < conf->qs->num; i++)
    {
        n_factor_t q_factors;
        ulong q = fmpz_get_ui(conf->qs->p + i);

        n_factor_init(&q_factors);
        n_factor(&q_factors, q - 1, 1);

        for (j = 0; j < q_factors.num; j++)
        {
            int state, pind;
            ulong p = q_factors.p[j];

            flint_printf("q = %wu; p = %wu\n", q, p);

            pind = _p_ind(conf, p);
            state = lambdas[pind];
            flint_printf("state = %wu\n", state);

            if (state == 0 && nmod4 == 1)
            {   
                flint_printf("nmod4 = 1\n");
                state = _is_gausspower_2q_equal_first(q, n);
                lambdas[pind] = state;
            }

            if (state == 0 && nmod4 == 3)
            {
                flint_printf("nmod4 = 3\n");
                state = _is_gausspower_2q_equal_second(q, n);
                lambdas[pind] = state;
            }

        }


    }

    for (i = 0; i < conf->rs.num; i++)
        flint_printf("%i ", lambdas[i]);
    flint_printf("\n");

/*
    for (i = 0; i < conf->qs->num; i++)
    {
        ulong q, r;
        n_factor_t q_factors;

        q = conf->qs->p[i];
        n_factor_init(&q_factors);
        n_factor(&q_factors, q, 1);

        for (j = 0; j <= q_factors.num; j++)
        {
            if (result == 0)
                break;

            ulong l = q_factors.p[j];
            for (k = 1; k <= q_factors.exp[j]; k++)
            {
                unity_zpq gauss, gausssigma, gausspow;
                fmpz_t qr, gcd;
                ulong nmod4;

                fmpz_init(qr);
                fmpz_init(gcd);
                r = n_pow(l, k);
                fmpz_set_ui(qr, q * r);

                fmpz_gcd(gcd, n, qr);
                if (fmpz_cmp_ui(gcd, 1) == 1)
                {
                    return 0;
                }

                unity_zpq_init(gauss, q, r, n);
                unity_zpq_init(gausssigma, q, r, n);
                unity_zpq_init(gausspow, q, r, n);

                unity_zpq_gauss_sum(gauss, q, r); 
                unity_zpq_gauss_sum_sigma_pow(gausssigma, q, r);

                unity_zpq_pow(gausspow, gauss, n);

                if (q_factors.p[j] == 2)
                {
                    nmod4 = fmpz_tdiv_ui(n, 4);
                    if (nmod4 == 1)
                    {
                        /* \tau = -1 
                    }

                    if (nmod4 == 3)
                    {
                        /* \tau generator  and \tau(r/2) = -1
                    }
                }
                else
                {
                    fmpz_t npow, lsquare;
                    fmpz_init(npow);
                    fmpz_init_set_ui(lsquare, l * l);
                    fmpz_powm_ui(npow, n, l - 1, lsquare);
                    if (fmpz_cmp_ui(npow, 1))
                    {
                        result = 0;
                    }
                    fmpz_clear(npow);
                    fmpz_clear(lsquare);
                }
              

                if (result != 0)
                    if (is_condition_two(gausssigma, gausspow) == 0)
                        result = 0;

                unity_zpq_clear(gauss);
                unity_zpq_clear(gausssigma);
                unity_zpq_clear(gausspow);
                fmpz_clear(qr);
                fmpz_clear(gcd);
            }
        }
    }
    */

    free(lambdas);
    result = is_prime_gauss_final_division(n, conf->s, conf->R);
    aprcl_config_clear(conf);
    return result;
}

