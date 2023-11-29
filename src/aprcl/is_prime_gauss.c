/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "aprcl.h"

/*
    Returns 1 if \tau^{\sigma_n-n}(\chi)=-1; otherwise returns 0.
    It is equal to check:
        (\chi(-1) * q)^((n - 1) / 2) congruent -1 mod n

    \tau is the Gauss sum;
    \chi is the character of conductor q and order 2 (quadratic character);
    \tau(\chi)^\sigma_n means \sigma_n(\tau(\chi)),
        there \sigma_n is the ring automorphism
*/
int
_aprcl_is_gausspower_2q_equal_first(ulong q, const fmpz_t n)
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

    \tau is the Gauss sum;
    \chi is the character of conductor q and order 2 (quadratic character);
    \tau(\chi)^\sigma_n means \sigma_n(\tau(\chi)),
        there \sigma_n is the ring automorphism
*/
int
_aprcl_is_gausspower_2q_equal_second(ulong q, const fmpz_t n)
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

    \tau is the Gauss sum;
    \chi is the character of conductor q and order r;
    \tau(\chi)^\sigma_n means \sigma_n(\tau(\chi)),
        there \sigma_n is the ring automorphism
*/
slong
_aprcl_is_gausspower_from_unity_p(ulong q, ulong r, const fmpz_t n)
{
    slong result;
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

primality_test_status
_aprcl_is_prime_gauss(const fmpz_t n, const aprcl_config config)
{
    int *lambdas;
    ulong i, j, k, nmod4;
    primality_test_status result;

    /*
        Condition (Lp) is satisfied iff:
        For each p | R we must show that for all prime r | n and all
        positive integers a there exists l such that:
            r^(p-1) congruent n^(l(p-1)) mod p^a

        If (Lp), then lambdas_p = 3.
    */
    lambdas = (int*) flint_malloc(sizeof(int) * config->rs.num);
    for (i = 0; i < config->rs.num; i++)
        lambdas[i] = 0;

    result = PROBABPRIME;

    /* nmod4 = n % 4 */
    nmod4 = fmpz_tdiv_ui(n, 4);

    /* for every prime q | s */
    for (i = 0; i < config->qs->num; i++)
    {
        n_factor_t q_factors;
        ulong q;
        if (result == COMPOSITE) break;

        q = fmpz_get_ui(config->qs->p + i);

        /* n == q, q - prime => n - prime */
        if (fmpz_equal_ui(n, q))
        {
            result = PRIME;
            break;
        }

        /* find prime factors of q - 1 */
        n_factor_init(&q_factors);
        n_factor(&q_factors, q - 1, 1);

        /* for every prime p | q - 1 */
        for (j = 0; j < q_factors.num; j++)
        {
            int state, pind;
            ulong p;
            if (result == COMPOSITE) break;

            p = q_factors.p[j];

            pind = _aprcl_p_ind(config, p);
            state = lambdas[pind];

            /*
                (Lp.a)
                if p == 2 and n = 1 mod 4 then (Lp) is equal to:
                    for quadratic character \chi (\tau(\chi))^(\sigma_n-n) = -1
            */
            if (p == 2 && state == 0 && nmod4 == 1)
            {
                if (_aprcl_is_gausspower_2q_equal_first(q, n) == 1)
                {
                    state = 3;
                    lambdas[pind] = state;
                }
            }

            /*
                (Lp.b)
                if p == 2, r = 2^k >= 4 and n = 3 mod 4 then (Lp) is equal to:
                    1) for quadratic character \chi
                        (\tau(\chi^(r / 2)))^(\sigma_n-n) = -1
                    2) for character \chi = \chi_{r, q}
                        (\tau(\chi))^(\sigma_n-n) is a generator of cyclic
                        group <\zeta_r>

                if 1) is true, then lambdas_p = 1
                if 2) is true, then lambdas_p = 2
                if 1) and 2) is true, then lambdas_p = 3
            */
            if (p == 2 && (state == 0 || state == 2) && nmod4 == 3)
            {
                if (_aprcl_is_gausspower_2q_equal_second(q, n) == 1)
                {
                    if (state == 2)
                        state = 3;
                    else
                        state = 1;
                    lambdas[pind] = state;
                }
            }

            /* for every prime power p^k | q - 1 */
            for (k = 1; k <= q_factors.exp[j]; k++)
            {
                int unity_power;
                ulong r;

                /* r = p^k */
                r = n_pow(p, k);

                /* if gcd(q*r, n) != 1 */
                if (aprcl_is_mul_coprime_ui_ui(q, r, n) == 0)
                {
                    result = COMPOSITE;
                    break;
                }

                /*
                    if exists z such that \tau(\chi^n) = \zeta_r^z*\tau^n(\chi)
                    unity_power = z; otherwise unity_power = -1
                */
                unity_power = _aprcl_is_gausspower_from_unity_p(q, r, n);

                /* if unity_power < 0 then n is composite */
                if (unity_power < 0)
                {
                    result = COMPOSITE;
                    break;
                }

                /*
                    (Lp.c)
                    if p > 2 then (Lp) is equal to:
                        (\tau(\chi))^(\sigma_n - n) is a generator of cyclic
                        group <\zeta_r>
                */
                if (p > 2 && state == 0 && unity_power > 0)
                {
                    ulong upow = unity_power;
                    /*
                        if gcd(r, unity_power) = 1 then
                        (\tau(\chi))^(\sigma_n - n) is a generator
                    */
                    if (n_gcd(r, upow) == 1)
                    {
                        state = 3;
                        lambdas[pind] = state;
                    }
                }

                /*
                    (Lp.b)
                    check 2) of (Lp) if p == 2 and nmod4 == 3
                */
                if (p == 2 && unity_power > 0
                    && (state == 0 || state == 1) && nmod4 == 3)
                {
                    ulong upow = unity_power;
                    if (n_gcd(r, upow) == 1)
                    {
                        if (state == 0)
                        {
                            state = 2;
                            lambdas[pind] = state;
                        }
                        if (state == 1)
                        {
                            state = 3;
                            lambdas[pind] = state;
                        }
                    }
                }
            }
        }
    }

    /*
        if for some p we have not proved (Lp)
        then n can be as prime or composite
    */
    if (result == PROBABPRIME)
        for (i = 0; i < config->rs.num; i++)
            if (lambdas[i] != 3)
                result = UNKNOWN;


    /* if n can be prime we do final division */
    if (result == UNKNOWN || result == PROBABPRIME)
    {
        int f_division;
        f_division = aprcl_is_prime_final_division(n, config->s, config->R);
        /* if (Lp) is true for all p | R and f_division == 1 then n - prime */
        if (result == PROBABPRIME && f_division == 1)
            result = PRIME;
        /*
            if we not prove (Lp) for some p and f_division == 1
            then n still can be prime
        */
        if (result == UNKNOWN && f_division == 1)
            result = PROBABPRIME;
        /* if f_division == 0 so we find the divisor of n, so n - composite */
        if (f_division == 0)
            result = COMPOSITE;
    }

    flint_free(lambdas);

    return result;
}

int
aprcl_is_prime_gauss_min_R(const fmpz_t n, ulong R)
{
    primality_test_status result;
    aprcl_config config;

    aprcl_config_gauss_init_min_R(config, n, R);
    result = _aprcl_is_prime_gauss(n, config);

    aprcl_config_gauss_clear(config);

    if (result == PRIME)
        return 1;
    return 0;
}

int
aprcl_is_prime_gauss(const fmpz_t n)
{
    ulong R;
    primality_test_status result;
    aprcl_config config;

    if (fmpz_cmp_ui(n, 2) < 0)
       return 0;

    aprcl_config_gauss_init_min_R(config, n, 180);
    result = _aprcl_is_prime_gauss(n, config);
    R = config->R;
    aprcl_config_gauss_clear(config);

    /*
        if result == PROBABPRIME it means that we have
        not proved (Lp) for some p (most likely we fail L.c step);
        we can try to use bigger R
    */
    if (result == PROBABPRIME)
    {
        R = R * 2;
        aprcl_config_gauss_init_min_R(config, n, R);
        result = _aprcl_is_prime_gauss(n, config);
        aprcl_config_gauss_clear(config);
    }

    if (result == PROBABPRIME)
    {
        R = R * 3;
        aprcl_config_gauss_init_min_R(config, n, R);
        result = _aprcl_is_prime_gauss(n, config);
        aprcl_config_gauss_clear(config);
    }

    if (result == PROBABPRIME)
    {
        R = R * 5;
        aprcl_config_gauss_init_min_R(config, n, R);
        result = _aprcl_is_prime_gauss(n, config);
        aprcl_config_gauss_clear(config);
    }

    if (result == PROBABPRIME || result == UNKNOWN)
    {
        flint_throw(FLINT_ERROR, "aprcl_is_prime_gauss: failed to prove n prime for n = %s\n", fmpz_get_str(NULL, 10, n));
    }

    if (result == PRIME)
        return 1;
    return 0;
}
