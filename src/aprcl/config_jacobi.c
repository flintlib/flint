/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef __GNUC__
# define log __builtin_log
#else
# include <math.h>
#endif

#include "fmpz.h"
#include "fmpz_factor.h"
#include "aprcl.h"

ulong
aprcl_R_value(const fmpz_t n)
{
    ulong bits = fmpz_bits(n);

    if (bits <= 17) return 6;           /*  2 * 3                       */
    if (bits <= 31) return 12;          /*  2^2 * 3                     */
    if (bits <= 54) return 36;          /*  2^2 * 3^2                   */
    if (bits <= 68) return 72;          /*  2^3 * 3^2                   */
    if (bits <= 101) return 180;        /*  2^2 * 3^2 * 5               */
    if (bits <= 127) return 360;        /*  2^3 * 3^2 * 5               */
    if (bits <= 152) return 720;        /*  2^4 * 3^2 * 5               */
    if (bits <= 204) return 1260;       /*  2^2 * 3^2 * 5 * 7           */
    if (bits <= 268) return 2520;       /*  2^3 * 3^2 * 5 * 7           */
    if (bits <= 344) return 5040;       /*  2^4 * 3^2 * 5 * 7           */
    if (bits <= 525) return 27720;      /*  2^3 * 3^2 * 5 * 7 * 11      */
    if (bits <= 650) return 55440;      /*  2^4 * 3^2 * 5 * 7 * 11      */
    if (bits <= 774) return 110880;     /*  2^5 * 3^2 * 5 * 7 * 11      */
    if (bits <= 1566) return 720720;    /*  2^4 * 3^2 * 5 * 7 * 11 * 13 */
    if (bits <= 1999) return 1441440;   /*  2^5 * 3^2 * 5 * 7 * 11 * 13 */
    if (bits <= 2096) return 1663200;
    if (bits <= 2165) return 1965600;
    if (bits <= 2321) return 2162160;
    if (bits <= 2377) return 2827440;
    if (bits <= 2514) return 3326400;
    if (bits <= 2588) return 3341520;
    if (bits <= 2636) return 3603600;
    if (bits <= 3028) return 4324320;
    if (bits <= 3045) return 5654880;
    if (bits <= 3080) return 6652800;
    if (bits <= 3121) return 6683040;
    if (bits <= 3283) return 7207200;
    if (bits <= 3491) return 8648640;   /*  2^6 * 3^3 * 5 * 7 * 11 * 13 */
    if (bits <= 3726) return 10810800;
    if (bits <= 3818) return 12972960;
    if (bits <= 3977) return 14414400;
    if (bits <= 4762) return 21621600;
    if (bits <= 5068) return 36756720;
    if (bits <= 5658) return 43243200;
    if (bits <= 5960) return 64864800;
    if (bits <= 6423) return 73513440;
    if (bits <= 6900) return 122522400;
    if (bits <= 9977) return 367567200;
    if (bits <= 12713) return 1396755360;

#if FLINT64
    /* 2^5 * 3^3 * 5^2 * 7 * 11 * 13 * 17 * 19 */
    return UWORD(6983776800);
#else
    flint_throw(FLINT_ERROR, "APRCL not supported for huge numbers on 32 bits\n");
    return 0;
#endif
}

static void
_aprcl_config_jacobi_reduce_s2(aprcl_config conf, const fmpz_t n)
{
    ulong i, j, q;
    double * w;
    n_factor_t q_factors;
    fmpz_t new_s, p;

    fmpz_init(new_s);
    fmpz_init(p);
    w = (double *) flint_malloc(sizeof(double) * conf->qs->num);

    for (i = 0; i < conf->qs->num; i++)
    {
        conf->qs_used[i] = 1;

        q = fmpz_get_ui(conf->qs->p + i);
        n_factor_init(&q_factors);
        n_factor(&q_factors, q - 1, 1);

        w[i] = 0;

        for (j = 0; j < q_factors.num; j++)
        {
            ulong p, euler_phi;

            p = q_factors.p[j];
            euler_phi = n_pow(p, q_factors.exp[j] - 1) * (p - 1);
            euler_phi = euler_phi * euler_phi;

            w[i] += euler_phi;
        }

        w[i] /= log((double) n_pow(q, conf->qs->exp[i]));
    }

    while (1)
    {
        double w_max;
        slong ind;

        w_max = -1;
        ind = -1;

        for (i = 0; i < conf->qs->num; i++)
        {
            if (conf->qs_used[i] == 0)
                continue;

            fmpz_pow_ui(p, conf->qs->p + i, conf->qs->exp[i]);
            fmpz_fdiv_q(new_s, conf->s, p);
            fmpz_mul(new_s, new_s, new_s);

            if (fmpz_cmp(new_s, n) > 0)
                if (w_max <= w[i])
                {
                    w_max = w[i];
                    ind = i;
                }
        }

        if (ind == -1)
            break;

        fmpz_pow_ui(p, conf->qs->p + ind, conf->qs->exp[ind]);
        fmpz_fdiv_q(new_s, conf->s, p);
        fmpz_set(conf->s, new_s);
        conf->qs_used[ind] = 0;
    }

    fmpz_clear(new_s);
    fmpz_clear(p);
    flint_free(w);
}

static void
_aprcl_config_jacobi_update(aprcl_config conf)
{
    ulong prime = 2;

    fmpz_set_ui(conf->s, 1);
    fmpz_factor_clear(conf->qs);
    fmpz_factor_init(conf->qs);
    conf->qs->sign = 1;

    _fmpz_factor_append_ui(conf->qs, prime, aprcl_p_power_in_q(conf->R, prime) + 2);
    fmpz_mul_ui(conf->s, conf->s, n_pow(prime, aprcl_p_power_in_q(conf->R, prime) + 2));

    prime = 3;
    while (2 * (prime - 1) <= conf->R)
    {
        if ((conf->R % (prime - 1)) == 0)
        {
            _fmpz_factor_append_ui(conf->qs, prime, aprcl_p_power_in_q(conf->R, prime) + 1);
            fmpz_mul_ui(conf->s, conf->s, n_pow(prime, aprcl_p_power_in_q(conf->R, prime) + 1));
        }
        prime++;
        while (n_is_prime(prime) == 0)
            prime++;
    }

    if (n_is_prime(conf->R + 1))
    {
        _fmpz_factor_append_ui(conf->qs, conf->R + 1, 1);
        fmpz_mul_ui(conf->s, conf->s, conf->R + 1);
    }
}

/* Computes s = \prod q^(k + 1) ; q - prime, q - 1 | R; q^k | R and q^(k + 1) not | R */
void
aprcl_config_jacobi_init(aprcl_config conf, const fmpz_t n)
{
    fmpz_init(conf->s);
    fmpz_factor_init(conf->qs);
    conf->R = aprcl_R_value(n);
    _aprcl_config_jacobi_update(conf);

    n_factor_init(&conf->rs);
    n_factor(&conf->rs, conf->R, 1);

    conf->qs_used = (int *) flint_malloc(sizeof(int) * conf->qs->num);
    _aprcl_config_jacobi_reduce_s2(conf, n);
}

void
aprcl_config_jacobi_clear(aprcl_config conf)
{
    fmpz_clear(conf->s);
    fmpz_factor_clear(conf->qs);
    flint_free(conf->qs_used);
}
