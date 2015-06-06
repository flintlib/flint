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

int is_condition_one()
{
    return 0;
}

/* if for some i from 0 to p-1 \tau(\chi^n) == \tau^n(\chi)*\zeta_p^i return 1; otherwise return 0 */
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

int is_prime_gauss(const fmpz_t n)
{
    ulong i, j, k;
    aprcl_config conf;
    aprcl_config_init(conf, n);

    for (i = 0; i < conf->qs->num; i++)
    {
        ulong q, r;
        n_factor_t q_factors;

        q = conf->qs->p[i];
        n_factor_init(&q_factors);
        n_factor(&q_factors, q, 1);

        for (j = 0; j <= q_factors.num; j++)
        {
            for (k = 1; k <= q_factors.exp[j]; k++)
            {
                unity_zpq gauss;
                fmpz_t qr, gcd;

                fmpz_init(qr);
                fmpz_init(gcd);
                r = n_pow(q_factors.p[j], k);
                fmpz_set_ui(qr, q * r);


                fmpz_gcd(gcd, n, qr);
                if (fmpz_cmp_ui(gcd, 1) == 1)
                {
                    return 0;
                }
/*
                if (is_condition_one() == 0)
                {
                    return 0;
                }
*/                
/*                unity_zpq_init(gauss, q, r, n);
                unity_zpq_gauss_sum(gauss, q, r);
                if (is_condition_two() == 0)
                {
                    return 0;
                }
*/
                unity_zpq_clear(gauss);
                fmpz_clear(qr);
                fmpz_clear(gcd);
            }
        }
    }


    aprcl_config_clear(conf);
    return 1;
}

