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

ulong
_R_value(const fmpz_t n)
{
    ulong bits = fmpz_bits(n);

    if (bits <= 101) return 180;
    if (bits <= 152) return 720;
    if (bits <= 204) return 1260;
    if (bits <= 268) return 2520;
    if (bits <= 344) return 5040;
    if (bits <= 525) return 27720;
    if (bits <= 774) return 98280;
    if (bits <= 1035) return 166320;
    if (bits <= 1566) return 720720;
    if (bits <= 2082) return 1663200;
    if (bits <= 3491) return 8648640;
    return 6983776800;
}

void _jacobi_config_update(aprcl_config conf)
{
    ulong prime = 2;

    fmpz_set_ui(conf->s, 1);
    fmpz_factor_clear(conf->qs);
    fmpz_factor_init(conf->qs);
    conf->qs->sign = 1;

    _fmpz_factor_append_ui(conf->qs, prime, p_power_in_q(conf->R, prime) + 2);
    fmpz_mul_ui(conf->s, conf->s, n_pow(prime, p_power_in_q(conf->R, prime) + 2));

    prime = 3;
    while (2 * (prime - 1) <= conf->R)
    {
        if ((conf->R % (prime - 1)) == 0)
        {
            _fmpz_factor_append_ui(conf->qs, prime, p_power_in_q(conf->R, prime) + 1);
            fmpz_mul_ui(conf->s, conf->s, n_pow(prime, p_power_in_q(conf->R, prime) + 1));
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
void jacobi_config_init(aprcl_config conf, const fmpz_t n)
{
    fmpz_init(conf->s);
    fmpz_factor_init(conf->qs);
    conf->R = _R_value(n);
    _jacobi_config_update(conf);

    n_factor_init(&conf->rs);
    n_factor(&conf->rs, conf->R, 1);
}

void jacobi_config_clear(aprcl_config conf)
{
    fmpz_clear(conf->s);
    fmpz_factor_clear(conf->qs);
}

