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

void _aprcl_config_update(aprcl_config conf)
{
    ulong prime = 2;

    fmpz_set_ui(conf->s, 1);
    fmpz_factor_clear(conf->qs);
    fmpz_factor_init(conf->qs);
    conf->qs->sign = 1;

    while (2 * (prime - 1) <= conf->R)
    {
        if ((conf->R % (prime - 1)) == 0)
        {
            _fmpz_factor_append_ui(conf->qs, prime, 1);
            fmpz_mul_ui(conf->s, conf->s, prime);
        }
        prime = n_nextprime(prime, 1);
    }
}

void aprcl_config_init(aprcl_config conf, const fmpz_t n)
{
    fmpz_t s2;

    fmpz_init_set_ui(s2, 0);
    fmpz_init(conf->s);
    fmpz_factor_init(conf->qs);
    conf->R = 1;

    while (fmpz_cmp(s2, n) <= 0)
    {
        conf->R += 1;
        _aprcl_config_update(conf);
        fmpz_mul(s2, conf->s, conf->s);
    }
}

void aprcl_config_clear(aprcl_config conf)
{
    fmpz_clear(conf->s);
    fmpz_factor_clear(conf->qs);
}


