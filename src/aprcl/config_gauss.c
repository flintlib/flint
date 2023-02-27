/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

static void
_aprcl_config_gauss_update(aprcl_config conf)
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
        prime++;
        while (n_is_prime(prime) == 0)
            prime++;
    }
}

void
aprcl_config_gauss_init(aprcl_config conf, const fmpz_t n)
{
    fmpz_t s2;

    fmpz_init_set_ui(s2, 0);
    fmpz_init(conf->s);
    fmpz_factor_init(conf->qs);
    conf->R = 1;

    while (fmpz_cmp(s2, n) <= 0)
    {
        conf->R += 1;
        _aprcl_config_gauss_update(conf);
        fmpz_mul(s2, conf->s, conf->s);
    }
    n_factor_init(&conf->rs);
    n_factor(&conf->rs, conf->R, 1);

    conf->qs_used = NULL;  /* not used */

    fmpz_clear(s2);
}

void
aprcl_config_gauss_init_min_R(aprcl_config conf, const fmpz_t n, ulong R)
{
    fmpz_t s2;

    fmpz_init_set_ui(s2, 0);
    fmpz_init(conf->s);
    fmpz_factor_init(conf->qs);
    conf->R = R - 1;

    while (fmpz_cmp(s2, n) <= 0)
    {
        conf->R += 1;
        _aprcl_config_gauss_update(conf);
        fmpz_mul(s2, conf->s, conf->s);
    }

    n_factor_init(&conf->rs);
    n_factor(&conf->rs, conf->R, 1);

    conf->qs_used = NULL;  /* not used */

    fmpz_clear(s2);
}

void
aprcl_config_gauss_clear(aprcl_config conf)
{
    fmpz_clear(conf->s);
    fmpz_factor_clear(conf->qs);
}
