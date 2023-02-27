/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

void
unity_zpq_pow(unity_zpq f, const unity_zpq g, const fmpz_t pow)
{
    unity_zpq value;
    fmpz_t power, rem;

    unity_zpq_init(value, f->q, f->p, fmpz_mod_ctx_modulus(f->ctx));
    fmpz_init_set(power, pow);
    fmpz_init(rem);

    unity_zpq_coeff_set_ui(f, 0, 0, 1);

    unity_zpq_copy(value, g);

    while (fmpz_is_zero(power) == 0)
    {
        unity_zpq temp_pow;

        fmpz_fdiv_r_2exp(rem, power, 1);
        if (fmpz_is_zero(rem) == 0)
        {
            unity_zpq temp;
            unity_zpq_init(temp, f->q, f->p, fmpz_mod_ctx_modulus(f->ctx));

            unity_zpq_mul(temp, f, value);
            unity_zpq_swap(f, temp);

            unity_zpq_clear(temp);
        }

        unity_zpq_init(temp_pow, f->q, f->p, fmpz_mod_ctx_modulus(f->ctx));
        unity_zpq_mul(temp_pow, value, value);
        unity_zpq_swap(value, temp_pow);
        fmpz_fdiv_q_2exp(power, power, 1);

        unity_zpq_clear(temp_pow);
    }


    fmpz_clear(power);
    fmpz_clear(rem);
    unity_zpq_clear(value);
}

void
unity_zpq_pow_ui(unity_zpq f, const unity_zpq g, ulong pow)
{
    fmpz_t p;

    fmpz_init_set_ui(p, pow);

    unity_zpq_pow(f, g, p);

    fmpz_clear(p);
}

