/*
    Copyright (C) 2019 Julian Rüth

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>

#include "arf.h"

static void
arf_set_fmpz_2exp_dump(arf_t x, const fmpz_t m, const fmpz_t e) {
    if (fmpz_is_zero(m)) {
        if (fmpz_get_si(e) == 0) arf_zero(x);
        else if (fmpz_get_si(e) == -1) arf_pos_inf(x);
        else if (fmpz_get_si(e) == -2) arf_neg_inf(x);
        else if (fmpz_get_si(e) == -3) arf_nan(x);
        else {
            /* Impossible to happen; all the special values have been treated above. */
            flint_abort();
        }
        return;
    }

    arf_set_fmpz_2exp(x, m, e);
}

int
arf_load_str(arf_t x, const char* data)
{
    fmpz_t mantissa, exponent;
    char * e_str;
    char * m_str;
    int err = 0;

    fmpz_init(mantissa);
    fmpz_init(exponent);

    e_str = strchr(data, ' ');
    if (e_str == NULL) return 1;

    m_str = (char*)flint_malloc(e_str - data + 1);
    strncpy(m_str, data, e_str - data);
    m_str[e_str - data] = '\0';
    e_str++;

    err = fmpz_set_str(mantissa, m_str, 16);

    flint_free(m_str);

    if (err)
    {
        fmpz_clear(exponent);
        fmpz_clear(mantissa);
        return err;
    }

    err = fmpz_set_str(exponent, e_str, 16);

    if (err)
    {
        fmpz_clear(exponent);
        fmpz_clear(mantissa);
        return err;
    }

    arf_set_fmpz_2exp_dump(x, mantissa, exponent);

    fmpz_clear(exponent);
    fmpz_clear(mantissa);

    return err;
}

int arf_load_file(arf_t x, FILE* stream)
{
    fmpz_t mantissa, exponent;
    __mpz_struct *mpz_mantissa, *mpz_exponent;
    int err;

    fmpz_init(mantissa);
    fmpz_init(exponent);

    mpz_mantissa = _fmpz_promote(mantissa);
    mpz_exponent = _fmpz_promote(exponent);

    err = 0;

    if (mpz_inp_str(mpz_mantissa, stream, 16) == 0)
        err = 1;

    if (!err && mpz_inp_str(mpz_exponent, stream, 16) == 0)
        err = 1;

    _fmpz_demote_val(mantissa);
    _fmpz_demote_val(exponent);

    if (!err)
        arf_set_fmpz_2exp_dump(x, mantissa, exponent);

    fmpz_clear(mantissa);
    fmpz_clear(exponent);

    return err;
}

