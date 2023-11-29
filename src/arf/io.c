/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2015 Arb authors
    Copyright (C) 2019 Julian Rüth

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string.h>
#include "arf.h"
#include "arb.h"

/* strings ********************************************************************/

static void
arf_set_fmpz_2exp_dump(arf_t x, const fmpz_t m, const fmpz_t e) {
    if (fmpz_is_zero(m)) {
        if (fmpz_get_si(e) == 0) arf_zero(x);
        else if (fmpz_get_si(e) == -1) arf_pos_inf(x);
        else if (fmpz_get_si(e) == -2) arf_neg_inf(x);
        else if (fmpz_get_si(e) == -3) arf_nan(x);
        else
        {
            /* Impossible to happen; all the special values have been treated above. */
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);
        }
        return;
    }

    arf_set_fmpz_2exp(x, m, e);
}

static void
arf_get_fmpz_2exp_dump(fmpz_t m, fmpz_t e, const arf_t x) {
    if (arf_is_special(x))
    {
        fmpz_zero(m);
        if (arf_is_zero(x)) fmpz_zero(e);
        else if (arf_is_pos_inf(x)) fmpz_set_si(e, -1);
        else if (arf_is_neg_inf(x)) fmpz_set_si(e, -2);
        else if (arf_is_nan(x)) fmpz_set_si(e, -3);
        else
        {
            /* Impossible to happen; all the special values have been treated above. */
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);
        }
        return;
    }

    arf_get_fmpz_2exp(m, e, x);
}

char * arf_get_str(const arf_t x, slong d)
{
    if (arf_is_special(x))
    {
        char * s = flint_malloc(5);

        if (arf_is_zero(x))
            strcpy(s, "0");
        else if (arf_is_pos_inf(x))
            strcpy(s, "+inf");
        else if (arf_is_neg_inf(x))
            strcpy(s, "-inf");
        else
            strcpy(s, "nan");

        return s;
    }
    else
    {
        arb_t t;
        *arb_midref(t) = *x;
        mag_init(arb_radref(t));  /* no need to free */
        return arb_get_str(t, FLINT_MAX(d, 1), ARB_STR_NO_RADIUS);
    }
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

char *
arf_dump_str(const arf_t x)
{
    size_t res_len;
    char * res;

    fmpz_t mantissa, exponent;

    fmpz_init(mantissa);
    fmpz_init(exponent);

    arf_get_fmpz_2exp_dump(mantissa, exponent, x);

    res_len = (fmpz_sgn(mantissa) < 0) + fmpz_sizeinbase(mantissa, 16) + 1
        + (fmpz_sgn(exponent) < 0) + fmpz_sizeinbase(exponent, 16);
    res = (char*)flint_malloc(res_len + 1);

    fmpz_get_str(res, 16, mantissa);
    strcat(res, " ");
    fmpz_get_str(res + strlen(res), 16, exponent);

    fmpz_clear(mantissa);
    fmpz_clear(exponent);

    if (strlen(res) != res_len)
        flint_throw(FLINT_ERROR, "(%s): strlen(res) != res_len\n", __func__);

    return res;
}

/* printing *******************************************************************/

void
arf_fprint(FILE * file, const arf_t x)
{
    if (arf_is_normal(x))
    {
        fmpz_t man, exp;

        fmpz_init(man);
        fmpz_init(exp);

        arf_get_fmpz_2exp(man, exp, x);

        flint_fprintf(file, "(");
        fmpz_fprint(file, man);
        flint_fprintf(file, " * 2^");
        fmpz_fprint(file, exp);
        flint_fprintf(file, ")");

        fmpz_clear(man);
        fmpz_clear(exp);
    }
    else
    {
        if (arf_is_zero(x)) flint_fprintf(file, "(0)");
        else if (arf_is_pos_inf(x)) flint_fprintf(file, "(+inf)");
        else if (arf_is_neg_inf(x)) flint_fprintf(file, "(-inf)");
        else flint_fprintf(file, "(nan)");
    }
}

void
arf_fprintd(FILE * file, const arf_t x, slong d)
{
    char * s = arf_get_str(x, d);

    fprintf(file, "%s", s);
    flint_free(s);
}

void arf_print(const arf_t x) { arf_fprint(stdout, x); }
void arf_printd(const arf_t y, slong d) { arf_fprintd(stdout, y, d); }

/* file I/O *******************************************************************/

int arf_load_file(arf_t x, FILE * stream)
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

int
arf_dump_file(FILE * stream, const arf_t x)
{
    int nwrite;
    char* data = arf_dump_str(x);
    nwrite = fputs(data, stream);
    if (nwrite == EOF)
        return nwrite;

    flint_free(data);
    return 0;
}
