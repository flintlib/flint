/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2010, 2012 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson
    Authored 2015 by Daniel S. Roche; US Government work in the public domain.

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_randbits(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
{
    if (bits <= SMALL_FMPZ_BITCOUNT_MAX)
    {
        _fmpz_demote(f);
        *f = n_randbits(state, bits);
        if (n_randint(state, 2))
            *f = -*f;
    }
    else
    {
        __mpz_struct *mf = _fmpz_promote(f);
        _flint_rand_init_gmp(state);
        mpz_urandomb(mf, state->gmp_state, bits);
        mpz_setbit(mf, bits - 1);

        if (n_randint(state, 2))
            mpz_neg(mf, mf);

        _fmpz_demote_val(f);
    }
}

void
fmpz_randm(fmpz_t f, flint_rand_t state, const fmpz_t m)
{
    flint_bitcnt_t bits = fmpz_bits(m);
    int sgn = fmpz_sgn(m);

    if (bits <= SMALL_FMPZ_BITCOUNT_MAX)
    {
        _fmpz_demote(f);
        *f =  (sgn >= 0) ? n_randint(state, *m) : - n_randint(state, -(*m));
    }
    else
    {
        __mpz_struct * mf = _fmpz_promote(f);

        _flint_rand_init_gmp(state);
        mpz_urandomm(mf, state->gmp_state, COEFF_TO_PTR(*m));
        if (sgn < 0)
            mpz_neg(mf, mf);
        _fmpz_demote_val(f);
    }
}

void fmpz_randprime(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits, int proved)
{
    if (bits <= SMALL_FMPZ_BITCOUNT_MAX)
    {
        _fmpz_demote(f);
        *f = n_randprime(state, bits, proved);
    }
    else
    {
        /* Here I would like to just call
         * fmpz_randbits(f, state, bits);
         * but it has different semantics from n_randbits,
         * and in particular may return integers with fewer bits.
         */
        __mpz_struct * mf = _fmpz_promote(f);
        _flint_rand_init_gmp(state);

        do
        {
            mpz_urandomb(mf, state->gmp_state, bits - 1);
            mpz_setbit(mf, bits - 1);

            fmpz_nextprime(f, f, proved);
        } while (fmpz_bits(f) != bits);
    }
}

void
fmpz_randtest(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
{
    ulong m;

    fmpz_randtest_unsigned(f, state, bits);

    m = n_randlimb(state);
    if (m & UWORD(1))
        fmpz_inplace_neg(f);
}

void
fmpz_randtest_unsigned(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
{
    ulong m;

    m    = n_randlimb(state);
    bits = n_randint(state, bits + 1);

    if (bits <= SMALL_FMPZ_BITCOUNT_MAX)
    {
        _fmpz_demote(f);
        if (m & UWORD(3))
            *f = n_randtest_bits(state, bits);
        else
        {
            m >>= 2;
            if (bits == 0)
                *f = 0;
            else if (bits < SMALL_FMPZ_BITCOUNT_MAX)
                *f = m & UWORD(1);
            else
                *f = COEFF_MAX;
        }
    }
    else
    {
        __mpz_struct * mf = _fmpz_promote(f);

        _flint_rand_init_gmp(state);
        mpz_rrandomb(mf, state->gmp_state, bits);
        _fmpz_demote_val(f);
    }
}

void
fmpz_randtest_not_zero(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
{
    if (bits == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_randtest_not_zero). bits == 0.\n");
    }

    fmpz_randtest(f, state, bits);
    if (fmpz_is_zero(f))
        fmpz_one(f);
}

void
fmpz_randtest_mod(fmpz_t f, flint_rand_t state, const fmpz_t m)
{
    fmpz_t t;

    fmpz_init(t);
    fmpz_randtest_unsigned(t, state, fmpz_bits(m) + 2);
    fmpz_mod(t, t, m);

    if (n_randlimb(state) & UWORD(1))
    {
        fmpz_sub(t, m, t);
        fmpz_sub_ui(t, t, UWORD(1));
    }

    fmpz_set(f, t);
    fmpz_clear(t);
}

void
fmpz_randtest_mod_signed(fmpz_t f, flint_rand_t state, const fmpz_t m)
{
    /* Randomly generate m/2 when included in the range */
    if ((n_randlimb(state) % 32 == 1) && (fmpz_fdiv_ui(m, 2) == 0))
    {
        fmpz_fdiv_q_ui(f, m, UWORD(2));
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_tdiv_q_ui(t, m, UWORD(2));
        fmpz_randtest_mod(t, state, t);
        if (n_randlimb(state) & UWORD(1))
        {
            fmpz_inplace_neg(t);
        }
        fmpz_set(f, t);
        fmpz_clear(t);
    }
}
