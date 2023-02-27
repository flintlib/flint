/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


static void _fmpz_binpow_fmpz(fmpz_t a, const fmpz * p, const fmpz_t e)
{
    ulong i;
    flint_bitcnt_t bits = fmpz_bits(e);

    fmpz_one(a);
    for (i = 0; i < bits; i++)
    {
        if (fmpz_tstbit(e, i))
        {
            fmpz_mul(a, a, p + 2*i);
        }
    }
}

/*

             scales+0          scales+n-1
--------               e1                  en
 \            (a1     )          (an      )
  \       c * (-- * x1)  * ... * (-- * x_n)
  /           (b1     )          (bn      )
 /
--------
  c, e

let fi and di be the min and max exponents of xi

rewrite as integer polynomial

--------
 \              e1-f1     d1-e1          en-f1     dn-en     e1          en
  \       c * a1      * b1      * ... * an      * bn      * x1 * ... * xn
  /
 /
--------
  c, e

times

  f1           fn
a1           an
---- * ... * ----
  d1           dn
b1           bn

*/

int _fmpq_mpoly_rescale(fmpq_t Acontent, fmpz * Acoeff, const fmpq_mpoly_t B,
                               const fmpq * scales, const fmpq_mpoly_ctx_t ctx)
{
    int success = 1;
    slong i, j, v;
    slong nvars, N;
    slong Blen;
    fmpz * Bcoeff;
    ulong * Bexp;
    fmpz * Buexp;
    flint_bitcnt_t Bbits, tbits;
    fmpz * emin;
    fmpz * emax;
    fmpz ** powertable;
    fmpz * lastexp;
    fmpz * lastpow;
    fmpz_t t, s, u;
    fmpq_t q;
    TMP_INIT;

    TMP_START;

    fmpz_init(t);
    fmpz_init(s);
    fmpz_init(u);
    fmpq_init(q);

    nvars = ctx->zctx->minfo->nvars;
    Blen = B->zpoly->length;
    Bcoeff = B->zpoly->coeffs;
    Bexp = B->zpoly->exps;
    Bbits = B->zpoly->bits;

    FLINT_ASSERT(Blen > 0);

    N = mpoly_words_per_exp(Bbits, ctx->zctx->minfo);

    Buexp = _fmpz_vec_init(nvars*Blen);

    i = 0;
    mpoly_get_monomial_ffmpz(Buexp + nvars*i, Bexp + N*i, Bbits, ctx->zctx->minfo);
    emin = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    emax = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (v = 0; v < nvars; v++)
    {
        fmpz_init_set(emin + v, Buexp + nvars*i + v);
        fmpz_init_set(emax + v, Buexp + nvars*i + v);
    }
    for (i = 1; i < Blen; i++)
    {
        mpoly_get_monomial_ffmpz(Buexp + nvars*i, Bexp + N*i, Bbits, ctx->zctx->minfo);
        for (v = 0; v < nvars; v++)
        {
            if (fmpz_cmp(emin + v, Buexp + nvars*i + v) > 0)
                fmpz_set(emin + v, Buexp + nvars*i + v);
            if (fmpz_cmp(emax + v, Buexp + nvars*i + v) < 0)
                fmpz_set(emax + v, Buexp + nvars*i + v);
        }
    }

    /* check that the powers in the powertable loop do not get too big */
    for (v = 0; v < nvars; v++)
    {
        if (_fmpz_pow_fmpz_is_not_feasible(fmpq_height_bits(scales + v), emax + v))
        {
            success = 0;
            goto cleanup_exp_bounds;
        }
    }

    powertable = (fmpz **) TMP_ALLOC(nvars*sizeof(fmpz *));
    for (v = 0; v < nvars; v++)
    {
        fmpz * row;

        tbits = fmpz_bits(emax + v);
        tbits = FLINT_MAX(UWORD(1), tbits);
        row = (fmpz *) flint_malloc(2*tbits*sizeof(fmpz));
        powertable[v] = row;

        j = 0;
        fmpz_init_set(row + 2*j + 0, fmpq_numref(scales + v));
        fmpz_init_set(row + 2*j + 1, fmpq_denref(scales + v));
        for (j = 1; j < tbits; j++)
        {
            fmpz_init(row + 2*j + 0);
            fmpz_init(row + 2*j + 1);
            fmpz_mul(row + 2*j + 0, row + 2*j - 2, row + 2*j - 2);
            fmpz_mul(row + 2*j + 1, row + 2*j - 1, row + 2*j - 1);
        }
    }

    lastexp = (fmpz *) TMP_ALLOC(2*nvars*sizeof(fmpz));
    lastpow = (fmpz *) TMP_ALLOC(2*nvars*sizeof(fmpz));
    for (v = 0; v < nvars; v++)
    {
        fmpz_init_set_ui(lastexp + 2*v + 0, 0);
        fmpz_init_set_ui(lastpow + 2*v + 0, 1);
        fmpz_init_set_ui(lastexp + 2*v + 1, 0);
        fmpz_init_set_ui(lastpow + 2*v + 1, 1);
    }

    fmpq_set(Acontent, B->content);
    for (v = 0; v < nvars; v++)
    {
        _fmpz_binpow_fmpz(fmpq_numref(q), powertable[v] + 0, emin + v);
        _fmpz_binpow_fmpz(fmpq_denref(q), powertable[v] + 1, emax + v);
        fmpq_mul(Acontent, Acontent, q);
    }
    for (i = 0; i < Blen; i++)
    {
        fmpz_set(Acoeff + i, Bcoeff + i);
        for (v = 0; v < nvars; v++)
        {
            /* power numerator of vth scale */
            fmpz_sub(t, Buexp + nvars*i + v, emin + v);
            FLINT_ASSERT(fmpz_sgn(t) >= 0);
            fmpz_sub(s, t, lastexp + 2*v + 0);
            if (fmpz_sgn(s) >= 0)
            {
                _fmpz_binpow_fmpz(u, powertable[v] + 0, s);
                fmpz_mul(lastpow + 2*v + 0, lastpow + 2*v + 0, u);
            }
            else
            {
                _fmpz_binpow_fmpz(lastpow + 2*v + 0, powertable[v] + 0, t);
            }
            fmpz_swap(lastexp + 2*v + 0, t);
            fmpz_mul(Acoeff + i, Acoeff + i, lastpow + 2*v + 0);

            /* power denominator of vth scale */
            fmpz_sub(t, emax + v, Buexp + nvars*i + v);
            FLINT_ASSERT(fmpz_sgn(t) >= 0);
            fmpz_sub(s, t, lastexp + 2*v + 1);
            if (fmpz_sgn(s) >= 0)
            {
                _fmpz_binpow_fmpz(u, powertable[v] + 1, s);
                fmpz_mul(lastpow + 2*v + 1, lastpow + 2*v + 1, u);
            }
            else
            {
                _fmpz_binpow_fmpz(lastpow + 2*v + 1, powertable[v] + 1, t);
            }
            fmpz_swap(lastexp + 2*v + 1, t);
            fmpz_mul(Acoeff + i, Acoeff + i, lastpow + 2*v + 1);
        }
    }

    /* clean up */
    for (v = 0; v < nvars; v++)
    {
        fmpz_clear(lastexp + 2*v + 0);
        fmpz_clear(lastpow + 2*v + 0);
        fmpz_clear(lastexp + 2*v + 1);
        fmpz_clear(lastpow + 2*v + 1);
    }
    for (v = 0; v < nvars; v++)
    {
        tbits = fmpz_bits(emax + v);
        tbits = FLINT_MAX(UWORD(1), tbits);
        for (j = 0; j < tbits; j++)
        {
            fmpz_clear(powertable[v] + 2*j + 0);
            fmpz_clear(powertable[v] + 2*j + 1);
        }
        flint_free(powertable[v]);
    }

cleanup_exp_bounds:

    for (v = 0; v < nvars; v++)
    {
        fmpz_clear(emin + v);
        fmpz_clear(emax + v);
    }

    _fmpz_vec_clear(Buexp, nvars*Blen);

    fmpz_clear(t);
    fmpz_clear(s);
    fmpz_clear(u);
    fmpq_clear(q);

    TMP_END;

    return success;
}

