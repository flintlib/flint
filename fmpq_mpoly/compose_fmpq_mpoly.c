/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"



void _fmpz_binpow_fmpz(fmpz_t a, const fmpz * p, const fmpz_t e)
{
    ulong i;
    mp_bitcnt_t bits = fmpz_bits(e);
/*
printf("computing "); fmpz_print(p); printf("^"); fmpz_print(e); printf("\n");
*/
    fmpz_one(a);
    for (i = 0; i < bits; i++)
    {
        if (fmpz_tstbit(e, i))
        {
            fmpz_mul(a, a, p + 2*i);
        }
    }
/*
printf("result "); fmpz_print(a); printf("\n");
*/
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

void _fmpq_mpoly_rescale(
    fmpq_t Acontent,
    fmpz * Acoeff,
    const fmpq_mpoly_t B,
    const fmpq * scales,
    const fmpq_mpoly_ctx_t ctx)
{
    slong i, j, v;
    slong nvars, N;
    slong Blen;
    fmpz * Bcoeff;
    ulong * Bexp;
    fmpz * Buexp;
    mp_bitcnt_t bits;
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
    bits = B->zpoly->bits;

    FLINT_ASSERT(Blen > 0);

    N = mpoly_words_per_exp(bits, ctx->zctx->minfo);

    Buexp = _fmpz_vec_init(nvars*Blen);

    i = 0;
    mpoly_get_monomial_ffmpz(Buexp + nvars*i, Bexp + N*i, bits, ctx->zctx->minfo);
    emin = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    emax = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (v = 0; v < nvars; v++)
    {
        fmpz_init_set(emin + v, Buexp + nvars*i + v);
        fmpz_init_set(emax + v, Buexp + nvars*i + v);
    }
    for (i = 1; i < Blen; i++)
    {
        mpoly_get_monomial_ffmpz(Buexp + nvars*i, Bexp + N*i, bits, ctx->zctx->minfo);
        for (v = 0; v < nvars; v++)
        {
            if (fmpz_cmp(emin + v, Buexp + nvars*i + v) > 0)
                fmpz_set(emin + v, Buexp + nvars*i + v);
            if (fmpz_cmp(emax + v, Buexp + nvars*i + v) < 0)
                fmpz_set(emax + v, Buexp + nvars*i + v);
        }
    }

/*
for (v = 0; v < nvars; v++)
{
flint_printf("emin[%wd]: ",v); fmpz_print(emin + v); printf("\n");
flint_printf("emax[%wd]: ",v); fmpz_print(emax + v); printf("\n");
}
*/

    powertable = (fmpz **) TMP_ALLOC(nvars*sizeof(fmpz *));
    for (v = 0; v < nvars; v++)
    {
        mp_bitcnt_t bits;
        fmpz * row;

        bits = FLINT_MAX(UWORD(1), fmpz_bits(emax + v));
        row = (fmpz *) TMP_ALLOC(2*bits*sizeof(fmpz));
        powertable[v] = row;

        j = 0;
        fmpz_init_set(row + 2*j + 0, fmpq_numref(scales + v));
        fmpz_init_set(row + 2*j + 1, fmpq_denref(scales + v));
/*
flint_printf("powtable[%wd] + %wd: ",v,2*j + 0); fmpz_print(row + 2*j + 0); printf("\n");
flint_printf("powtable[%wd] + %wd: ",v,2*j + 1); fmpz_print(row + 2*j + 1); printf("\n");
*/
        for (j = 1; j < bits; j++)
        {
            fmpz_init(row + 2*j + 0);
            fmpz_init(row + 2*j + 1);
            fmpz_mul(row + 2*j + 0, row + 2*j - 2, row + 2*j - 2);
            fmpz_mul(row + 2*j + 1, row + 2*j - 1, row + 2*j - 1);
/*
flint_printf("powtable[%wd] + %wd: ",v,2*j + 0); fmpz_print(row + 2*j + 0); printf("\n");
flint_printf("powtable[%wd] + %wd: ",v,2*j + 1); fmpz_print(row + 2*j + 1); printf("\n");
*/
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
/*
flint_printf("q[%wd]: ",v); fmpq_print(q); printf("\n");
*/
        fmpq_mul(Acontent, Acontent, q);
    }
/*
printf("Acontent: "); fmpq_print(Acontent); printf("\n");
*/
    for (i = 0; i < Blen; i++)
    {
/*
flint_printf("\ni: %wd\n",i);
*/
        fmpz_set(Acoeff + i, Bcoeff + i);
/*
flint_printf("1 Acoeff[%wd]: ",i); fmpz_print(Acoeff + i); printf("\n");
*/
        for (v = 0; v < nvars; v++)
        {
/*
flint_printf("v: %wd\n",v);
*/
            /* power numerator of vth scale */
            fmpz_sub(t, Buexp + nvars*i + v, emin + v);
            FLINT_ASSERT(fmpz_sgn(t) >= 0);
            fmpz_sub(s, t, lastexp + 2*v + 0);
/*
flint_printf("lastpow + 2*v + 0: "); fmpz_print(lastpow + 2*v + 0); printf("\n");
*/
            if (fmpz_sgn(s) >= 0)
            {
                _fmpz_binpow_fmpz(u, powertable[v] + 0, s);
                fmpz_mul(lastpow + 2*v + 0, lastpow + 2*v + 0, u);
            }
            else
            {
                _fmpz_binpow_fmpz(lastpow + 2*v + 0, powertable[v] + 0, t);
            }
/*
flint_printf("lastpow + 2*v + 0: "); fmpz_print(lastpow + 2*v + 0); printf("\n");
*/
            fmpz_swap(lastexp + 2*v + 0, t);
/*
flint_printf("2 Acoeff[%wd]: ",i); fmpz_print(Acoeff + i); printf("\n");
flint_printf("lastpow + 2*v + 0: "); fmpz_print(lastpow + 2*v + 0); printf("\n");
*/
            fmpz_mul(Acoeff + i, Acoeff + i, lastpow + 2*v + 0);
/*
flint_printf("3 Acoeff[%wd]: ",i); fmpz_print(Acoeff + i); printf("\n");
*/
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
/*
flint_printf("4 Acoeff[%wd]: ",i); fmpz_print(Acoeff + i); printf("\n");
*/
            fmpz_mul(Acoeff + i, Acoeff + i, lastpow + 2*v + 1);
/*
flint_printf("5 Acoeff[%wd]: ",i); fmpz_print(Acoeff + i); printf("\n");
*/
        }
/*
flint_printf("6 Acoeff[%wd]: ",i); fmpz_print(Acoeff + i); printf("\n");
*/
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
        mp_bitcnt_t bits;
        bits = FLINT_MAX(UWORD(1), fmpz_bits(emax + v));
        for (j = 0; j < bits; j++)
        {
            fmpz_clear(powertable[v] + 2*j + 0);
            fmpz_clear(powertable[v]+ 2*j + 1);
        }
    }
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
}


void fmpq_mpoly_compose_fmpq_mpoly(fmpq_mpoly_t A,
                   const fmpq_mpoly_t B, fmpq_mpoly_struct * const * C,
                     const fmpq_mpoly_ctx_t ctxB, const fmpq_mpoly_ctx_t ctxAC)
{
    slong i;
    fmpq * scales;
    fmpz_mpoly_struct ** Czpoly;
    fmpz_mpoly_t newB;
    slong nvarsB = ctxB->zctx->minfo->nvars;
    TMP_INIT;

    if (fmpq_mpoly_is_zero(B, ctxB))
    {
        fmpq_poly_zero(A);
        return;
    }

    TMP_START;

    Czpoly = (fmpz_mpoly_struct **) TMP_ALLOC(nvarsB*sizeof(fmpz_mpoly_struct *));

    /*
        scale B by the contents of the polynomials in C
        We are only borrowing B to feed it to fmpz_mpoly_compose.
        There might be zero coeffs in newBcoeffs, but fmpz_mpoly_compose
            should have no problem with zero coeffs.
    */
    scales = (fmpq *) TMP_ALLOC(nvarsB*sizeof(fmpq));
    for (i = 0; i < nvarsB; i++)
    {
        Czpoly[i] = C[i]->zpoly;
        /* we are only borrowing the content of each of the C[i] */
        *(scales + i) = *C[i]->content;
    }
    *newB = *B->zpoly;
    newB->coeffs = _fmpz_vec_init(B->zpoly->length);
    _fmpq_mpoly_rescale(A->content, newB->coeffs, B, scales, ctxB);
/*
printf("newB: "); fmpz_mpoly_print_pretty(newB, NULL, ctxB->zctx); printf("\n");
*/
    fmpz_mpoly_compose_fmpz_mpoly(A->zpoly, newB, Czpoly, ctxB->zctx, ctxAC->zctx);
    fmpq_mpoly_canonicalise(A, ctxAC);

    _fmpz_vec_clear(newB->coeffs, B->zpoly->length);

    TMP_END;
}


