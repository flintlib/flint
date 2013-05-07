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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_radix_init(fmpz **Rpow, fmpz **Rinv, 
                    const fmpz *R, long lenR, long k, 
                    const fmpz_t invL, const fmpz_t p)
{
    const long degR = lenR - 1;
    long i;
    fmpz_t invLP;
    fmpz *W;

    fmpz_init_set(invLP, invL);
    W = flint_malloc((1L << (k - 1)) * degR * sizeof(fmpz));

    _fmpz_vec_set(Rpow[0], R, lenR);
    for (i = 1; i < k; i++)
    {
        _fmpz_mod_poly_sqr(Rpow[i], Rpow[i - 1], degR * (1L << (i - 1)) + 1, p);
    }

    for (i = 0; i < k; i++)
    {
        const long lenQ = (1L << i) * degR;
        long j;

        /* W := rev{Rpow[i], lenQ} */
        for (j = 0; j < lenQ; j++)
        {
            W[j] = Rpow[i][lenQ - j];
        }

        _fmpz_mod_poly_inv_series_newton(Rinv[i], W, lenQ, invLP, p);

        /* invLP := inv{lead{R^{2^i}}} */
        if (i != k - 1)
        {
            fmpz_mul(invLP, invLP, invLP);
            fmpz_mod(invLP, invLP, p);
        }
    }

    fmpz_clear(invLP);
    flint_free(W);
}

void fmpz_mod_poly_radix_init(fmpz_mod_poly_radix_t D, 
                              const fmpz_mod_poly_t R, long degF)
{
    const long degR = R->length - 1;

    if (degF < degR)
    {
        D->k = 0;
        D->degR = degR;
    }
    else
    {
        const long N = degF / degR;
        const long k = FLINT_BIT_COUNT(N);     /* k := ceil{log{N+1}} */
        const long lenV = degR * ((1L << k) - 1) + k;
        const long lenW = degR * ((1L << k) - 1);

        long i;

        D->V = _fmpz_vec_init(lenV + lenW);
        D->W = D->V + lenV;

        D->Rpow = flint_malloc(k * sizeof(fmpz *));
        D->Rinv = flint_malloc(k * sizeof(fmpz *));

        for (i = 0; i < k; i++)
        {
            D->Rpow[i] = D->V + (degR * ((1L << i) - 1) + i);
            D->Rinv[i] = D->W + (degR * ((1L << i) - 1));
        }

        fmpz_init(&(D->invL));
        fmpz_invmod(&(D->invL), R->coeffs + degR, &(R->p));

        _fmpz_mod_poly_radix_init(D->Rpow, D->Rinv, R->coeffs, degR + 1, 
                                  k, &(D->invL), &(R->p));

        D->k = k;
        D->degR = degR;
    }
}

void fmpz_mod_poly_radix_clear(fmpz_mod_poly_radix_t D)
{
    if (D->k)
    {
        const long degR = D->degR;
        const long k    = D->k;
        const long lenV = degR * ((1L << k) - 1) + k;
        const long lenW = degR * ((1L << k) - 1);

        _fmpz_vec_clear(D->V, lenV + lenW);
        flint_free(D->Rpow);
        flint_free(D->Rinv);
        fmpz_clear(&(D->invL));
    }
}

void _fmpz_mod_poly_radix(fmpz **B, const fmpz *F, fmpz **Rpow, fmpz **Rinv, 
                          long degR, long k, long i, fmpz *W, const fmpz_t p)
{
    if (i == -1)
    {
        _fmpz_vec_set(B[k], F, degR);
    }
    else
    {
        const long lenQ = (1L << i) * degR;

        fmpz *Frev = W;
        fmpz *Q    = W + lenQ;
        fmpz *S    = W;

        _fmpz_poly_reverse(Frev, F + lenQ, lenQ, lenQ);
        _fmpz_mod_poly_mullow(Q, Frev, lenQ, Rinv[i], lenQ, p, lenQ);
        _fmpz_poly_reverse(Q, Q, lenQ, lenQ);

        _fmpz_mod_poly_radix(B, Q, Rpow, Rinv, degR, k + (1L << i), i-1, W, p);

        _fmpz_mod_poly_mullow(S, Rpow[i], lenQ, Q, lenQ, p, lenQ);
        _fmpz_mod_poly_sub(S, F, lenQ, S, lenQ, p);

        _fmpz_mod_poly_radix(B, S, Rpow, Rinv, degR, k, i-1, W + lenQ, p);
    }
}

void fmpz_mod_poly_radix(fmpz_mod_poly_struct **B, 
                         const fmpz_mod_poly_t F, 
                         const fmpz_mod_poly_radix_t D)
{
    const long lenF = F->length;
    const long degF = F->length - 1;
    const long degR = D->degR;
    const long N    = degF / degR;

    if (N == 0)
    {
        fmpz_mod_poly_set(B[0], F);
    }
    else
    {
        const long k    = FLINT_BIT_COUNT(N);     /* k := ceil{log{N+1}}    */
        const long lenG = (1L << k) * degR;       /* Padded size            */
        const long t    = (lenG - 1) / degR - N;  /* Extra {degR}-blocks    */

        fmpz *G;                                  /* Padded copy of F       */
        fmpz *T;                                  /* Additional B[i]        */
        fmpz **C;                                 /* Enlarged version of B  */
        fmpz *W;                                  /* Temporary space        */

        long i;

        if (lenF < lenG)
        {
            G = flint_malloc(lenG * sizeof(fmpz));
            for (i = 0; i < lenF; i++)
                G[i] = F->coeffs[i];
            flint_mpn_zero((mp_ptr) G + lenF, lenG - lenF);

            T = t ? _fmpz_vec_init(t * degR) : NULL;
        }
        else  /* lenF == lenG */
        {
            G = F->coeffs;
            T = NULL;
        }

        C = flint_malloc((N + 1 + t) * sizeof(fmpz *));
        for (i = 0; i <= N; i++)
        {
            fmpz_mod_poly_fit_length(B[i], degR);
            C[i] = B[i]->coeffs;
        }
        for (i = 0; i < t; i++)
        {
            C[N + 1 + i] = T + i * degR;
        }

        W = _fmpz_vec_init(lenG);

        _fmpz_mod_poly_radix(C, G, D->Rpow, D->Rinv, degR, 0, k-1, W, &(F->p));

        _fmpz_vec_clear(W, lenG);

        for (i = 0; i <= N; i++)
        {
            _fmpz_mod_poly_set_length(B[i], degR);
            _fmpz_mod_poly_normalise(B[i]);
        }
        flint_free(C);
        if (lenF < lenG)
        {
            flint_free(G);
        }
        if (t) 
        {
            _fmpz_vec_clear(T, t * degR);
        }
    }
}

