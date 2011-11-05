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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_change_radix_precomp(fmpz **Rpow, fmpz **Rinv, 
                    const fmpz *R, long lenR, long k, 
                    const fmpz_t invL, const fmpz_t p)
{
    const long degR = lenR - 1;
    long i;
    fmpz_t invLP;
    fmpz *W;

    fmpz_init_set(invLP, invL);
    W = malloc((1L << (k - 1)) * degR * sizeof(fmpz));

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
    free(W);
}

void _fmpz_mod_poly_change_radix_recursive(fmpz **B, 
                    const fmpz *F, fmpz **Rpow, fmpz **Rinv, long degR, 
                    long k, long i, fmpz *W, const fmpz_t p)
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

        _fmpz_mod_poly_change_radix_recursive(B, Q, Rpow, Rinv, degR, 
                                              k + (1L << i), i-1, W, p);

        _fmpz_mod_poly_mullow(S, Rpow[i], lenQ, Q, lenQ, p, lenQ);
        _fmpz_mod_poly_sub(S, F, lenQ, S, lenQ, p);

        _fmpz_mod_poly_change_radix_recursive(B, S, Rpow, Rinv, degR, 
                                              k, i-1, W + lenQ, p);
    }
}

void _fmpz_mod_poly_change_radix(fmpz **B, const fmpz *F, 
                                 const fmpz *R, long lenR, long k, 
                                 const fmpz_t invL, const fmpz_t p)
{
    const long degR = lenR - 1;
    const long lenV = degR * ((1L << k) - 1) + k;
    const long lenW = degR * ((1L << k) - 1);
    const long lenX = degR * (1L << k);

    long i;
    fmpz **Rpow, **Rinv;
    fmpz *V, *W, *X;

    V = _fmpz_vec_init(lenV + lenW + lenX);
    W = V + lenV;
    X = W + lenW;

    Rpow = malloc(k * sizeof(fmpz *));
    Rinv = malloc(k * sizeof(fmpz *));

    for (i = 0; i < k; i++)
    {
        Rpow[i] = V + (degR * ((1L << i) - 1) + i);
        Rinv[i] = W + (degR * ((1L << i) - 1));
    }

    _fmpz_mod_poly_change_radix_precomp(Rpow, Rinv, R, lenR, k, invL, p);

    _fmpz_mod_poly_change_radix_recursive(B, F, Rpow, Rinv, degR, 0, k-1, X, p);

    free(Rpow);
    free(Rinv);

    _fmpz_vec_clear(V, lenV + lenW + lenX);
}

void fmpz_mod_poly_change_radix(fmpz_mod_poly_struct **B, 
                                const fmpz_mod_poly_t F, 
                                const fmpz_mod_poly_t R)
{
    const long lenF = F->length;
    const long degF = F->length - 1;
    const long degR = R->length - 1;
    const long N    = degF / degR;

    if (N == 0)
    {
        fmpz_mod_poly_set(B[0], F);
    }
    else
    {
        const long k    = FLINT_BIT_COUNT(N);     /* k := ceil{log{N+1}} */
        const long lenG = (1L << k) * degR;       /* Padded size */
        const long t    = (lenG - 1) / degR - N;  /* Extra {degR}-blocks */

        fmpz *G, *T, **C;
        fmpz_t invL;
        long i;

        if (lenF < lenG)
        {
            G = malloc(lenG * sizeof(fmpz));
            for (i = 0; i < lenF; i++)
                G[i] = F->coeffs[i];
            mpn_zero((mp_ptr) G + lenF, lenG - lenF);

            T = t ? _fmpz_vec_init(t * degR) : NULL;
        }
        else  /* lenF == lenG */
        {
            G = F->coeffs;
            T = NULL;
        }

        C = malloc((N + 1 + t) * sizeof(fmpz *));
        for (i = 0; i <= N; i++)
        {
            fmpz_mod_poly_fit_length(B[i], degR);
            C[i] = B[i]->coeffs;
        }
        for (i = 0; i < t; i++)
        {
            C[N + 1 + i] = T + i * degR;
        }

        fmpz_init(invL);
        fmpz_invmod(invL, R->coeffs + degR, &(R->p));

        _fmpz_mod_poly_change_radix(C, G, R->coeffs, R->length, k, invL, &(R->p));

        fmpz_clear(invL);

        for (i = 0; i <= N; i++)
        {
            _fmpz_mod_poly_set_length(B[i], degR);
            _fmpz_mod_poly_normalise(B[i]);
        }
        free(C);
        if (lenF < lenG)
        {
            free(G);
        }
        if (t) 
        {
            _fmpz_vec_clear(T, t * degR);
        }
    }
}

