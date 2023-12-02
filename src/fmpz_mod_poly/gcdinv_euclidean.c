/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

slong _fmpz_mod_poly_gcdinv_euclidean(fmpz *G, fmpz *S,
                                   const fmpz *A, slong lenA,
                                   const fmpz *B, slong lenB,
                                   const fmpz_t invA, const fmpz_mod_ctx_t ctx)
{
	_fmpz_vec_zero(G, lenA);
    _fmpz_vec_zero(S, lenB - 1);

    if (lenA == 1)
    {
        fmpz_set(G + 0, A + 0);
        fmpz_one(S + 0);
        return 1;
    }
    else
    {
        fmpz *Q, *R;
        slong i, lenQ, lenR;
	    TMP_INIT;
	    TMP_START;

        FMPZ_VEC_TMP_INIT(Q, 2*lenB);
        R = Q + lenB;

        _fmpz_mod_poly_divrem(Q, R, B, lenB, A, lenA, invA, ctx);
        lenR = lenA - 1;
        FMPZ_VEC_NORM(R, lenR);

        if (lenR == 0)
        {
            _fmpz_vec_set(G, A, lenA);
            fmpz_one(S + 0);
            FMPZ_VEC_TMP_CLEAR(Q, 2*lenB);
			TMP_END;
            return lenA;
        }
        else if (lenR == 1)
		{
			lenQ = lenB - lenA + 1;
			FMPZ_VEC_NORM(Q, lenQ);

			_fmpz_vec_swap(G, R, lenR);
            _fmpz_vec_swap(S, Q, lenQ);
		    _fmpz_mod_vec_neg(S, S, lenQ, ctx);

            FMPZ_VEC_TMP_CLEAR(Q, 2*lenB);
			TMP_END;

			return 1;
		}
        else
        {
            fmpz_t inv;
            fmpz *D, *U1, *U2, *V3, *W;
            slong lenD, lenU1, lenU2, lenV3, lenW;

            fmpz_init(inv);
            FMPZ_VEC_TMP_INIT(W, 3*lenB + 2*lenA);
            D  = W  + lenB;
            U1 = D  + lenA;
            U2 = U1 + lenB;
            V3 = U2 + lenB;

			lenQ = lenB - lenA + 1;
			FMPZ_VEC_NORM(Q, lenQ);

            fmpz_set_ui(U1, 1);
			lenU1 = 1;
            _fmpz_vec_set(D, A, lenA);
            lenD = lenA;
            _fmpz_mod_vec_neg(U2, Q, lenQ, ctx);
            lenU2 = lenQ;
            lenV3 = 0;
            FMPZ_VEC_SWAP(V3, lenV3, R, lenR);

            do {
                fmpz_invmod(inv, V3 + (lenV3 - 1), fmpz_mod_ctx_modulus(ctx));
                _fmpz_mod_poly_divrem_basecase(Q, D, D, lenD, V3, lenV3, inv, ctx);
                lenQ = lenD - lenV3 + 1;
                lenD = lenV3 - 1;
                FMPZ_VEC_NORM(D, lenD);

                if (lenV3 != 0)
				{
				    if (lenU2 >= lenQ)
                        _fmpz_mod_poly_mul(W, U2, lenU2, Q, lenQ, ctx);
                    else
                        _fmpz_mod_poly_mul(W, Q, lenQ, U2, lenU2, ctx);
                    lenW = lenQ + lenU2 - 1;

                    _fmpz_mod_poly_sub(U1, U1, lenU1, W, lenW, ctx);
                    lenU1 = FLINT_MAX(lenU1, lenW);
                    FMPZ_VEC_NORM(U1, lenU1);
                }

                FMPZ_VEC_SWAP(U1, lenU1, U2, lenU2);
				FMPZ_VEC_SWAP(V3, lenV3, D, lenD);
            } while (lenV3 != 0);

            _fmpz_vec_swap(G, D, lenD);
            _fmpz_vec_swap(S, U1, lenU1);

            for (i = 0; i < lenU1; i++)
            {
                if (fmpz_sgn(S + i) < 0)
                    fmpz_add(S + i, S + i, fmpz_mod_ctx_modulus(ctx));
            }

            FMPZ_VEC_TMP_CLEAR(W, 3*lenB + 2*lenA);
            FMPZ_VEC_TMP_CLEAR(Q, 2*lenB);

		    fmpz_clear(inv);

            TMP_END;

			return lenD;
        }
	}
}

void fmpz_mod_poly_gcdinv_euclidean(fmpz_mod_poly_t G, fmpz_mod_poly_t S,
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                      const fmpz_mod_ctx_t ctx)
{
    const slong lenA = A->length, lenB = B->length;
    fmpz_t inv;

	if (lenB < 2)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_gcdinv). lenB < 2.\n");
    }
    if (lenA >= lenB)
    {
        fmpz_mod_poly_t T;

        fmpz_mod_poly_init(T, ctx);
        fmpz_mod_poly_rem(T, A, B, ctx);
        fmpz_mod_poly_gcdinv_euclidean(G, S, T, B, ctx);
        fmpz_mod_poly_clear(T, ctx);
        return;
    }

    fmpz_init(inv);
    if (lenA == 0)
	{
        fmpz_mod_poly_zero(G, ctx);
        fmpz_mod_poly_zero(S, ctx);
    }
    else  /* lenB >= lenA >= 1 */
    {
        fmpz *g, *s;
        slong lenG;

        if (G == A || G == B)
        {
            g = _fmpz_vec_init(FLINT_MIN(lenA, lenB));
        }
        else
        {
            fmpz_mod_poly_fit_length(G, FLINT_MIN(lenA, lenB), ctx);
            g = G->coeffs;
        }
        if (S == A || S == B)
        {
            s = _fmpz_vec_init(lenB);
        }
        else
        {
            fmpz_mod_poly_fit_length(S, lenB, ctx);
            s = S->coeffs;
        }

        fmpz_invmod(inv, fmpz_mod_poly_lead(A, ctx), fmpz_mod_ctx_modulus(ctx));
        lenG = _fmpz_mod_poly_gcdinv_euclidean(g, s, A->coeffs, lenA,
                              B->coeffs, lenB, inv, ctx);

        if (G == A || G == B)
        {
            _fmpz_vec_clear(G->coeffs, G->alloc);
            G->coeffs = g;
            G->alloc  = FLINT_MIN(lenA, lenB);
        }
        if (S == A || S == B)
        {
            _fmpz_vec_clear(S->coeffs, S->alloc);
            S->coeffs = s;
            S->alloc  = lenB;
        }

        _fmpz_mod_poly_set_length(G, lenG);
        _fmpz_mod_poly_set_length(S, FLINT_MAX(lenB - lenG, 1));
        _fmpz_mod_poly_normalise(S);

        if (!fmpz_is_one(fmpz_mod_poly_lead(G, ctx)))
        {
            fmpz_invmod(inv, fmpz_mod_poly_lead(G, ctx),
                                                    fmpz_mod_ctx_modulus(ctx));
            fmpz_mod_poly_scalar_mul_fmpz(G, G, inv, ctx);
            fmpz_mod_poly_scalar_mul_fmpz(S, S, inv, ctx);
        }
        fmpz_clear(inv);
    }
}

