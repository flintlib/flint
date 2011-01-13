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

    Copyright (C) 2010, 2011 William Hart

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void
_nmod_poly_divrem_basecase_1(mp_ptr Q, mp_ptr R, mp_ptr W,
                             mp_srcptr A, long A_len, mp_srcptr B, long B_len,
                             nmod_t mod)
{
    mp_limb_t lead_inv = n_invmod(B[B_len - 1], mod.n);
    long coeff = A_len - 1;
    mp_ptr coeff_Q = Q - B_len + 1;

    mp_ptr R1 = W;
    
    mpn_copyi(R1, A, A_len);

    while (coeff + 1 >= B_len)
    {
        R1[coeff] = n_mod2_preinv(R1[coeff], mod.n, mod.ninv);

        while ((coeff + 1 >= B_len) && (R1[coeff] == 0L))
        {
            coeff_Q[coeff--] = 0L;
            if (coeff + 1 >= B_len)
                R1[coeff] = n_mod2_preinv(R1[coeff], mod.n, mod.ninv);
        }

        if (coeff + 1 >= B_len)
        {
            mp_limb_t c, * R_sub;

            coeff_Q[coeff] =
                n_mulmod2_preinv(R1[coeff], lead_inv, mod.n, mod.ninv);

            c = n_negmod(coeff_Q[coeff], mod.n);

            R_sub = R1 + coeff - B_len + 1;
            if (B_len > 1)
                mpn_addmul_1(R_sub, B, B_len - 1, c);

            coeff--;
        }
    }

    while (coeff + 1 > 0)
    {
        R[coeff] = n_mod2_preinv(R1[coeff], mod.n, mod.ninv);
        coeff--;
    }
}

void
_nmod_poly_divrem_basecase_2(mp_ptr Q, mp_ptr R, mp_ptr W,
                             mp_srcptr A, long A_len, mp_srcptr B, long B_len,
                             nmod_t mod)
{
    long coeff, i;
    mp_limb_t lead_inv = n_invmod(B[B_len - 1], mod.n);
    mp_ptr coeff_Q, B2, R2;

    B2 = W;
    for (i = 0; i < B_len - 1; i++)
    {
        B2[2 * i] = B[i];
        B2[2 * i + 1] = 0;
    }

    R2 = W + 2*(B_len - 1);
    for (i = 0; i < A_len; i++)
    {
        R2[2 * i] = A[i];
        R2[2 * i + 1] = 0;
    }

    coeff = A_len - 1;
    coeff_Q = Q - B_len + 1;

    while (coeff + 1 >= B_len)
    {
        mp_limb_t r_coeff;
        r_coeff =
            n_ll_mod_preinv(R2[2 * coeff + 1], R2[2 * coeff], mod.n, mod.ninv);

        while ((coeff + 1 >= B_len) && (r_coeff == 0L))
        {
            coeff_Q[coeff--] = 0L;
            if (coeff + 1 >= B_len)
                r_coeff =
                    n_ll_mod_preinv(R2[2 * coeff + 1], R2[2 * coeff], mod.n,
                                    mod.ninv);
        }

        if (coeff + 1 >= B_len)
        {
            mp_limb_t c, * R_sub;

            coeff_Q[coeff] =
                n_mulmod2_preinv(r_coeff, lead_inv, mod.n, mod.ninv);

            c = n_negmod(coeff_Q[coeff], mod.n);

            R_sub = R2 + 2 * (coeff - B_len + 1);
            if (B_len > 1)
                mpn_addmul_1(R_sub, B2, 2 * B_len - 2, c);

            coeff--;
        }
    }

    while (coeff + 1 > 0)
    {
        R[coeff] =
            n_ll_mod_preinv(R2[2 * coeff + 1], R2[2 * coeff], mod.n, mod.ninv);
        coeff--;
    }
}

void
_nmod_poly_divrem_basecase_3(mp_ptr Q, mp_ptr R, mp_ptr W,
                             mp_srcptr A, long A_len, mp_srcptr B, long B_len,
                             nmod_t mod)
{
    long coeff, i;
    mp_limb_t lead_inv = n_invmod(B[B_len - 1], mod.n);
    mp_limb_t r_coeff;
    mp_ptr B3, R3, coeff_Q;

    B3 = W;
    for (i = 0; i < B_len - 1; i++)
    {
        B3[3 * i] = B[i];
        B3[3 * i + 1] = 0;
        B3[3 * i + 2] = 0;
    }

    R3 = W + 3*(B_len - 1);
    for (i = 0; i < A_len; i++)
    {
        R3[3 * i] = A[i];
        R3[3 * i + 1] = 0;
        R3[3 * i + 2] = 0;
    }

    coeff = A_len - 1;
    coeff_Q = Q - B_len + 1;

    while (coeff + 1 >= B_len)
    {
        r_coeff =
            n_lll_mod_preinv(R3[3 * coeff + 2], R3[3 * coeff + 1],
                             R3[3 * coeff], mod.n, mod.ninv);

        while ((coeff + 1 >= B_len) && (r_coeff == 0L))
        {
            coeff_Q[coeff--] = 0L;
            if (coeff + 1 >= B_len)
                r_coeff =
                    n_lll_mod_preinv(R3[3 * coeff + 2], R3[3 * coeff + 1],
                                     R3[3 * coeff], mod.n, mod.ninv);
        }

        if (coeff + 1 >= B_len)
        {
            mp_limb_t c, * R_sub;

            coeff_Q[coeff] =
                n_mulmod2_preinv(r_coeff, lead_inv, mod.n, mod.ninv);

            c = n_negmod(coeff_Q[coeff], mod.n);

            R_sub = R3 + 3 * (coeff - B_len + 1);
            if (B_len > 1)
                mpn_addmul_1(R_sub, B3, 3 * B_len - 3, c);

            coeff--;
        }
    }

    while (coeff + 1 > 0)
    {
        R[coeff] =
            n_lll_mod_preinv(R3[3 * coeff + 2], R3[3 * coeff + 1],
                             R3[3 * coeff], mod.n, mod.ninv);
        coeff--;
    }
}

void
_nmod_poly_divrem_basecase(mp_ptr Q, mp_ptr R, mp_ptr W,
                           mp_srcptr A, long A_len, mp_srcptr B, long B_len,
                           nmod_t mod)
{
    long bits =
        2 * (FLINT_BITS - mod.norm) + FLINT_BIT_COUNT(A_len - B_len + 1);

    if (bits <= FLINT_BITS)
        _nmod_poly_divrem_basecase_1(Q, R, W, A, A_len, B, B_len, mod);
    else if (bits <= 2 * FLINT_BITS)
        _nmod_poly_divrem_basecase_2(Q, R, W, A, A_len, B, B_len, mod);
    else
        _nmod_poly_divrem_basecase_3(Q, R, W, A, A_len, B, B_len, mod);
}

void
nmod_poly_divrem_basecase(nmod_poly_t Q, nmod_poly_t R, const nmod_poly_t A,
                          const nmod_poly_t B)
{
    mp_ptr Q_coeffs, R_coeffs, W;
    nmod_poly_t t1, t2;
    long Alen, Blen;

    Blen = B->length;

    if (Blen == 0)
    {
        printf("Exception: division by zero in nmod_poly_divrem_basecase\n");
        abort();
    }

    Alen = A->length;

    if (Alen < Blen)
    {
        nmod_poly_set(R, A);
        nmod_poly_zero(Q);

        return;
    }

    if (Q == A || Q == B)
    {
        nmod_poly_init2_preinv(t1, B->mod.n, B->mod.ninv,
                               Alen - Blen + 1);
        Q_coeffs = t1->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Q, Alen - Blen + 1);
        Q_coeffs = Q->coeffs;
    }

    if (R == A || R == B)
    {
        nmod_poly_init2_preinv(t2, B->mod.n, B->mod.ninv, Blen - 1);
        R_coeffs = t2->coeffs;
    }
    else
    {
        nmod_poly_fit_length(R, Blen - 1);
        R_coeffs = R->coeffs;
    }

    W = _nmod_vec_init(NMOD_DIVREM_BC_ITCH(Alen, Blen, A->mod));
    
    _nmod_poly_divrem_basecase(Q_coeffs, R_coeffs, W, A->coeffs, Alen,
                               B->coeffs, Blen, B->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(Q, t1);
        nmod_poly_clear(t1);
    }
        
    Q->length = Alen - Blen + 1;

    if (R == A || R == B)
    {
        nmod_poly_swap(R, t2);
        nmod_poly_clear(t2);
    }
        
    R->length = Blen - 1;

    _nmod_vec_free(W);
    _nmod_poly_normalise(Q);
    _nmod_poly_normalise(R);
}
