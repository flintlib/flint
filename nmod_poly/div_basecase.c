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

    Copyright (C) 2011 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void
_nmod_poly_div_basecase_1(mp_ptr Q, mp_ptr W,
                             mp_srcptr A, len_t A_len, mp_srcptr B, len_t B_len,
                             nmod_t mod)
{
    mp_limb_t lead_inv = n_invmod(B[B_len - 1], mod.n);
    len_t len, coeff = A_len - B_len;
    
    mp_ptr R1 = W;
    mp_srcptr Btop = B + B_len - 1;
    
    flint_mpn_copyi(R1, A + B_len - 1, A_len - B_len + 1);

    while (coeff >= 0)
    {
        R1[coeff] = n_mod2_preinv(R1[coeff], mod.n, mod.ninv);

        while (coeff >= 0 && R1[coeff] == 0L)
        {
            Q[coeff--] = 0L;
            if (coeff >= 0)
                R1[coeff] = n_mod2_preinv(R1[coeff], mod.n, mod.ninv);
        }

        if (coeff >= 0)
        {
            mp_limb_t c, * R_sub;

            Q[coeff] =
                n_mulmod2_preinv(R1[coeff], lead_inv, mod.n, mod.ninv);

            c = n_negmod(Q[coeff], mod.n);

            len = FLINT_MIN(B_len - 1, coeff);
            R_sub = R1 + coeff - len;
            if (len > 0)
                mpn_addmul_1(R_sub, Btop - len, len, c);

            coeff--;
        }
    }
}

void
_nmod_poly_div_basecase_2(mp_ptr Q, mp_ptr W,
                             mp_srcptr A, len_t A_len, mp_srcptr B, len_t B_len,
                             nmod_t mod)
{
    len_t coeff, i, len;
    mp_limb_t lead_inv = n_invmod(B[B_len - 1], mod.n);
    mp_ptr B2, R2;
    mp_srcptr Btop;
    
    B2 = W;
    for (i = 0; i < B_len - 1; i++)
    {
        B2[2 * i] = B[i];
        B2[2 * i + 1] = 0;
    }
    Btop = B2 + 2*(B_len - 1);

    R2 = W + 2*(B_len - 1);
    for (i = 0; i < A_len - B_len + 1; i++)
    {
        R2[2 * i] = A[B_len + i - 1];
        R2[2 * i + 1] = 0;
    }

    coeff = A_len - B_len;
    
    while (coeff >= 0)
    {
        mp_limb_t r_coeff;
        r_coeff =
            n_ll_mod_preinv(R2[2 * coeff + 1], R2[2 * coeff], mod.n, mod.ninv);

        while (coeff >= 0 && r_coeff == 0L)
        {
            Q[coeff--] = 0L;
            if (coeff >= 0)
                r_coeff =
                    n_ll_mod_preinv(R2[2 * coeff + 1], R2[2 * coeff], mod.n,
                                    mod.ninv);
        }

        if (coeff >= 0)
        {
            mp_limb_t c, * R_sub;

            Q[coeff] =
                n_mulmod2_preinv(r_coeff, lead_inv, mod.n, mod.ninv);

            c = n_negmod(Q[coeff], mod.n);

            len = FLINT_MIN(B_len - 1, coeff);
            R_sub = R2 + 2 * (coeff - len);
            if (len > 0)
                mpn_addmul_1(R_sub, Btop - 2*len, 2 * len, c);

            coeff--;
        }
    }
}

void
_nmod_poly_div_basecase_3(mp_ptr Q, mp_ptr W,
                             mp_srcptr A, len_t A_len, mp_srcptr B, len_t B_len,
                             nmod_t mod)
{
    len_t coeff, i, len;
    mp_limb_t lead_inv = n_invmod(B[B_len - 1], mod.n);
    mp_limb_t r_coeff;
    mp_ptr B3, R3;
    mp_srcptr Btop;
    
    B3 = W;
    for (i = 0; i < B_len - 1; i++)
    {
        B3[3 * i] = B[i];
        B3[3 * i + 1] = 0;
        B3[3 * i + 2] = 0;
    }
    Btop = B3 + 3*(B_len - 1);

    R3 = W + 3*(B_len - 1);
    for (i = 0; i < A_len - B_len + 1; i++)
    {
        R3[3 * i] = A[B_len + i - 1];
        R3[3 * i + 1] = 0;
        R3[3 * i + 2] = 0;
    }

    coeff = A_len - B_len;
    
    while (coeff >= 0)
    {
        r_coeff =
            n_lll_mod_preinv(R3[3 * coeff + 2], R3[3 * coeff + 1],
                             R3[3 * coeff], mod.n, mod.ninv);

        while (coeff >= 0 && r_coeff == 0L)
        {
            Q[coeff--] = 0L;
            if (coeff >= 0)
                r_coeff =
                    n_lll_mod_preinv(R3[3 * coeff + 2], R3[3 * coeff + 1],
                                     R3[3 * coeff], mod.n, mod.ninv);
        }

        if (coeff >= 0)
        {
            mp_limb_t c, * R_sub;

            Q[coeff] =
                n_mulmod2_preinv(r_coeff, lead_inv, mod.n, mod.ninv);

            c = n_negmod(Q[coeff], mod.n);

            len = FLINT_MIN(B_len - 1, coeff);
            R_sub = R3 + 3 * (coeff - len);
            if (len > 0)
                mpn_addmul_1(R_sub, Btop - 3*len, 3 * len, c);

            coeff--;
        }
    }
}

void
_nmod_poly_div_basecase(mp_ptr Q, mp_ptr W,
                           mp_srcptr A, len_t A_len, mp_srcptr B, len_t B_len,
                           nmod_t mod)
{
    len_t bits =
        2 * (FLINT_BITS - mod.norm) + FLINT_BIT_COUNT(A_len - B_len + 1);

    if (bits <= FLINT_BITS)
        _nmod_poly_div_basecase_1(Q, W, A, A_len, B, B_len, mod);
    else if (bits <= 2 * FLINT_BITS)
        _nmod_poly_div_basecase_2(Q, W, A, A_len, B, B_len, mod);
    else
        _nmod_poly_div_basecase_3(Q, W, A, A_len, B, B_len, mod);
}

void
nmod_poly_div_basecase(nmod_poly_t Q, const nmod_poly_t A,
                          const nmod_poly_t B)
{
    mp_ptr Q_coeffs, W;
    nmod_poly_t t1;
    len_t Alen, Blen;

    Blen = B->length;

    if (Blen == 0)
    {
        printf("Exception (nmod_poly_div_base). Division by zero.\n");
        abort();
    }

    Alen = A->length;

    if (Alen < Blen)
    {
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

    W = _nmod_vec_init(NMOD_DIV_BC_ITCH(Alen, Blen, A->mod));
    
    _nmod_poly_div_basecase(Q_coeffs, W, A->coeffs, Alen,
                               B->coeffs, Blen, B->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(Q, t1);
        nmod_poly_clear(t1);
    }
    
    Q->length = Alen - Blen + 1;

    _nmod_vec_clear(W);
    _nmod_poly_normalise(Q);
}
