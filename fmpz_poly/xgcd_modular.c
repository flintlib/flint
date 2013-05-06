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
    Copyright (C) 2011 Sebastian Pancratz
   
******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "mpn_extras.h"

void _fmpz_poly_xgcd_modular(fmpz_t r, fmpz * s, fmpz * t, 
                             const fmpz * poly1, long len1, 
                             const fmpz * poly2, long len2)
{
    mp_ptr G, S, T, A, B, T1, T2;
    fmpz_t prod;
    int stabilised = 0, first;
    mp_limb_t p;
    mp_bitcnt_t s_bits = 0, t_bits = 0;

    /* Compute resultant of input polys */
    _fmpz_poly_resultant(r, poly1, len1, poly2, len2);

    if (fmpz_is_zero(r)) 
        return;

    fmpz_init(prod);
    fmpz_one(prod);

    _fmpz_vec_zero(s, len2);
    _fmpz_vec_zero(t, len1);

    p = (1UL << (FLINT_BITS - 1));

    G = _nmod_vec_init(4 * len1 + 5 * len2 - 2);
    S = G + len2;
    T = S + len2;
    A = T + len1;
    B = A + len1;
    T1 = B + len2;
    T2 = T1 + (len1 + len2 - 1);

    _nmod_vec_zero(S, len2 + len1); /* S = T = 0 */

    first = 1;

    for (;;) 
    {
        mp_limb_t R;
        nmod_t mod;

        /* Get next prime */
        p = n_nextprime(p, 0);

        /* Resultant mod p */
        R = fmpz_fdiv_ui(r, p);

        /* If p divides resultant or either leading coeff, discard p */
        if ((fmpz_fdiv_ui(poly1 + len1 - 1, p) == 0L) || 
            (fmpz_fdiv_ui(poly2 + len2 - 1, p) == 0L) || (R == 0))
            continue;

        nmod_init(&mod, p);

        /* Reduce polynomials modulo p */
        _fmpz_vec_get_nmod_vec(A, poly1, len1, mod);
        _fmpz_vec_get_nmod_vec(B, poly2, len2, mod);

        if (stabilised) /* CRT has stabilised, probably don't need more xgcds */
        {
            long tlen;

            /* Multiply out A*S + B*T to see if it is R mod p */
            _fmpz_vec_get_nmod_vec(S, s, len2, mod);
            _fmpz_vec_get_nmod_vec(T, t, len1, mod);

            _nmod_poly_mul(T1, A, len1, S, len2, mod); 
            _nmod_poly_mul(T2, T, len1, B, len2, mod);
            _nmod_vec_add(T1, T1, T2, len1 + len2 - 1, mod);
            tlen = len1 + len2 - 1;
            FMPZ_VEC_NORM(T1, tlen);

            if (tlen == 1 && T1[0] == R) /* It is, so this prime is good */
                fmpz_mul_ui(prod, prod, p);
            else
                stabilised = 0; /* It's not, keep going with xgcds */   
        }

        if (!stabilised) /* Need to keep computing xgcds mod p */
        {
            mp_limb_t RGinv;

            /* Compute xgcd mod p */
            _nmod_poly_xgcd(G, S, T, A, len1, B, len2, mod);
            RGinv = n_invmod(G[0], mod.n);
            RGinv = n_mulmod2_preinv(RGinv, R, mod.n, mod.ninv);

            /* Scale appropriately */
            _nmod_vec_scalar_mul_nmod(S, S, len2, RGinv, mod);
            _nmod_vec_scalar_mul_nmod(T, T, len1, RGinv, mod);

            if (first) /* First time around set s and t to S and T */
            {
                _fmpz_vec_set_nmod_vec(s, S, len2, mod);
                _fmpz_vec_set_nmod_vec(t, T, len1, mod);
                fmpz_set_ui(prod, p);

                stabilised = 1; /* Optimise the case where one prime is enough */
                first = 0;
            }
            else /* Otherwise do CRT */
            {
                mp_bitcnt_t new_s_bits, new_t_bits;

                _fmpz_poly_CRT_ui(s, s, len2, prod, S, len2, mod.n, mod.ninv, 1);
                _fmpz_poly_CRT_ui(t, t, len1, prod, T, len1, mod.n, mod.ninv, 1);
                fmpz_mul_ui(prod, prod, p);

                /* Check to see if CRT has stabilised */
                new_s_bits = FLINT_ABS(_fmpz_vec_max_bits(s, len2));
                new_t_bits = FLINT_ABS(_fmpz_vec_max_bits(t, len1));

                stabilised = (s_bits == new_s_bits && t_bits == new_t_bits);

                s_bits = new_s_bits;
                t_bits = new_t_bits;
            }
        }

        if (stabilised) 
        {
            long bound1, bound2, bound;

            bound1 = FLINT_BIT_COUNT(len2) 
                    + FLINT_ABS(_fmpz_vec_max_bits(poly1, len1)) 
                    + FLINT_ABS(_fmpz_vec_max_bits(s, len2));
            bound2 = FLINT_BIT_COUNT(len2) 
                    + FLINT_ABS(_fmpz_vec_max_bits(poly2, len2)) 
                    + FLINT_ABS(_fmpz_vec_max_bits(t, len1));

            bound = 4 + FLINT_MAX(fmpz_bits(r), FLINT_MAX(bound1, bound2));

            if (fmpz_bits(prod) > bound)
                break;
        }
    }

    _nmod_vec_clear(G);
    fmpz_clear(prod);
}

void
fmpz_poly_xgcd_modular(fmpz_t r, fmpz_poly_t s, fmpz_poly_t t,
                       const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    if (poly1->length < poly2->length)
    {
        fmpz_poly_xgcd_modular(r, t, s, poly2, poly1);
    } else /* len1 >= len2 >= 0 */
    {
        const long len1 = poly1->length;
        const long len2 = poly2->length;
        fmpz *S, *T;
        fmpz_poly_t temp1, temp2;

        if (len1 == 0 || len2 == 0) 
        {
            fmpz_zero(r);
        } 
        else /* len1 >= len2 >= 1 */
        {
            if (s == poly1 || s == poly2)
            {
                fmpz_poly_init2(temp1, len2);
                S = temp1->coeffs;
            }
            else
            {
                fmpz_poly_fit_length(s, len2);
                S = s->coeffs;
            }
    
            if (t == poly1 || t == poly2)
            {
                fmpz_poly_init2(temp2, len1);
                T = temp2->coeffs;
            }
            else
            {
                fmpz_poly_fit_length(t, len1);
                T = t->coeffs;
            }
     
            _fmpz_poly_xgcd_modular(r, S, T, poly1->coeffs, len1,
                                        poly2->coeffs, len2);
            
            if (s == poly1 || s == poly2)
            {
                fmpz_poly_swap(s, temp1);
                fmpz_poly_clear(temp1);
            }

            if (t == poly1 || t == poly2)
            {
                fmpz_poly_swap(t, temp2);
                fmpz_poly_clear(temp2);
            }
 
            _fmpz_poly_set_length(s, len2);
            _fmpz_poly_normalise(s);
 
            _fmpz_poly_set_length(t, len1);
            _fmpz_poly_normalise(t);
        }
    }
}

