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

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "mpn_extras.h"


void _fmpz_poly_gcd_modular(fmpz * res, const fmpz * poly1, len_t len1, 
                                        const fmpz * poly2, len_t len2)
{
    mp_bitcnt_t bits1, bits2, nb1, nb2, bits_small, pbits, curr_bits = 0, new_bits;   
    fmpz_t ac, bc, hc, d, g, l, eval_A, eval_B, eval_GCD, modulus;
    fmpz * A, * B, * Q, * lead_A, * lead_B;
    mp_ptr a, b, h;
    mp_limb_t p, h_inv, g_mod;
    nmod_t mod;
    len_t i, n, n0, unlucky, hlen, bound;
    int g_pm1;

    fmpz_init(ac);
    fmpz_init(bc);
    fmpz_init(d);

    /* compute gcd of content of poly1 and poly2 */
    _fmpz_vec_content(ac, poly1, len1);
    _fmpz_vec_content(bc, poly2, len2);
    fmpz_gcd(d, ac, bc);

    /* special case, one of the polys is a constant */
    if (len2 == 1) /* if len1 == 1 then so does len2 */
    {
        fmpz_set(res, d);

        fmpz_clear(ac);
        fmpz_clear(bc);
        fmpz_clear(d);
        return;
    }

    /* divide poly1 and poly2 by their content */
    A = _fmpz_vec_init(len1);
    B = _fmpz_vec_init(len2);
    _fmpz_vec_scalar_divexact_fmpz(A, poly1, len1, ac);
    _fmpz_vec_scalar_divexact_fmpz(B, poly2, len2, bc);
    fmpz_clear(ac);
    fmpz_clear(bc);

    /* get bound on size of gcd coefficients */
    lead_A = A + len1 - 1;
    lead_B = B + len2 - 1;

    bits1 = _fmpz_vec_max_bits(A, len1); bits1 = FLINT_ABS(bits1);
    bits2 = _fmpz_vec_max_bits(B, len2); bits2 = FLINT_ABS(bits2);

    fmpz_init(l);
   
    if (len1 < 64 && len2 < 64) /* compute the squares of the 2-norms */
    {
        fmpz_set_ui(l, 0);
        for (i = 0; i < len1; i++)
            fmpz_addmul(l, A + i, A + i);
        nb1 = fmpz_bits(l);
        fmpz_set_ui(l, 0);
        for (i = 0; i < len2; i++)
            fmpz_addmul(l, B + i, B + i);
        nb2 = fmpz_bits(l);
    } else /* approximate to save time */
    {
        nb1 = 2*bits1 + FLINT_BIT_COUNT(len1);
        nb2 = 2*bits2 + FLINT_BIT_COUNT(len2);
    }

    /* get gcd of leading coefficients */
    fmpz_init(g);
    fmpz_gcd(g, lead_A, lead_B);
    fmpz_mul(l, lead_A, lead_B);

    g_pm1 = fmpz_is_pm1(g);
   
    /* evaluate -A at -1 */
    fmpz_init(eval_A);
    for (i = 0; i < len1; i++)
    {
        if (i & 1) fmpz_add(eval_A, eval_A, A + i);
        else fmpz_sub(eval_A, eval_A, A + i);
    }

    /* evaluate -B at -1 */
    fmpz_init(eval_B);
    for (i = 0; i < len2; i++)
    {
        if (i & 1) fmpz_add(eval_B, eval_B, B + i);
        else fmpz_sub(eval_B, eval_B, B + i);
    }

    /* compute the gcd of eval(-A, -1) and eval(-B, -1) */
    fmpz_init(eval_GCD);
    fmpz_gcd(eval_GCD, eval_A, eval_B);

    /* compute a heuristic bound after which we should begin checking if we're done */
    bits_small = FLINT_MAX(fmpz_bits(eval_GCD), fmpz_bits(g));
    if (bits_small < 2L) bits_small = 2;

    fmpz_clear(eval_GCD);
    fmpz_clear(eval_A);
    fmpz_clear(eval_B);

    /* set size of first prime */
    pbits = FLINT_BITS - 1;
    p = (1UL<<pbits);

    fmpz_init(modulus);
    fmpz_init(hc);

    Q = _fmpz_vec_init(len1);

    /* make space for polynomials mod p */
    a = _nmod_vec_init(len1);
    b = _nmod_vec_init(len2);
    h = _nmod_vec_init(len2);

    /* zero entire output */
    _fmpz_vec_zero(res, len2);

    n = len2; 
    /* current bound on length of result 
      the bound we use is from section 6 of 
       http://cs.nyu.edu/~yap/book/alge/ftpSite/l4.ps.gz 
    */
    n0 = len1 - 1;
    bound = (n0 + 3)*FLINT_MAX(nb1, nb2) + (n0 + 1); /* initialise bound */
    unlucky = 0;

    for (;;)
    {
        /* get new prime and initialise modulus */
        p = n_nextprime(p, 0);
        if (fmpz_fdiv_ui(l, p) == 0)
        {
            unlucky += pbits;
            continue;
        }
        nmod_init(&mod, p);

        /* reduce polynomials modulo p */
        _fmpz_vec_get_nmod_vec(a, A, len1, mod);
        _fmpz_vec_get_nmod_vec(b, B, len2, mod);

        /* compute gcd over Z/pZ */
        hlen = _nmod_poly_gcd(h, a, len1, b, len2, mod);

        if (hlen == 1) /* gcd is 1 */
        {
            fmpz_one(res);
            _fmpz_vec_zero(res + 1, len2 - 1);
            break; 
        }

        if (hlen > n + 1) /* discard this prime */
        {
            unlucky += pbits;
            continue;
        }

        /* scale new polynomial mod p appropriately */
        if (g_pm1) _nmod_poly_make_monic(h, h, hlen, mod);
        else
        {
            h_inv = n_invmod(h[hlen - 1], mod.n);
            g_mod = fmpz_fdiv_ui(g, mod.n);
            h_inv = n_mulmod2_preinv(h_inv, g_mod, mod.n, mod.ninv);
            _nmod_vec_scalar_mul_nmod(h, h, hlen, h_inv, mod);
        }

        if (hlen <= n) /* we have a new bound on size of result */
        {
            unlucky += fmpz_bits(modulus);
            _fmpz_vec_set_nmod_vec(res, h, hlen, mod);
            _fmpz_vec_zero(res + hlen, len2 - hlen);

            if (g_pm1)
            {
                /* are we done? */
                if (_fmpz_poly_divides(Q, B, len2, res, hlen) &&
                    _fmpz_poly_divides(Q, A, len1, res, hlen))
                break;
            }
            else
            {
                if (pbits + unlucky >= bound) /* if we reach the bound with one prime */
                { 
                    _fmpz_vec_content(hc, res, hlen);

                   /* divide by content */
                   _fmpz_vec_scalar_divexact_fmpz(res, res, hlen, hc);
                   break;
                }

                if (pbits >= bits_small) /* if one prime is already big enough to check */
                {
                    /* divide by content */
                    _fmpz_vec_content(hc, res, hlen);

                    /* correct sign of leading term */
                    if (fmpz_sgn(res + hlen - 1) < 0)
                        fmpz_neg(hc, hc);

                    _fmpz_vec_scalar_divexact_fmpz(res, res, hlen, hc);

                    /* are we done? */
                    if (_fmpz_poly_divides(Q, B, len2, res, hlen) &&
                        _fmpz_poly_divides(Q, A, len1, res, hlen))
                        break;

                    /* no, so multiply by content again */
                    _fmpz_vec_scalar_mul_fmpz(res, res, hlen, hc);
                }
            }

            curr_bits = FLINT_ABS(_fmpz_vec_max_bits(res, hlen));
            fmpz_set_ui(modulus, p);
            n = hlen - 1; /* if we reach this we have a new bound on length of result */
            continue;
        }

        _fmpz_poly_CRT_ui(res, res, hlen, modulus, h, hlen, mod.n, mod.ninv, 1);
        fmpz_mul_ui(modulus, modulus, mod.n);

        new_bits = _fmpz_vec_max_bits(res, hlen);
        new_bits = FLINT_ABS(new_bits);

        if (new_bits == curr_bits || fmpz_bits(modulus) >= bits_small)
        {
            if (!g_pm1)
            {
                _fmpz_vec_content(hc, res, hlen);

                /* correct sign of leading term */
                if (fmpz_sgn(res + hlen - 1) < 0)
                    fmpz_neg(hc, hc);

                /* divide by content */
                _fmpz_vec_scalar_divexact_fmpz(res, res, hlen, hc);      
            }

            if (fmpz_bits(modulus) + unlucky >= bound)
                break;

            /* are we done? */
            if (_fmpz_poly_divides(Q, B, len2, res, hlen) &&
                _fmpz_poly_divides(Q, A, len1, res, hlen))
                break;

            if (!g_pm1) 
            {            
                /* no, so multiply by content again */
                _fmpz_vec_scalar_mul_fmpz(res, res, hlen, hc);
            }
        }

        curr_bits = new_bits;
    }

    fmpz_clear(modulus);
    fmpz_clear(g); 
    fmpz_clear(l); 
    fmpz_clear(hc);

    _nmod_vec_clear(a);
    _nmod_vec_clear(b);
    _nmod_vec_clear(h); 

    /* finally multiply by content */
    _fmpz_vec_scalar_mul_fmpz(res, res, hlen, d);

    fmpz_clear(d);
    _fmpz_vec_clear(A, len1);
    _fmpz_vec_clear(B, len2);
    _fmpz_vec_clear(Q, len1);
}

void
fmpz_poly_gcd_modular(fmpz_poly_t res, const fmpz_poly_t poly1,
              const fmpz_poly_t poly2)
{
    if (poly1->length < poly2->length)
    {
        fmpz_poly_gcd_modular(res, poly2, poly1);
    }
    else /* len1 >= len2 >= 0 */
    {
        const len_t len1 = poly1->length;
        const len_t len2 = poly2->length;
        
        if (len1 == 0) /* len1 = len2 = 0 */
        {
            fmpz_poly_zero(res);
        } 
        else if (len2 == 0) /* len1 > len2 = 0 */
        {
            if (fmpz_sgn(poly1->coeffs + (len1 - 1)) > 0)
                fmpz_poly_set(res, poly1);
            else
                fmpz_poly_neg(res, poly1);
        }
        else /* len1 >= len2 >= 1 */
        {
            /* underscore function automatically aliases */
            fmpz_poly_fit_length(res, len2);
                
            _fmpz_poly_gcd_modular(res->coeffs, poly1->coeffs, len1,
                                    poly2->coeffs, len2);
     
    
            _fmpz_poly_set_length(res, len2);
            _fmpz_poly_normalise(res);
        }
    }
}

