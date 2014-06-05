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


void _fmpz_poly_resultant_modular(fmpz_t res, const fmpz * poly1, slong len1, 
                                        const fmpz * poly2, slong len2)
{
    mp_bitcnt_t bits1, bits2, bound, pbits, curr_bits = 0;   
    fmpz_t ac, bc, l, modulus;
    fmpz * A, * B, * lead_A, * lead_B;
    mp_ptr a, b;
    mp_limb_t p, r;
    nmod_t mod;
    
    /* special case, one of the polys is zero */
    if (len2 == 0) /* if len1 == 1 then so does len2 */
    {
       fmpz_zero(res);

       return;
    }
    
    /* special case, one of the polys is a constant */
    if (len2 == 1) /* if len1 == 1 then so does len2 */
    {
        fmpz_pow_ui(res, poly2, len1 - 1);

        return;
    }
    
    fmpz_init(ac);
    fmpz_init(bc);
    
    /* compute content of poly1 and poly2 */
    _fmpz_vec_content(ac, poly1, len1);
    _fmpz_vec_content(bc, poly2, len2);
    
    /* divide poly1 and poly2 by their content */
    A = _fmpz_vec_init(len1);
    B = _fmpz_vec_init(len2);
    _fmpz_vec_scalar_divexact_fmpz(A, poly1, len1, ac);
    _fmpz_vec_scalar_divexact_fmpz(B, poly2, len2, bc);
    
    /* get product of leading coefficients */
    fmpz_init(l);
    
    lead_A = A + len1 - 1;
    lead_B = B + len2 - 1;
    fmpz_mul(l, lead_A, lead_B);

    /* set size of first prime */
    pbits = FLINT_BITS - 1;
    p = (UWORD(1)<<pbits);

    /* get bound on size of resultant */
    bits1 = FLINT_ABS(_fmpz_vec_max_bits(A, len1)); 
    bits2 = FLINT_ABS(_fmpz_vec_max_bits(B, len2));

    /* bound on number of bits of 2(deg1 + deg2)! */
    bound = (len1 + len2 - 1)*FLINT_BIT_COUNT((10*(len1 + len2 - 1) + 26)/27) + 3;

    /* Upper bound Hadamard bound */
    bound += (len1 - 1)*bits2 + (len2 - 1)*bits1;
    
    fmpz_init(modulus);
    fmpz_set_ui(modulus, 1);
    fmpz_zero(res);

    /* make space for polynomials mod p */
    a = _nmod_vec_init(len1);
    b = _nmod_vec_init(len2);
    
    while (curr_bits < bound)
    {
        /* get new prime and initialise modulus */
        p = n_nextprime(p, 0);
        if (fmpz_fdiv_ui(l, p) == 0)
            continue;
        
        curr_bits += pbits;

        nmod_init(&mod, p);

        /* reduce polynomials modulo p */
        _fmpz_vec_get_nmod_vec(a, A, len1, mod);
        _fmpz_vec_get_nmod_vec(b, B, len2, mod);

        /* compute resultant over Z/pZ */
        r = _nmod_poly_resultant(a, len1, b, len2, mod);

        fmpz_CRT_ui(res, res, modulus, r, mod.n, 1);
        fmpz_mul_ui(modulus, modulus, mod.n);
    }

    fmpz_clear(modulus);
    
    _nmod_vec_clear(a);
    _nmod_vec_clear(b);
    
    /* finally multiply by powers of content */
    fmpz_pow_ui(l, ac, len2 - 1);
    fmpz_mul(res, res, l);
    fmpz_pow_ui(l, bc, len1 - 1);
    fmpz_mul(res, res, l);
    
    fmpz_clear(l); 
    
    _fmpz_vec_clear(A, len1);
    _fmpz_vec_clear(B, len2);

    fmpz_clear(ac);
    fmpz_clear(bc);
}

void
fmpz_poly_resultant_modular(fmpz_t res, const fmpz_poly_t poly1,
              const fmpz_poly_t poly2)
{
    if (poly1->length < poly2->length)
    {
        fmpz_poly_resultant_modular(res, poly2, poly1);
        if ((poly1->length & 1) == 0 && (poly2->length & 1) == 0)
           fmpz_neg(res, res);
    }
    else /* len1 >= len2 >= 0 */
    {
        _fmpz_poly_resultant_modular(res, poly1->coeffs, poly1->length, 
           poly2->coeffs, poly2->length);
    }
}

