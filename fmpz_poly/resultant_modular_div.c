/*
    Copyright (C) 2014 William Hart
    Copyright (C) 2015 Claus Fieker

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "mpn_extras.h"


void _fmpz_poly_resultant_modular_div(fmpz_t res, 
        const fmpz * poly1, slong len1, 
        const fmpz * poly2, slong len2, const fmpz_t divisor, slong nbits)
{
    flint_bitcnt_t pbits;
    slong i, num_primes;
    fmpz_comb_t comb;
    fmpz_comb_temp_t comb_temp;
    fmpz_t ac, bc, l, modulus, div, la, lb;
    fmpz * A, * B, * lead_A, * lead_B;
    mp_ptr a, b, rarr, parr;
    mp_limb_t p, d;
    nmod_t mod;

    if (fmpz_is_zero(divisor))
    {
        fmpz_zero(res);
        return;
    }

    /* special case, one of the polys is a constant */
    if (len2 == 1) /* if len1 == 1 then so does len2 */
    {
        fmpz_pow_ui(res, poly2, len1 - 1);
        fmpz_divexact(res, res, divisor);

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


    fmpz_init(l);
    /* we have, originally 
    // res(p1, p2) = r d 
    // with nbits(r) <= nbits
    // Now
    // res(p1, p2) = ac^(len2-1) bc^(len1-1) res(A, B)
    // So we need to split ac^(len2-1) bc^(len1-1) = xy
    // such that d mod x == 0 and gcd(ac^... /x, y) = 1
    // Then we need to compute res(A,B)/(d/x) 
    // res(p1, p2) =  x y res(A,B) = r d
    // and
    // r = x y res(A,B)/d = y res(A, B)/(d/x)
    // The length of res(A, B) shrinks by length(y)
    */

    if (!fmpz_is_one(ac))
    {
        fmpz_pow_ui(l, ac, len2 - 1);
        fmpz_init(div);
        fmpz_init(la);
        fmpz_gcd(div, l, divisor); /* div = gcd(ac^(len2-1), divisor) */
        fmpz_divexact(la, l, div); /* la = ac^(len2 -1)/gcd           */
        fmpz_divexact(div, divisor, div); /*div /= gcd                */
        nbits = nbits - fmpz_bits(la) + 1;
    } else {
        fmpz_init_set(div, divisor);
    }
    
    if (!fmpz_is_one(bc))
    {
        fmpz_init(lb);
        fmpz_pow_ui(lb, bc, len1 - 1);
        fmpz_gcd(l, lb, div);
        fmpz_divexact(lb, lb, l);
        fmpz_divexact(div, div, l);
        nbits = nbits - fmpz_bits(lb) + 1;
    }

    
    /* get product of leading coefficients */
    lead_A = A + len1 - 1;
    lead_B = B + len2 - 1;
    fmpz_mul(l, lead_A, lead_B);

    fmpz_init(modulus);
    fmpz_set_ui(modulus, 1);
    fmpz_zero(res);

    /* make space for polynomials mod p */
    a = _nmod_vec_init(len1);
    b = _nmod_vec_init(len2);

    pbits = FLINT_BITS - 1;
    p = (UWORD(1)<<pbits);

    nbits = FLINT_MAX((slong)0, nbits);

    num_primes = (nbits+pbits-1)/pbits;
    if (num_primes <= 0) num_primes = 1;

    parr = _nmod_vec_init(num_primes);
    rarr = _nmod_vec_init(num_primes);

    
    for(i=0; i< num_primes; )
    {
        /* get new prime and initialise modulus */
        p = n_nextprime(p, 0);

        if (fmpz_fdiv_ui(l, p) == 0) 
            continue;
        d = fmpz_fdiv_ui(div, p);
        if (d==0)
            continue;
        d = n_invmod(d, p);

        nmod_init(&mod, p);

        /* reduce polynomials modulo p */
        _fmpz_vec_get_nmod_vec(a, A, len1, mod);
        _fmpz_vec_get_nmod_vec(b, B, len2, mod);

        /* compute resultant over Z/pZ */
        rarr[i] = _nmod_poly_resultant(a, len1, b, len2, mod);
        rarr[i] = n_mulmod2_preinv(rarr[i], d, mod.n, mod.ninv);
        parr[i++] = p;
    }

    fmpz_comb_init(comb, parr, num_primes);
    fmpz_comb_temp_init(comb_temp, comb);
    
    fmpz_multi_CRT_ui(res, rarr, comb, comb_temp, 1);
        
    fmpz_clear(modulus);
    fmpz_comb_temp_clear(comb_temp);
    fmpz_comb_clear(comb);
        
    _nmod_vec_clear(a);
    _nmod_vec_clear(b);

    _nmod_vec_clear(parr);
    _nmod_vec_clear(rarr);
    
    /* finally multiply by powers of content */
    if (!fmpz_is_one(ac))
    {
        fmpz_mul(res, res, la);
        fmpz_clear(la);
    }
    
    if (!fmpz_is_one(bc))
    {
        fmpz_mul(res, res, lb);
        fmpz_clear(lb);
    }

    fmpz_clear(l); 
    fmpz_clear(div); 
    
    _fmpz_vec_clear(A, len1);
    _fmpz_vec_clear(B, len2);

    fmpz_clear(ac);
    fmpz_clear(bc);
}

void
fmpz_poly_resultant_modular_div(fmpz_t res, const fmpz_poly_t poly1,
              const fmpz_poly_t poly2, const fmpz_t divisor, slong nbits)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    
    if (len1 == 0 || len2 == 0 || fmpz_is_zero(divisor))
        fmpz_zero(res);
    else if (len1 >= len2)
        _fmpz_poly_resultant_modular_div(res, poly1->coeffs, len1,
                                          poly2->coeffs, len2, divisor, nbits);
    else
    {
        _fmpz_poly_resultant_modular_div(res, poly2->coeffs, len2, 
                                          poly1->coeffs, len1, divisor, nbits);  
        if ((len1 > 1) && (!(len1 & WORD(1)) & !(len2 & WORD(1))))
            fmpz_neg(res, res);
    }
}

