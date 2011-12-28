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

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "mpn_extras.h"


#define __mul(C, lenC, A, lenA, B, lenB)                        \
do {                                                            \
    if ((lenA) != 0 && (lenB) != 0)                             \
    {                                                           \
        if ((lenA) >= (lenB))                                   \
            _nmod_poly_mul((C), (A), (lenA), (B), (lenB), mod); \
        else                                                    \
            _nmod_poly_mul((C), (B), (lenB), (A), (lenA), mod); \
        (lenC) = (lenA) + (lenB) - 1;                           \
    }                                                           \
    else                                                        \
    {                                                           \
        (lenC) = 0;                                             \
    }                                                           \
} while (0)

/*
    Set res = s*f + t*g where res is the gcd of f and g
Aliasing of res with f and g is permitted
 */
void nmod_poly_xgcd_hgcd(nmod_poly_t res, nmod_poly_t s, nmod_poly_t t, 
                         const nmod_poly_t f, const nmod_poly_t g)
{
    const nmod_t mod = f->mod;
    const mp_limb_t p = mod.n;

    int sign;
    mp_limb_t a;

    if (f->length == 0)
    {
        nmod_poly_fit_length(t, 1);
        if (g->length == 0)
        {
            nmod_poly_zero(res);
            t->coeffs[0] = 1L;
        } 
        else 
        {
            mp_limb_t Z = n_invmod(g->coeffs[g->length - 1], p);
            nmod_poly_scalar_mul_nmod(res, g, Z);
            t->coeffs[0] = Z;
        }
        t->length = 1;
        nmod_poly_zero(s);
        return;
    }

    if (g->length == 0)
    {
        nmod_poly_fit_length(s, 1);
        mp_limb_t Z = n_invmod(f->coeffs[f->length - 1], p);
        nmod_poly_scalar_mul_nmod(res, f, Z);
        s->coeffs[0] = Z;
        s->length = 1;
        nmod_poly_zero(t);
        return;
    }

    if (f->length == 1)
    {
        a = n_invmod(f->coeffs[0], p);
        nmod_poly_set_coeff_ui(s, 0, a);
        s->length = 1;
        nmod_poly_set_coeff_ui(res, 0, 1L);
        res->length = 1;
        nmod_poly_zero(t);
        return;
    }

    if (g->length == 1)
    {
        a = n_invmod(g->coeffs[0], p);
        nmod_poly_set_coeff_ui(t, 0, a);
        t->length = 1;
        nmod_poly_set_coeff_ui(res, 0, 1L);
        res->length = 1;
        nmod_poly_zero(s);
        return;
    }
   	
	long CUTOFF;
	long bits = FLINT_BIT_COUNT(p);
	if (bits <= 8) CUTOFF = NMOD_POLY_SMALL_GCD_CUTOFF;
	else CUTOFF = NMOD_POLY_GCD_CUTOFF;
	
	nmod_poly_t h, j, q, r, u0, u1, temp, temp2;
	nmod_poly_init(q, p);
   nmod_poly_init(r, p);

	/* g = 0*f + 1*g, r = 1*f - q * g, s is 0, t = 1 */
	nmod_poly_divrem(q, r, f, g); 

   nmod_poly_zero(s);
    nmod_poly_one(t);

	if (r->length == 0)
	{
		/* t is already initialised, s is already 0 */
		mp_limb_t Z = n_invmod(g->coeffs[g->length - 1], p);
		t->coeffs[0] = Z;
		nmod_poly_scalar_mul_nmod(res, g, Z);

	   nmod_poly_clear(q);
	   nmod_poly_clear(r);
		return;
	}

    mp_ptr R[4];
    long lenR[4];

    R[0] = _nmod_vec_init(g->length);
    R[1] = _nmod_vec_init(g->length);
    R[2] = _nmod_vec_init(g->length);
    R[3] = _nmod_vec_init(g->length);

    nmod_poly_init(j, p);
	nmod_poly_init(h, p);
	nmod_poly_init(u0, p);
	nmod_poly_init(u1, p);
	nmod_poly_init(temp, p);
	nmod_poly_init(temp2, p);

    nmod_poly_fit_length(j, g->length);
    nmod_poly_fit_length(h, g->length);
	
	/*
	  Let R = (a b)
	          (c d) then
			 
			 (h j) = (g r)(d -c)
                       (-b a)  if sign > 0 

		     j = -c*g + a*r = -c*(s*f + ?*g) + a*(t*f + ?*g) 
			    = (a*t - c*s)*f + ?*g = a*t*f at this point as s = 0
			 i.e. send t -> a
			  h = d*g - b*r = d*(s*f + ?*g) - b*(t*f + ?*g)
			    = (d*s - b*t)*f + ?*g = -b*t*f at this point as s = 0
			 i.e. send s-> -b
	*/
printf("%ld %ld\n", g->length, r->length);fflush(stdout);
    sign = _nmod_poly_hgcd(R, lenR, h->coeffs, &(h->length), j->coeffs, &(j->length), g->coeffs, g->length, r->coeffs, r->length, mod);

    nmod_poly_fit_length(s, lenR[1]);
    _nmod_vec_neg(s->coeffs, R[1], lenR[1], mod);
    s->length = lenR[1];

    nmod_poly_fit_length(t, lenR[0]);
    _nmod_vec_set(t->coeffs, R[0], lenR[0]);
    t->length = lenR[0];

    if (sign < 0L) 
    {
        nmod_poly_neg(s, s);
        nmod_poly_neg(t, t);
    }

	while (j->length != 0)
	{
      /* r = h - q * j = s*f - q*t*f + ?*g
		   j = t*f + ?*g
			i.e. s->t, t ->s - q*t
		*/
		nmod_poly_divrem(q, r, h, j);
	   nmod_poly_mul(temp, q, t);
		nmod_poly_swap(s, t);
		nmod_poly_sub(t, t, temp);

		if (r->length == 0)
	   {
		   
			/*
			   now res = s*f + ?*g
				so compute ? = (res - s*f)/g
			*/
			nmod_poly_mul(temp, s, f);
			nmod_poly_sub(t, j, temp);
			nmod_poly_div(t, t, g);
		   
			mp_limb_t Z = n_invmod(j->coeffs[j->length - 1], p);
			nmod_poly_scalar_mul_nmod(s, s, Z);
			nmod_poly_scalar_mul_nmod(t, t, Z);
			nmod_poly_scalar_mul_nmod(res, j, Z);

            _nmod_vec_clear(R[0]);	      
            _nmod_vec_clear(R[1]);	      
            _nmod_vec_clear(R[2]);	      
            _nmod_vec_clear(R[3]);	      

	      nmod_poly_clear(u0);
	      nmod_poly_clear(u1);
         nmod_poly_clear(j);
	      nmod_poly_clear(h);
    	   nmod_poly_clear(q);
	      nmod_poly_clear(r);
	      nmod_poly_clear(temp);
	      nmod_poly_clear(temp2);
       	return;
	   }

		if (j->length < CUTOFF)
	   {
		   if ((res == f) || (res == g))
				nmod_poly_xgcd_euclidean(temp2, u0, u1, j, r);
			else
				nmod_poly_xgcd_euclidean(res, u0, u1, j, r);
			/* 
			   we have res = u0*j + u1*r
			   and j = s*f + ?*g
				    r = t*f + ?*g
				i.e. s -> u0*s + u1*t
			*/
			nmod_poly_mul(s, s, u0);
			nmod_poly_mul(temp, t, u1);
			nmod_poly_add(s, s, temp);
			
			/*
			   now res = s*f + ?*g
				so compute ? = (res - s*f)/g
			*/
			nmod_poly_mul(temp, s, f);
			if ((res == f) || (res == g))
			   nmod_poly_sub(t, temp2, temp);
			else 
            nmod_poly_sub(t, res, temp);
			nmod_poly_div(t, t, g);

			/* 
			   Note res and hence also s and t are normalised correctly
			   i.e. res is monic already
			*/
			
			if ((res == f) || (res == g)) 
				nmod_poly_set(res, temp2);
			

            _nmod_vec_clear(R[0]);
            _nmod_vec_clear(R[1]);
            _nmod_vec_clear(R[2]);
            _nmod_vec_clear(R[3]);

      	nmod_poly_clear(u0);
      	nmod_poly_clear(u1);
	      nmod_poly_clear(j);
	      nmod_poly_clear(h);
   	   nmod_poly_clear(q);
	      nmod_poly_clear(r);
      	nmod_poly_clear(temp);
	      nmod_poly_clear(temp2);
			
      	return;
	   }

    sign = _nmod_poly_hgcd(R, lenR, h->coeffs, &(h->length), j->coeffs, &(j->length), j->coeffs, j->length, r->coeffs, r->length, mod);

		/*
		    j' = -c*j + a*r = -c*(s*f + ?*g) + a*(t*f + ?*g) 
			    = (a*t - c*s)*f + ?*g 
			 i.e. send t -> a*t - c*s
			 h' = d*j - b*r = d*(s*f + ?*g) - b*(t*f + ?*g)
			    = (d*s - b*t)*f + ?*g
			 i.e. send s-> d*s - b*t
		*/

        nmod_poly_fit_length(temp, lenR[1] + t->length - 1);
        __mul(temp->coeffs, temp->length, R[1], lenR[1], t->coeffs, t->length);


        nmod_poly_fit_length(temp2, s->length + lenR[3] - 1);
        __mul(temp2->coeffs, temp2->length, s->coeffs, s->length, R[3], lenR[3]);


		if (sign > 0L) nmod_poly_sub(s, temp2, temp);
		else nmod_poly_sub(s, temp, temp2);

        nmod_poly_fit_length(temp, lenR[2] + s->length - 1);
        __mul(temp->coeffs, temp->length, R[2], lenR[2], s->coeffs, s->length);

        nmod_poly_fit_length(temp2, t->length + lenR[0] - 1);
        __mul(temp2->coeffs, temp2->length, t->coeffs, t->length, R[0], lenR[0]);

		if (sign > 0L) nmod_poly_sub(t, temp2, temp);
		else nmod_poly_sub(t, temp, temp2);
	}

	/*
	   now res = s*f + ?*g
		so compute ? = (res - s*f)/g
	*/
	nmod_poly_mul(temp, s, f);
	nmod_poly_sub(t, h, temp);
	nmod_poly_div(t, t, g);

	mp_limb_t Z = n_invmod(h->coeffs[h->length - 1], p);
	nmod_poly_scalar_mul_nmod(s, s, Z);
   nmod_poly_scalar_mul_nmod(t, t, Z);
   nmod_poly_scalar_mul_nmod(res, h, Z);

    _nmod_vec_clear(R[0]);
    _nmod_vec_clear(R[1]);
    _nmod_vec_clear(R[2]);
    _nmod_vec_clear(R[3]);

	nmod_poly_clear(u0);
	nmod_poly_clear(u1);
	nmod_poly_clear(j);
	nmod_poly_clear(h);
	nmod_poly_clear(q);
	nmod_poly_clear(r);
	nmod_poly_clear(temp);
	nmod_poly_clear(temp2);
}

