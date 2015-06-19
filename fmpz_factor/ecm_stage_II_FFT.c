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

    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

/* Outer wrapper for ECM 
   makes calls to stage I and stage II (one) */

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"


int
fmpz_factor_ecm_stage_II_FFT(fmpz_t f, mp_limb_t B1, mp_limb_t B2, mp_limb_t P,
                             fmpz_t n, ecm_t ecm_inf)
{
    fmpz_t g, tim, Qx, Qz, Rx, Rz, Qdx, Qdz, a, b;
    mp_limb_t mmin, mmax, maxj, mdiff;
    int i, j, ret;

    mp_limb_t num_roots;

    fmpz * arrx, * arrz;

    fmpz * roots, * roots2, * evals;
    fmpz_poly_struct ** tree;

    mmin = (B1 + (P/2)) / P;
    mmax = ((B2 - P/2) + P - 1)/P;      /* ceil */
    maxj = (P + 1)/2; 
    mdiff = mmax - mmin + 1;

    arrx = _fmpz_vec_init(maxj + 1);
    arrz = _fmpz_vec_init(maxj + 1);

    fmpz_init(tim);
    fmpz_init(Qx);
    fmpz_init(Qz);
    fmpz_init(Qdx);
    fmpz_init(Qdz);
    fmpz_init(Rx);
    fmpz_init(Rz);
    fmpz_init(a);
    fmpz_init(b);

    fmpz_init_set_ui(g, 1);

    ret = 0;

    /* arr[1] = Q0 */
    fmpz_set(arrx + 1, ecm_inf->x);  
    fmpz_set(arrz + 1, ecm_inf->z);

    /* arr[2] = 2Q0 */
    fmpz_factor_ecm_double(arrx + 2, arrz + 2, arrx + 1, arrz + 1, n, ecm_inf);

    /* arr[3] = 3Q0 */
    fmpz_factor_ecm_add(arrx + 3, arrz + 3, arrx + 3, arrz + 3, arrx + 1, arrz + 1, 
                        arrx + 1, arrz + 1, n, ecm_inf);

    /* For each odd j (j > 3) , compute j * Q0 [x0 :: z0] */

    /* We are adding 2Q0 every time. Need to calculate all j's 
       as (j - 2)Q0 is required for (j + 2)Q0 */

    for (j = 5; j <= maxj; j += 2)
    {
        /* jQ0 = (j - 2)Q0 + 2Q0 
           Differnce is (j - 4)Q0 */

        fmpz_factor_ecm_add(arrx + j, arrz + j, arrx + j - 2, arrz + j - 2, 
                            arrx + 2, arrz + 2, arrx + j - 4, arrz + j - 4, 
                            n, ecm_inf);
    }


    num_roots = (maxj + 1)/2;

	roots = _fmpz_vec_init(num_roots);
	evals = _fmpz_vec_init(num_roots);
	roots2 = _fmpz_vec_init(mdiff);

	/* roots has the small steps (the j's) */

	for (i = 0; i < num_roots; i++)
		fmpz_set(roots + i, arrx + 2*i + 1);
    
    /* Q = D * Q_0 */
    fmpz_set_ui(tim, P);
    fmpz_factor_ecm_mul_montgomery_ladder(Qx, Qz, ecm_inf->x, ecm_inf->z, tim, n, ecm_inf);
    
    /* R = mmin * Q */
    fmpz_set_ui(tim, mmin);
    fmpz_factor_ecm_mul_montgomery_ladder(Rx, Rz, Qx, Qz, tim, n, ecm_inf);

    /* Qd = (mmin - 1) * Q */
    fmpz_set_ui(tim, mmin - 1);
    fmpz_factor_ecm_mul_montgomery_ladder(Qdx, Qdz, Qx, Qz, tim, n, ecm_inf);
                
    /* main stage II step */

    for (i = mmin; i <= mmax; i ++)
    {
    	/* roots2 has the giant steps (the points of evaluations) */

    	fmpz_set(roots2 + i - mmin, Rx);
        fmpz_set(a, Rx);
        fmpz_set(b, Rz);

        /* R = R + Q    
           difference is stored in Qd, initially (Mmin - 1)Q */

        fmpz_factor_ecm_add(Rx, Rz, Rx, Rz, Qx, Qz, Qdx, Qdz, n, ecm_inf);

        fmpz_set(Qdx, a);
        fmpz_set(Qdz, b);

    }

    /* building product tree for *mdiff* number of giant steps */
    tree = _fmpz_mod_poly_tree_alloc(mdiff);
    _fmpz_mod_poly_tree_build(tree, roots2, mdiff, n);

    /* Evaluating the poly with small steps as coefficients (roots) at all big step points */
    _fmpz_mod_poly_evaluate_fmpz_vec_fast_precomp(evals, roots, num_roots, tree, mdiff, n);

    /* Checking for gcd's */
    
	for (i = 0; i < num_roots; i++)
	{
		fmpz_gcd(f, n, evals + i);

		if (!fmpz_is_zero(f) && !fmpz_is_one(f))
		{
			ret = 1;
			break;
		}
	}

	_fmpz_vec_clear(roots, num_roots);
	_fmpz_vec_clear(evals, num_roots);
	_fmpz_vec_clear(roots2, mdiff);
	_fmpz_mod_poly_tree_free(tree, mdiff);

    fmpz_clear(tim);
    fmpz_clear(Qx);
    fmpz_clear(Qz);
    fmpz_clear(Qdx);
    fmpz_clear(Qdz);
    fmpz_clear(Rx);
    fmpz_clear(Rz);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(g);

    _fmpz_vec_clear(arrx, maxj + 1);
    _fmpz_vec_clear(arrz, maxj + 1);

    return ret;
}

