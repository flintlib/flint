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

void
create_subproduct_tree(fmpz * prod, const fmpz * arr, mp_limb_t len, fmpz_t mod)
{
    int i, j;
    mp_limb_t uppow, uplen, treelen;
    fmpz * tree;

    uppow = FLINT_CLOG2(len);
    uplen = 1 << uppow;
    treelen = (uplen << 1) - 2;

    tree = _fmpz_vec_init(treelen);

    _fmpz_vec_set(tree, arr, len);

    for (i = len; i < uplen; i++)
        fmpz_set_ui(tree + i, 1);

    /* Build subproduct tree */

    for (i = uplen; i < treelen; i++)
    {
        fmpz_mul(tree + i, tree + ((i - uplen) << 1), tree + ((i - uplen) << 1) + 1);
        fmpz_mod(tree + i, tree + i, mod);
    }

    for (i = 0; i < treelen; i += 2)
        fmpz_swap(tree + i, tree + i + 1);

    /* Go back, multiplying */

    for (i = treelen - 1; i >= uplen; i --)
    {
        fmpz_mul(tree + ((i - uplen) << 1), tree + ((i - uplen) << 1), tree + i);
        fmpz_mod(tree + ((i - uplen) << 1), tree + ((i - uplen) << 1), mod);

        fmpz_mul(tree + ((i - uplen) << 1) + 1, tree + ((i - uplen) << 1) + 1, tree + i);
        fmpz_mod(tree + ((i - uplen) << 1) + 1, tree + ((i - uplen) << 1) + 1, mod);

    }

    _fmpz_vec_set(prod, tree, len);
    _fmpz_vec_clear(tree, treelen);
}

int
fmpz_factor_ecm_stage_II_one_FFT(fmpz_t f, mp_limb_t B1, mp_limb_t B2, mp_limb_t P,
                           fmpz_t n, ecm_t ecm_inf)
{

    fmpz_t g, tim, Qx, Qz, Rx, Rz, Qdx, Qdz, a, b;
    mp_limb_t mmin, mmax, maxj, mdiff;
    int i, ret;
    mp_limb_t num_roots;

    fmpz * arrx, * arrz, * prod, * zs;
    fmpz * roots2, * evals, * poly;
    fmpz_poly_struct ** tree;

    mmin = (B1 + (P/2)) / P;
    mmax = ((B2 - P/2) + P - 1)/P;      /* ceil */
    maxj = (P + 1)/2; 
    mdiff = mmax - mmin + 1;
    num_roots = (maxj + 1)/2;
    ret = 0;

    arrx = _fmpz_vec_init(num_roots);
    arrz = _fmpz_vec_init(num_roots);
    poly = _fmpz_vec_init(num_roots + 1);
    prod = _fmpz_vec_init(num_roots + mdiff);
    zs = _fmpz_vec_init(num_roots + mdiff);
    evals = _fmpz_vec_init(mdiff);
    roots2 = _fmpz_vec_init(mdiff);
    tree = _fmpz_mod_poly_tree_alloc(mdiff);

    fmpz_init(Qx);
    fmpz_init(Qz);
    fmpz_init(Qdx);
    fmpz_init(Qdz);
    fmpz_init(Rx);
    fmpz_init(Rz);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(tim);

    fmpz_init_set_ui(g, 1);

    /* arr[0] = Q0 */
    fmpz_set(arrx, ecm_inf->x);  
    fmpz_set(arrz, ecm_inf->z);

    /* [Qx :: Qz] = 2Q0 */
    fmpz_factor_ecm_double(Qx, Qz, arrx, arrz, n, ecm_inf);

    /* arr[1] = 3Q0 */
    fmpz_factor_ecm_add(arrx + 1, arrz + 1, Qx, Qz, arrx, arrz, arrx, arrz,
                          n, ecm_inf);

    /* For each odd i (i > 3) , compute i * Q0 [x0 :: z0] */
    /* Note, stored consecutively in arrx and arrz */

    /* We are adding 2Q0 every time. Need to calculate all i's 
       as (i - 1)Q0 is required for (i + 1)Q0 */

    for (i = 2; i < num_roots; i += 1)
    {
        /* iQ0 = (i - 2)Q0 + 2Q0 
           Differnce is (i - 4)Q0 */

        fmpz_factor_ecm_add(arrx + i, arrz + i, arrx + i - 1, arrz + i - 1,
                              Qx, Qz, arrx + i - 2, arrz + i - 2, n, ecm_inf);

    }

    /* vector zs [0 : num_roots] has z coordinates of small steps */

    for (i = 0; i < num_roots; i++)
        fmpz_set(zs + i, arrz + i);

    /* Q = D * Q_0 */
    fmpz_set_ui(tim, P);
    fmpz_factor_ecm_mul_montgomery_ladder(Qx, Qz, ecm_inf->x, ecm_inf->z,
                                            tim, n, ecm_inf);

    /* R = mmin * Q */
    fmpz_set_ui(tim, mmin);
    fmpz_factor_ecm_mul_montgomery_ladder(Rx, Rz, Qx, Qz, tim, n, ecm_inf);

    /* Qd = (mmin - 1) * Q */
    fmpz_set_ui(tim, mmin - 1);
    fmpz_factor_ecm_mul_montgomery_ladder(Qdx, Qdz, Qx, Qz, tim, n, ecm_inf);

    /* main stage II step */


    for (i = 0; i < mdiff; i ++)
    {
        /* roots2 has the x coordinate of giant steps (the points of evaluations) */
        /* vector zs [num_roots  : mdiff] has z coordinates of giant steps */

        fmpz_set(roots2 + i, Rx);
        fmpz_set(zs + num_roots + i, Rz);

        fmpz_set(a, Rx);
        fmpz_set(b, Rz);

        /* R = R + Q    
           difference is stored in Qd, initially (Mmin - 1)Q */

        fmpz_factor_ecm_add(Rx, Rz, Rx, Rz, Qx, Qz, Qdx, Qdz, n, ecm_inf);

        fmpz_set(Qdx, a);
        fmpz_set(Qdz, b);

    }

    /* zs has all the z coordinates, creating subproduct tree */

    create_subproduct_tree(prod, zs, num_roots + mdiff, n);

    for (i = 0; i < num_roots; i++)
    {
        fmpz_mul(arrx + i, prod + i, arrx + i);
        fmpz_mod(arrx + i, arrx + i, n);
    }

    for (i = 0; i < mdiff; i++)
    {
        fmpz_mul(roots2 + i, prod + num_roots + i, roots2 + i);
        fmpz_mod(roots2 + i, roots2 + i, n);
    }   

    /* create poly from roots */

    _fmpz_mod_poly_product_roots_fmpz_vec(poly, arrx, num_roots, n);

    /* building product tree for *mdiff* number of giant steps */

    _fmpz_mod_poly_tree_build(tree, roots2, mdiff, n);

    /* Evaluating the poly with small steps as coefficients (roots) at all big step points */

    _fmpz_mod_poly_evaluate_fmpz_vec_fast_precomp(evals, poly, num_roots + 1, tree, mdiff, n);

    /* Checking for gcd's */
    
    for (i = 0; i < mdiff; i++)
    {
        fmpz_gcd(f, evals + i, n);

        if (!fmpz_is_one(f) && fmpz_cmp(f, n))
        {
            ret = 1;
            goto cleanup;
        }   
    }

    cleanup:

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

    _fmpz_vec_clear(arrx, num_roots);
    _fmpz_vec_clear(arrz, num_roots);
    _fmpz_vec_clear(poly, num_roots + 1);
    _fmpz_vec_clear(prod, num_roots + mdiff);
    _fmpz_vec_clear(zs, num_roots + mdiff);
    _fmpz_vec_clear(evals, mdiff);
    _fmpz_vec_clear(roots2, mdiff);
    _fmpz_mod_poly_tree_free(tree, mdiff);

    return ret;
}
