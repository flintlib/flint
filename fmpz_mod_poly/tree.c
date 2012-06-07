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

    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 William Hart

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

fmpz_poly_struct ** _fmpz_mod_poly_tree_alloc(long len)
{
    fmpz_poly_struct ** tree = NULL;

    long i, j, height = FLINT_CLOG2(len);

    tree = flint_malloc(sizeof(fmpz_poly_struct *) * (height + 1));
    
    tree[0] = flint_malloc(len*sizeof(fmpz_poly_struct));
    
    for (j = 0; j < len; j++)
        fmpz_poly_init(tree[0] + j);

    len = (1<<(FLINT_FLOG2(len)));
        
    for (i = 1; i <= height; i++, len /= 2)
    {
        tree[i] = flint_malloc(len*sizeof(fmpz_poly_struct));

        for (j = 0; j < len; j++)
           fmpz_poly_init(tree[i] + j);
    }

    return tree;
}

void _fmpz_mod_poly_tree_free(fmpz_poly_struct ** tree, long len)
{
    long i, j, height = FLINT_CLOG2(len);

    for (j = 0; j < len; j++)
        fmpz_poly_clear(tree[0] + j);

    flint_free(tree[0]);

    len = (1<<(FLINT_FLOG2(len)));
    
    for (i = 1; i <= height; i++, len /= 2)
    {
        for (j = 0; j < len; j++)
            fmpz_poly_clear(tree[i] + j);

        flint_free(tree[i]);
    }

    flint_free(tree);
}

void
_fmpz_mod_poly_tree_build(fmpz_poly_struct ** tree, const fmpz * roots, long len, const fmpz_t mod)
{
    long height, i, jump, j, k;
    long len2 = (1<<(FLINT_FLOG2(len)));
    long * group = malloc(len2*sizeof(long));
    
    if (len == 0)
        return;

    height = FLINT_CLOG2(len);
    
    /* zeroth level, (x-a) */
    for (i = 0; i < len; i++)
    {
        fmpz_poly_set_coeff_ui(tree[0] + i,  1, 1);
        if (!fmpz_is_zero(roots + i))
           fmpz_sub((tree[0] + i)->coeffs, mod, roots + i);
    }
    
    group[0] = len;
    jump = len2;
    while (jump > 1)
    {
        for (i = 0; i < len2; i += jump)
        {
           long full = group[i];
           long half = (group[i] + 1)/2;
           group[i] = half;
           group[i + jump/2] = full - half;
        }
        jump /= 2;
    }
    
    j = 0;
    for (i = 0; i < len2; i++)
    {
        if (group[i] == 1)
        {
            fmpz_poly_set(tree[1] + i, tree[0] + j);
            j++;
        } else
        {
            fmpz_poly_mul(tree[1] + i, tree[0] + j, tree[0] + j + 1);
            _fmpz_vec_scalar_mod_fmpz((tree[1] + i)->coeffs, (tree[1] + i)->coeffs, (tree[1] + i)->length - 1, mod);
        }
    }
    
    len2 /= 2;
    for (k = 1; k < height - 1; k++)
    {
        for (i = 0; i < len2; i++)
        {
            fmpz_poly_mul(tree[k + 1] + i, tree[k] + 2*i, tree[k] + 2*i + 1);
            _fmpz_vec_scalar_mod_fmpz((tree[k + 1] + i)->coeffs, (tree[k + 1] + i)->coeffs, (tree[k + 1] + i)->length - 1, mod);

        }
        len2 /= 2;
    }
}
