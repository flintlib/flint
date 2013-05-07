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

    Copyright (C) 2012 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <float.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"
#include "profiler.h"

int
main(void)
{
    long len, iter, i;
    fmpz_mod_poly_t f;
    fmpz * roots;
    fmpz_poly_struct ** tree;
    fmpz_t N;
    timeit_t t;
    char N_str[201] = "29799904256775982671863388319999573561548825027149399972531599612392671227006866151136667908641695103422986028076864929902803267437351318167549013218980573566942647077444419419003164546362008247462049";
    
    flint_rand_t state;
    flint_randinit(state);
       
    len = 36865;

    fmpz_init(N);
    fmpz_set_str(N, N_str, 10);

    roots = _fmpz_vec_init(len);
    
    for (i = 0; i < len; i++)
       fmpz_randm(roots + i, state, N);

    fmpz_mod_poly_init(f, N);
    
    /*
        Time tree 
    */
    
    timeit_start(t);
                
    for (iter = 0; iter < 10; iter++)
    {
       tree = _fmpz_mod_poly_tree_alloc(len);
       _fmpz_mod_poly_tree_build(tree, roots, len, N);
       _fmpz_mod_poly_tree_free(tree, len);
    }

    timeit_stop(t);

    printf("len = %ld, time = %ldms\n", len, ((long) t->cpu)/10);
        
    fmpz_mod_poly_clear(f);
    
    _fmpz_vec_clear(roots, len);

    fmpz_clear(N);
    
    flint_randclear(state);
}
