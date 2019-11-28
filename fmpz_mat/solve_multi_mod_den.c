/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpq_mat.h"

mp_limb_t fmpz_mat_find_good_prime_and_solve(nmod_mat_t Xmod,
		                 nmod_mat_t Amod, nmod_mat_t Bmod,
                const fmpz_mat_t A, const fmpz_mat_t B, const fmpz_t det_bound)
{
    mp_limb_t p;
    fmpz_t tested;

    p = UWORD(1) << NMOD_MAT_OPTIMAL_MODULUS_BITS;
    fmpz_init(tested);
    fmpz_one(tested);

    while (1)
    {
        p = n_nextprime(p, 0);
        _nmod_mat_set_mod(Xmod, p);
	_nmod_mat_set_mod(Amod, p);
        _nmod_mat_set_mod(Bmod, p);
        fmpz_mat_get_nmod_mat(Amod, A);
	fmpz_mat_get_nmod_mat(Bmod, B);
        if (nmod_mat_solve(Xmod, Amod, Bmod))
            break;
        fmpz_mul_ui(tested, tested, p);
        if (fmpz_cmp(tested, det_bound) > 0)
        {
            p = 0;
            break;
        }
    }

    fmpz_clear(tested);
    return p;
}

void
_fmpz_mat_solve_multi_mod_den(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B,
                     nmod_mat_t Xmod, nmod_mat_t Amod, nmod_mat_t Bmod,
		                   mp_limb_t p, const fmpz_t N, const fmpz_t D)
{
    fmpz_t bound, pprod, t, u, dmul;
    fmpz_mat_t x, d, AX, Bden;
    fmpq_mat_t x_q;
    fmpz_t xknum, xkden;
    slong i, j, k, jstart = 0, kstart = 0, n, cols;
    int stabilised; /* has CRT stabilised */

    n = A->r;
    cols = B->c;

    fmpz_init(bound);
    fmpz_init(pprod);
    fmpz_init(xknum);
    fmpz_init(xkden);
    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(dmul);

    fmpz_mat_init(Bden, B->r, B->c);
    fmpz_mat_init(AX, B->r, B->c);
    fmpz_mat_init(x, n, cols);
    fmpz_mat_init_set(d, B);

    fmpq_mat_init(x_q, n, cols);

    /* Compute bound for the needed modulus. TODO: if one of N and D
       is much smaller than the other, we could use a tighter bound (i.e. 2ND).
       This would require the ability to forward N and D to the
       CRT routine.
     */
    if (fmpz_cmpabs(N, D) < 0)
        fmpz_mul(bound, D, D);
    else
        fmpz_mul(bound, N, N);
    fmpz_mul_ui(bound, bound, UWORD(2));  /* signs */

    fmpz_set_ui(pprod, p);
    fmpz_mat_set_nmod_mat(x, Xmod);

    i = 1; /* working with i primes */
    
    while (fmpz_cmp(pprod, bound) <= 0)
    {
	stabilised = 1;
        fmpz_one(dmul);

        /* check if stabilised */
        for (j = jstart; j < x->r && stabilised; j++)
	{
           for (k = kstart; k < x->c && stabilised; k++)
           {
	      fmpz_mul(t, dmul, fmpz_mat_entry(x, j, k));
	      fmpz_fdiv_qr(u, t, t, pprod);
			      
              /* set stabilised to success of reconstruction */
              if ((stabilised = _fmpq_reconstruct_fmpz(xknum, xkden, t, pprod)))
              {
                 /* save starting point for next time */
		 jstart = j;
		 kstart = k + 1;
		 
		 if (kstart == x->c)
	            kstart = 0, jstart = j + 1;
		 
		 fmpz_mul(xkden, xkden, dmul);
                 fmpz_set(dmul, xkden);
	      }     
           }
        }

        /* full matrix stabilisation check */
	if (stabilised)
        {
           stabilised = fmpq_mat_set_fmpz_mat_mod_fmpz(x_q, x, pprod);
	   if (stabilised)
           {
	      fmpq_mat_get_fmpz_mat_matwise(X, den, x_q);

              fmpz_mat_mul(AX, A, X);
              fmpz_mat_scalar_mul_fmpz(Bden, B, den);

	      if (fmpz_mat_equal(AX, Bden))
	         goto multi_mod_done;
	   }
        }
	i++;

        while (1)
        {
           p = n_nextprime(p, 1);

           _nmod_mat_set_mod(Xmod, p);
           _nmod_mat_set_mod(Amod, p);
           _nmod_mat_set_mod(Bmod, p);
           fmpz_mat_get_nmod_mat(Amod, A);
           fmpz_mat_get_nmod_mat(Bmod, B);
           if (nmod_mat_solve(Xmod, Amod, Bmod))
              break;
        }

        fmpz_mat_CRT_ui(x, x, pprod, Xmod, 0); 

        fmpz_mul_ui(pprod, pprod, p);
    }

    /* TODO can be changed to one step */
    fmpq_mat_set_fmpz_mat_mod_fmpz(x_q, x, pprod);
    fmpq_mat_get_fmpz_mat_matwise(X, den, x_q);

multi_mod_done:

    fmpz_clear(xknum);
    fmpz_clear(xkden);
    fmpz_clear(bound);
    fmpz_clear(pprod);
    fmpz_clear(dmul);
    fmpz_clear(u);
    fmpz_clear(t);

    fmpq_mat_clear(x_q);

    fmpz_mat_clear(AX);
    fmpz_mat_clear(Bden);
    fmpz_mat_clear(x);
    fmpz_mat_clear(d);
}

int
fmpz_mat_solve_multi_mod_den(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    nmod_mat_t Xmod, Amod, Bmod;
    fmpz_t N, D;
    mp_limb_t p;

    if (!fmpz_mat_is_square(A))
    {
        flint_printf("Exception (fmpz_mat_solve_multi_mod_den). Non-square system matrix.\n");
        flint_abort();
    }

    if (fmpz_mat_is_empty(A) || fmpz_mat_is_empty(B))
        return 1;

    fmpz_init(N);
    fmpz_init(D);
    fmpz_mat_solve_bound(N, D, A, B);

    nmod_mat_init(Amod, A->r, A->c, 1);
    nmod_mat_init(Bmod, B->r, B->c, 1);
    nmod_mat_init(Xmod, B->r, B->c, 1);

    p = fmpz_mat_find_good_prime_and_solve(Xmod, Amod, Bmod, A, B, D);
    if (p != 0)
        _fmpz_mat_solve_multi_mod_den(X, den, A, B, Xmod, Amod, Bmod, p, N, D);

    nmod_mat_clear(Xmod);
    nmod_mat_clear(Bmod);
    nmod_mat_clear(Amod);
    fmpz_clear(N);
    fmpz_clear(D);

    return p != 0;
}
