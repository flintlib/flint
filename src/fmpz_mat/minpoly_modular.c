/*
    Copyright (C) 2013 Sebastian Pancratz
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>

#include "ulong_extras.h"
#include "nmod_mat.h"
#include "nmod_poly.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"

#define MINPOLY_M_LOG2E  1.44269504088896340736  /* log2(e) */

static inline long double _log2(const long double x)
{
    return log(x) * MINPOLY_M_LOG2E;
}

slong _fmpz_mat_minpoly_small(fmpz * rop, const fmpz_mat_t op)
{
    slong len = 0;

    if (op->r == 0)
    {
        fmpz_one(rop + 0);
        len = 1;
    }
    else if (op->r == 1)
    {
        fmpz_one(rop + 1);
        fmpz_neg(rop + 0, &(op->rows[0][0]));
        len = 2;
    }

    return len;
}

void _fmpz_mat_bound_ovals_of_cassini(fmpz_t b, const fmpz_mat_t op)
{
   slong n = op->r, i, j;
   fmpz * v1;
   fmpz_t t, q, r1, r2;

   fmpz_init(t);
   fmpz_init(q);
   fmpz_init(r1);
   fmpz_init(r2);

   v1 = _fmpz_vec_init(n);

   /* |A| [1,1,...,1]^T */
   for (i = 0; i < n; i++)
   {
      for (j = 0; j < n; j++)
      {
         if (fmpz_sgn(fmpz_mat_entry(op, i, j)) >= 0)
            fmpz_add(v1 + i, v1 + i, fmpz_mat_entry(op, i, j));
         else
            fmpz_sub(v1 + i, v1 + i, fmpz_mat_entry(op, i, j));
      }
   }

   for (i = 0; i < n; i++)
   {
      fmpz_zero(t);

      /* q_i */
      fmpz_abs(t, fmpz_mat_entry(op, i, i));

      if (fmpz_cmp(t, q) > 0)
         fmpz_set(q, t);

      /* r_i */
      fmpz_sub(t, v1 + i, t);

      if (fmpz_cmp(t, r2) > 0)
      {
         fmpz_swap(t, r2);

         if (fmpz_cmp(r2, r1) > 0)
            fmpz_swap(r2, r1);
      }
   }

   fmpz_mul(r1, r1, r2);

   fmpz_sqrtrem(b, r2, r1);

   if (!fmpz_is_zero(r2))
      fmpz_add_ui(b, b, 1);

   fmpz_add(b, b, q);

   _fmpz_vec_clear(v1, n);
   fmpz_clear(r1);
   fmpz_clear(r2);
   fmpz_clear(t);
   fmpz_clear(q);
}

slong _fmpz_mat_minpoly_modular(fmpz * rop, const fmpz_mat_t op)
{
    const slong n = op->r;
    slong len = 0, oldlen = 0;

    if (n < 2)
    {
        return _fmpz_mat_minpoly_small(rop, op);
    }
    else
    {
        /*
            If $A$ is an $n \times n$ matrix with spectral radius
            bound by b, the coefficients of the minimal polynomial have
            at most $\ceil{d\log_2(b)}$ bits if $d \leq b$. Otherwise
            if $d > 0$ it has at most
            $\min{\ceil{d/2\log_2(bd)}, \ceil{d\log_2(2b)}}$
            bits, where $d$ is the degree of the minimal polynomial.
            See Lemma 3.1 in Dumas, "Bounds on the coefficients of the
            characteristic and minimal polynomials", 2007.
        */
        slong bound;
        double b1, b2, b3, bb;

        slong pbits  = FLINT_BITS - 1, i, j;
        mp_limb_t p = (UWORD(1) << pbits);
        ulong * P, * Q;

        fmpz_mat_t v1, v2, v3;
        fmpz * rold;
        fmpz_t m;

        if (fmpz_mat_is_zero(op))
        {
           fmpz_set_ui(rop + 0, 1);
           return 1;
        }

        /* Determine the bound in bits */
        {
            fmpz_t b;

            fmpz_init(b);

            _fmpz_mat_bound_ovals_of_cassini(b, op);
            bb = fmpz_get_d(b);
            if (bb <= 1.0)
               bb = 1.0;

            b1 = _log2(bb);
            b2 = _log2(bb*n)/2;
            b3 = _log2(bb*2);

            if (b3 < b2)
               b2 = b3;

            bound = n <= bb ? (slong) ceil(n*b1) : (slong) ceil(n*b2);
            bound = bound + 1;

            fmpz_clear(b);
        }

        P = (ulong *) flint_calloc(n, sizeof(ulong));
        Q = (ulong *) flint_calloc(n, sizeof(ulong));
        rold = (fmpz *) _fmpz_vec_init(n + 1);
        fmpz_mat_init(v1, n, 1);
        fmpz_mat_init(v2, n, 1);
        fmpz_mat_init(v3, n, 1);

        fmpz_init_set_ui(m, 1);

        oldlen = 0;
        len = 0;

        for ( ; fmpz_bits(m) <= bound; )
        {
            nmod_mat_t mat;
            nmod_poly_t poly;

            p = n_nextprime(p, 0);

            nmod_mat_init(mat, n, n, p);
            nmod_poly_init(poly, p);

            for (i = 0; i < n; i++)
               P[i] = 0;

            fmpz_mat_get_nmod_mat(mat, op);
            nmod_mat_minpoly_with_gens(poly, mat, P);

            len = poly->length;

            if (oldlen != 0 && len > oldlen)
            {
               /* all previous primes were bad, discard */

               fmpz_one(m);
               oldlen = len;

               for (i = 0; i < n + 1; i++)
                  fmpz_zero(rop + i);

               for (i = 0; i < n; i++)
                  Q[i] = 0;
            } else if (len < oldlen)
            {
               /* this prime was bad, skip */

               nmod_mat_clear(mat);
               nmod_poly_clear(poly);

               continue;
            }

            for (i = 0; i < n; i++)
               Q[i] |= P[i];

            _fmpz_poly_CRT_ui(rop, rop, n + 1, m, poly->coeffs,
                              poly->length, poly->mod.n, poly->mod.ninv, 1);

            fmpz_mul_ui(m, m, p);

            /* check if stabilised */
            for (i = 0; i < len; i++)
            {
               if (!fmpz_equal(rop + i, rold + i))
                  break;
            }

            for (j = 0; j < len; j++)
               fmpz_set(rold + j, rop + j);

            if (i == len) /* stabilised */
            {
               for (i = 0; i < n; i++)
               {
                  if (Q[i] == 1)
                  {
                     fmpz_mat_zero(v1);
                     fmpz_mat_zero(v3);

                     fmpz_set_ui(fmpz_mat_entry(v1, i, 0), 1);

                     for (j = 0; j < len; j++)
                     {
                        fmpz_mat_scalar_mul_fmpz(v2, v1, rop + j);
                        fmpz_mat_add(v3, v3, v2);

                        if (j != len - 1)
                        {
                           fmpz_mat_mul(v2, op, v1);
                           fmpz_mat_swap(v1, v2);
                        }
                     }

                     /* check f(A)v = 0 */
                     for (j = 0; j < n; j++)
                     {
                        if (!fmpz_is_zero(v3->rows[j] + 0))
                            break;
                     }

                     if (j != n)
                        break;
                  }
               }

               /* if f(A)v = 0 for all generators v, we are done */
               if (i == n)
               {
                  nmod_mat_clear(mat);
                  nmod_poly_clear(poly);
                  break;
               }
            }

            nmod_mat_clear(mat);
            nmod_poly_clear(poly);
        }

        flint_free(P);
        flint_free(Q);
        fmpz_mat_clear(v2);
        fmpz_mat_clear(v1);
        fmpz_mat_clear(v3);
        fmpz_clear(m);
        _fmpz_vec_clear(rold, n);
    }

    return len;
}

void fmpz_mat_minpoly_modular(fmpz_poly_t cp, const fmpz_mat_t mat)
{
    slong len;

    fmpz_poly_fit_length(cp, mat->r + 1);

    len = _fmpz_mat_minpoly_modular(cp->coeffs, mat);

    _fmpz_poly_set_length(cp, len);
}
