/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2020 Claus Fieker

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string.h>

#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "padic.h"
#include "qadic.h"


extern int flint_conway_polynomials [];

void qadic_ctx_init(qadic_ctx_t ctx,
                           const fmpz_t p, slong d, slong min, slong max, 
                           const char *var, enum padic_print_mode mode)
{
    unsigned int position;
    flint_rand_t state;
    fmpz_mod_poly_t poly;
    slong i, j;
    fmpz_mod_ctx_t ctxp;

    if (fmpz_cmp_ui(p, 109987) <= 0)
    {  
      for (position = 0; flint_conway_polynomials[position] != 0;
                               position += 3 + flint_conway_polynomials[position + 1])
      {
          /* Different prime? */
          if (fmpz_cmp_ui(p, flint_conway_polynomials[position]))
              continue;

          /* Same degree? */
          if (d == flint_conway_polynomials[position + 1])
          {
              /* Find number of non-zero coefficients */
              ctx->len = 1;

              for (i = 0; i < d; i++)
              {
                  if (flint_conway_polynomials[position + 2 + i])
                      ctx->len ++;
              }

              ctx->a = _fmpz_vec_init(ctx->len);
              ctx->j = flint_malloc(ctx->len*sizeof(slong));

              /* Copy the polynomial */
              j = 0;

              for (i = 0; i < d; i++)
              {
                  int coeff = flint_conway_polynomials[position + 2 + i];

                  if (coeff)
                  {
                      fmpz_set_ui(ctx->a + j, coeff);
                      ctx->j[j] = i;
                      j++;
                  }
              }

              fmpz_set_ui(ctx->a + j, 1);
              ctx->j[j] = d;

              /* Complete the initialisation of the context */
              padic_ctx_init(&ctx->pctx, p, min, max, mode);

              ctx->var = flint_malloc(strlen(var) + 1);
              strcpy(ctx->var, var);

              return;
          }
      }
    }  

    flint_randinit(state);

    fmpz_mod_ctx_init(ctxp, p);
    fmpz_mod_poly_init2(poly, d + 1, ctxp);
    
    fmpz_mod_poly_randtest_sparse_irreducible(poly, state, d + 1, ctxp);
    
    flint_randclear(state);

    /* Find number of non-zero coefficients */
    ctx->len = 1;

    for (i = 0; i < d; i++)
    {
       if (!fmpz_is_zero(poly->coeffs + i))
          ctx->len ++;
    }

    ctx->a = _fmpz_vec_init(ctx->len);
    ctx->j = flint_malloc(ctx->len*sizeof(slong));

    /* Copy the polynomial */
    j = 0;

    for (i = 0; i < d; i++)
    {
       if (!fmpz_is_zero(poly->coeffs+i))
       {
           fmpz_set(ctx->a + j, poly->coeffs + i);
           ctx->j[j] = i;
           j++;
       }
    }

    fmpz_set_ui(ctx->a + j, 1);
    ctx->j[j] = d;

    /* Complete the initialisation of the context */
    padic_ctx_init(&ctx->pctx, p, min, max, mode);

    ctx->var = flint_malloc(strlen(var) + 1);
    strcpy(ctx->var, var);

    fmpz_mod_poly_clear(poly, ctxp);
    fmpz_mod_ctx_clear(ctxp);
}

