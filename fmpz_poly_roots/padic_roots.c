#include <stdio.h>

#include "padic.h"
#include "padic_poly.h"
#include "fmpz_poly_roots.h"


void
fmpz_poly_roots_padic_init2 (fmpz_poly_roots_padic_t roots, slong n)
{
  slong j;
  
  roots->x0 = flint_malloc (sizeof (padic_t) * n);
  
  for (j = 0; j < n; j++)
    {
      padic_init (roots->x0 + j);
    }
  
  roots->multiplicity = flint_malloc (sizeof (slong) * n);
  roots->num = n;
  roots->alloc = n;
}

void
fmpz_poly_roots_padic_clear (fmpz_poly_roots_padic_t roots)
{
  slong j;
  
  for (j = 0; j < roots->alloc; j++)
    {
      padic_clear (roots->x0 + j);
    }
  
  flint_free (roots->x0);
  flint_free (roots->multiplicity);
}

char*
fmpz_poly_roots_padic_get_str (fmpz_poly_roots_padic_t roots, padic_ctx_t pctx)
{
  char * buffer = NULL;
  size_t buffer_size = 0;
  FILE *out = open_memstream(&buffer, &buffer_size);
  
  slong j;
  
  for (j = 0; j < roots->num; j++)
    {
      padic_fprint(out, roots->x0 + j, pctx);
      flint_fprintf (out, " %wd\n", roots->multiplicity[j]);
    }
  
  fclose(out);
  
  return buffer;
}

int
fmpz_poly_roots_padic_fprint (FILE *file, fmpz_poly_roots_padic_t roots, padic_ctx_t pctx)
{
  slong j;
  
  for (j = 0; j < roots->num; j++)
    {
      padic_fprint(file, roots->x0 + j, pctx);
      flint_fprintf (file, " %wd\n", roots->multiplicity[j]);
    }

  return 1;
}

int
fmpz_poly_roots_padic_print (fmpz_poly_roots_padic_t roots, padic_ctx_t pctx)
{
  fmpz_poly_roots_padic_fprint(stdout, roots, pctx);
  return 1;
}

static void
padic_hensel_iteration (fmpz_poly_t poly, padic_t x, padic_ctx_t ctx)
{
  slong j;
  padic_t tmp, y0, y1;
  
  padic_init (tmp);
  padic_init (y0);
  padic_init (y1);
  
  do
    {
      /* Horner evaluation of poly and poly' at x */
      padic_set_fmpz (y0, poly->coeffs + poly->length - 1, ctx);
      padic_zero (y1);
      for (j = poly->length - 2; j >= 0; j--)
	{
	  padic_mul (y1, y1, x, ctx);
	  padic_add (y1, y1, y0, ctx);
	  padic_mul (y0, y0, x, ctx);
	  padic_set_fmpz (tmp, poly->coeffs + j, ctx);
	  padic_add (y0, y0, tmp, ctx);
	}
      /* Newton step: x -> x - poly / poly' */
      padic_inv (y1, y1, ctx);
      padic_mul (y1, y1, y0, ctx);
      padic_sub (x, x, y1, ctx);
    }
  while (padic_val (y0));
  
  padic_clear (tmp);
  padic_clear (y0);
  padic_clear (y1);
}

void
fmpz_poly_roots_padic (fmpz_poly_roots_padic_t roots, fmpz_poly_t poly,
		       padic_ctx_t pctx)
{
  slong j;
  fq_ctx_t fctx;
  fmpz_poly_roots_fq_t froots;
  
  fq_ctx_init (fctx, pctx->p, 1, "a");
  
  fmpz_poly_roots_fq (froots, poly, fctx);
  fmpz_poly_roots_padic_init2 (roots, froots->num);
  
  for (j = 0; j < roots->num; j++)
    {
      *(roots->multiplicity + j) = *(froots->multiplicity + j);
      padic_set_fmpz (roots->x0 + j, (froots->x0 + j)->coeffs, pctx);
      padic_hensel_iteration (poly, roots->x0 + j, pctx);
    }
  
  fq_ctx_clear (fctx);
  fmpz_poly_roots_fq_clear (froots, fctx);
}
