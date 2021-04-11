#include "qadic.h"
#include "fmpz_poly_roots.h"

void
fmpz_poly_roots_qadic_init2 (fmpz_poly_roots_qadic_t roots, slong n)
{
  slong j;
  
  roots->x0 = flint_malloc (sizeof (qadic_t) * n);
  
  for (j = 0; j < n; j++)
    {
      qadic_init (roots->x0 + j);
    }
  roots->multiplicity = flint_malloc (sizeof (slong) * n);
  roots->num = n;
  roots->alloc = n;
}

void
fmpz_poly_roots_qadic_clear (fmpz_poly_roots_qadic_t roots)
{
  slong j;
  
  for (j = 0; j < roots->alloc; j++)
    {
      qadic_clear (roots->x0 + j);
    }
  
  flint_free (roots->x0);
  flint_free (roots->multiplicity);
}

char *
fmpz_poly_roots_qadic_get_str_pretty (fmpz_poly_roots_qadic_t roots, qadic_ctx_t qctx)
{
  char * buffer = NULL;
  size_t buffer_size = 0;
  FILE *out = open_memstream(&buffer, &buffer_size);
  slong j;
  
  for (j = 0; j < roots->num; j++)
    {
      qadic_fprint_pretty (out, roots->x0 + j, qctx);
      flint_fprintf (out, " %wd\n", roots->multiplicity[j]);
    }

  fclose(out);

  return buffer;
}


int
fmpz_poly_roots_qadic_fprint_pretty (FILE *file, fmpz_poly_roots_qadic_t roots, qadic_ctx_t qctx)
{
  slong j;
  
  for (j = 0; j < roots->num; j++)
    {
      qadic_fprint_pretty (file, roots->x0 + j, qctx);
      flint_fprintf (file, " %wd\n", roots->multiplicity[j]);
    }

  return 1;
}

int
fmpz_poly_roots_qadic_print_pretty (fmpz_poly_roots_qadic_t roots, qadic_ctx_t qctx)
{
  fmpz_poly_roots_qadic_fprint_pretty (stdout, roots, qctx);
  
  return 1;
}

static void
qadic_hensel_iteration (fmpz_poly_t poly, qadic_t x, qadic_ctx_t ctx)
{
  slong j;
  
  qadic_t tmp, y0, y1;
  
  qadic_init (tmp);
  qadic_init (y0);
  qadic_init (y1);
  
  do
    {
      /* Horner evaluation of poly and poly' at x */
      padic_poly_set_fmpz (y0, poly->coeffs + poly->length - 1, &(ctx->pctx));
      qadic_zero (y1);
      for (j = poly->length - 2; j >= 0; j--)
	{
	  qadic_mul (y1, y1, x, ctx);
	  qadic_add (y1, y1, y0, ctx);
	  qadic_mul (y0, y0, x, ctx);
	  padic_poly_set_fmpz (tmp, poly->coeffs + j, &(ctx->pctx));
	  qadic_add (y0, y0, tmp, ctx);
	}
      /* Newton step: x -> x - poly / poly' */
      qadic_inv (y1, y1, ctx);
      qadic_mul (y1, y1, y0, ctx);
      qadic_sub (x, x, y1, ctx);
    }
  
  while (qadic_val (y0));
  
  qadic_clear (tmp);
  qadic_clear (y0);
  qadic_clear (y1);
}

void
fmpz_poly_roots_qadic (fmpz_poly_roots_qadic_t roots, fmpz_poly_t poly,
		       qadic_ctx_t qctx)
{
  slong j, k;
  
  fq_ctx_t fctx;
  fmpz_poly_roots_fq_t froots;
  qadic_t q, a;
  
  fq_ctx_init (fctx, (qctx->pctx).p, qadic_ctx_degree (qctx), "a");
  fmpz_poly_roots_fq (froots, poly, fctx);
  fmpz_poly_roots_qadic_init2 (roots, froots->num);
 
  qadic_init (q);
  qadic_init (a);
  qadic_gen (q, qctx);
  
  for (j = 0; j < roots->num; j++)
    {
      *(roots->multiplicity + j) = *(froots->multiplicity + j);
      padic_poly_set_fmpz (roots->x0 + j,
			   (froots->x0 + j)->coeffs + (froots->x0 +
						       j)->length - 1,
			   &(qctx->pctx));
      for (k = (froots->x0 + j)->length - 2; k >= 0; k--)
	{
	  qadic_mul (roots->x0 + j, roots->x0 + j, q, qctx);
	  padic_poly_set_fmpz (a, (froots->x0 + j)->coeffs + k,
			       &(qctx->pctx));
	  qadic_add (roots->x0 + j, roots->x0 + j, a, qctx);
	}
      qadic_hensel_iteration (poly, roots->x0 + j, qctx);
    }
  
  qadic_clear (a);
  qadic_clear (q);
  
  fq_ctx_clear (fctx);
  fmpz_poly_roots_fq_clear (froots, fctx);
}
