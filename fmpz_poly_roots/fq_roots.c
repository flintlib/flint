#include "flint.h"
#include "fq_vec.h"
#include "fq_poly.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "fmpz_poly_roots.h"

void
fmpz_poly_roots_fq_init2 (fmpz_poly_roots_fq_t roots, slong n, fq_ctx_t fctx)
{
  roots->x0 = _fq_vec_init (n, fctx);
  roots->multiplicity = flint_malloc (sizeof (slong) * n);
  roots->num = n;
  roots->alloc = n;
}

void
fmpz_poly_roots_fq_clear (fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)
{
  _fq_vec_clear (roots->x0, roots->alloc, fctx);
  flint_free (roots->multiplicity);
}


char *
fmpz_poly_roots_fq_get_str (fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)
{
  char * buffer = NULL;
  size_t buffer_size = 0;
  FILE * out = open_memstream(&buffer, &buffer_size);
  
  _fq_vec_fprint (out, roots->x0, roots->num, fctx);
  fclose(out);

  return buffer;
}

void
fmpz_poly_roots_fq_print (fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)
{
  _fq_vec_print (roots->x0, roots->num, fctx);
}

char *
fmpz_poly_roots_fq_get_str_pretty (fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)
{
  char * buffer = NULL;
  size_t buffer_size = 0;
  FILE * out = open_memstream(&buffer, &buffer_size);
  slong j;
  
  for (j = 0; j < roots->num; j++)
    {
      fq_fprint_pretty (out, roots->x0 + j, fctx);
      flint_fprintf (out, " %wd\n", roots->multiplicity[j]);
    }

  fclose(out);
  
  return buffer;
}

int
fmpz_poly_roots_fq_fprint_pretty (FILE* file, fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)
{
  slong j;
  
  for (j = 0; j < roots->num; j++)
    {
      fq_fprint_pretty (file, roots->x0 + j, fctx);
      flint_fprintf (file, " %wd\n", roots->multiplicity[j]);
    }

  return 1;
}

int
fmpz_poly_roots_fq_print_pretty (fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)
{
  fmpz_poly_roots_fq_fprint_pretty(stdout, roots, fctx);
  return 1;
}

void
fmpz_poly_roots_fq (fmpz_poly_roots_fq_t roots, fmpz_poly_t poly, fq_ctx_t fctx)
{
  slong j, k, num;

  fmpz_mod_ctx_t fmctx;
  fmpz_mod_poly_t mpoly;
  fq_poly_factor_t f;
  fq_t lead;
  fq_poly_t fpoly;
  
  fmpz_mod_ctx_init(fmctx, fq_ctx_prime(fctx));
  fmpz_mod_poly_init (mpoly, fmctx);
  fmpz_mod_poly_set_fmpz_poly (mpoly, poly, fmctx);

  fq_init (lead, fctx);
  fq_poly_factor_init (f, fctx);
  fq_poly_init (fpoly, fctx);
  fq_poly_set_fmpz_mod_poly (fpoly, mpoly, fctx);
  fq_poly_factor(f, lead, fpoly, fctx);
  fq_clear(lead, fctx);
  
  fmpz_mod_poly_clear (mpoly, fmctx);
  fmpz_mod_ctx_clear(fmctx);
  fq_poly_clear (fpoly, fctx);
 
  num = 0;
  
  for (j = 0; j < f->num; j++)
    {
      if (fq_poly_degree (f->poly + j, fctx) == 1)
	num++;
    }
  
  fmpz_poly_roots_fq_init2 (roots, num, fctx);
  
  k = 0;
  
  for (j = 0; j < f->num; j++)
    {
      if (fq_poly_degree (f->poly + j, fctx) == 1)
	{
	  fq_poly_get_coeff (roots->x0 + k, f->poly + j, 0, fctx);
	  fq_neg(roots->x0 + k, roots->x0 + k, fctx);
	  roots->multiplicity[k] = f->exp[j];
	  k++;
	}
    }
}
