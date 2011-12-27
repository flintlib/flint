#include <stdio.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void fmpz_poly_factor_squarefree(fmpz_poly_factor_t fac, fmpz_poly_t F)
{
   fmpz_poly_content(&fac->c, F);

   fmpz_poly_t f;
   fmpz_poly_init(f);

   fmpz_poly_scalar_divexact_fmpz(f, F, &fac->c);

   fmpz_poly_factor_clear(fac);
   fmpz_poly_factor_init(fac);

   if (f->length == 1)
   {
      fmpz_poly_clear(f);
      return;
   }

   fmpz_poly_t d, v, w, s, t1;
   fmpz_poly_init(d);
   fmpz_poly_init(v);
   fmpz_poly_init(w);
   fmpz_poly_init(s);
   fmpz_poly_init(t1);

   fmpz_poly_derivative(t1, f);
   fmpz_poly_gcd(d, f, t1);

   if (d->length == 1)
   {
      fmpz_poly_factor_insert(fac, f, 1);

      fmpz_poly_clear(d);
      fmpz_poly_clear(v);
      fmpz_poly_clear(w);
      fmpz_poly_clear(s);
      fmpz_poly_clear(t1);
      fmpz_poly_clear(f);
   
	  return;
   }

   long i;

   fmpz_poly_div(v, f, d);
   fmpz_poly_div(w, t1, d);
   

   for (i = 1; ; i++)
   {
      fmpz_poly_derivative(t1, v);
      fmpz_poly_sub(s, w, t1);
      
	  if (s->length == 0)
	  {
         if (v->length > 1)
            fmpz_poly_factor_insert(fac, v, i);
         
		 fmpz_poly_clear(d);
         fmpz_poly_clear(v);
         fmpz_poly_clear(w);
         fmpz_poly_clear(s);
         fmpz_poly_clear(t1);
         fmpz_poly_clear(f);
         return;
      }

      fmpz_poly_gcd(d, v, s);
      fmpz_poly_div(v, v, d);
      fmpz_poly_div(w, s, d);

      if (d->length > 1)
         fmpz_poly_factor_insert(fac, d, i);
   }

   fmpz_poly_clear(f);
   fmpz_poly_clear(d);
   fmpz_poly_clear(v);
   fmpz_poly_clear(w);
   fmpz_poly_clear(s);
   fmpz_poly_clear(t1);
   
   return;
}

int main(void)
{
    /* Part I ****************************************************************/
    {
        fmpz_poly_t f;
        fmpz_poly_factor_t facs;

        fmpz_poly_init(f);
        fmpz_poly_factor_init(facs);

        fmpz_poly_set_str(f, "5  190713877264 0 354248 0 1");

        fmpz_poly_factor_squarefree(facs, f);
        
        fmpz_poly_factor_print(facs);
        printf(" was facs\n");
        
        fmpz_poly_clear(f);
        fmpz_poly_factor_clear(facs);
    }
    return EXIT_SUCCESS;
}

/*
int 
fmpz_poly_hensel_checker(fmpz_poly_factor_t lifted_fac, fmpz_poly_t f, fmpz_t P)
{
    if (lifted_fac->num == 0)
    {
        if (F->length > 0)
            return 0;
        else
            return 1;
    }

    fmpz_mod_poly_t product, temp;
    fmpz_mod_poly_init(product, P);
    fmpz_mod_poly_init(temp, P);

    fmpz_t lead_coeff;
    fmpz_init(lead_coeff);

    fmpz_poly_to_fmpz_mod_poly(product, lifted_fac->p[0]);

    long i;
    for (i = 1; i < lifted_fac->num; i++)
    {
        fmpz_poly_to_fmpz_mod_poly(temp, lifted_fac->p[i]);
        fmpz_mod_poly_mul(product, product, temp);
    }

    fmpz_poly_to_fmpz_mod_poly(temp, F);

    fmpz_set(lead_coeff, F->coeffs + F->length - 1);

    fmpz_mod_poly_scalar_mul(product, product, lead_coeff);

    _fmpz_mod_poly_reduce_coeffs(product);
    _fmpz_mod_poly_reduce_coeffs(temp);
   
    int res = fmpz_mod_poly_equal(product, temp);
   
    fmpz_clear(lead_coeff);
    fmpz_mod_poly_clear(temp);
    fmpz_mod_poly_clear(product);
   
    return res;
}
*/
