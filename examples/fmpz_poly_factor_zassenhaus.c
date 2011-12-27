#include <stdio.h>
#include "math.h"
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void fmpz_poly_factor_sq_fr_prim_internal(fmpz_poly_factor_t final_fac, 
								          ulong exp, fmpz_poly_t f, ulong cutoff)
{
    if (f->length <= 1)
        return;
   
    if (f->length == 2)
    {
        fmpz_poly_factor_insert(final_fac, f, exp);
        return;
    }

    ulong len = f->length;
    fmpz_t lc;
    fmpz_init(lc);
    fmpz_set(lc, f->coeffs + len - 1);
    ulong M_bits = 0;
    M_bits = M_bits + fmpz_bits(lc) + FLINT_ABS(fmpz_poly_max_bits(f)) + len + (long) ceil(log2((double) len));
    nmod_poly_t F, F_d, F_sbo, F_tmp;

    int tryme = 1;
    ulong p = 2UL;
    long i, num_primes;
    nmod_poly_factor_t fac, temp_fac;

    ulong min_p = p; 
    unsigned long lead_coeff, min_lead_coeff;
    long r, min_r;
    min_r = len;
    i = 0;

    nmod_poly_factor_init(fac);
    nmod_poly_init(F, p);

    for (num_primes = 1; num_primes < 3; num_primes++)
    {
        tryme  = 1;
        for ( ; i < 200 && tryme == 1; i++)
	    {
            nmod_poly_init(F_tmp, p);
            nmod_poly_init(F_d, p);
            nmod_poly_init(F_sbo, p);
            nmod_poly_clear(F);
            nmod_poly_init(F, p);

            fmpz_poly_to_nmod_poly(F_tmp, f);
            if (F_tmp->length < f->length)
		    {
                p = z_nextprime( p, 0);
                nmod_poly_clear(F_d);
                nmod_poly_clear(F_sbo);
                nmod_poly_clear(F_tmp);
                continue;
            }

            nmod_poly_derivative(F_d, F_tmp);
            nmod_poly_gcd(F_sbo, F_tmp, F_d);


            if (zmod_poly_is_one(F_sbo))
                tryme = 0;
            else
		    {
                p = z_nextprime(p, 0);
                nmod_poly_clear(F_tmp);
            } 
            nmod_poly_clear(F_d);
            nmod_poly_clear(F_sbo);
        }
   
	    if (i == 200)
	    {
            printf("FLINT Warning: wasn't square_free after 200 primes, maybe an error\n");

            nmod_poly_clear(F);
            nmod_poly_factor_clear(fac);
		    fmpz_clear(lc);
         
		    return;
        }

        nmod_poly_factor_init(temp_fac);
        lead_coeff = nmod_poly_factor(temp_fac, F_tmp);

        r = temp_fac->num_factors;

        if (r <= min_r)
        {
            min_r = r;
            min_p = p;
            nmod_poly_factor_clear(fac);
            nmod_poly_factor_init(fac);
            nmod_poly_factor_concat(fac, temp_fac);
            nmod_poly_set(F, F_tmp);
            min_lead_coeff = lead_coeff;
        }

        p = z_nextprime(p, 0);
        nmod_poly_clear(F_tmp);
        nmod_poly_factor_clear(temp_fac);
    }

    p = min_p;
    r = fac->num_factors;
    lead_coeff = min_lead_coeff;

    int use_Hoeij_Novocin = 0;
    int solved_yet = 0;
   
    fmpz_mat_t M;
    ulong bit_r = FLINT_MAX(r, 20);
    long U_exp = (long) ceil(log2((double) bit_r));;
   
    if (r*3 > f->length)
        U_exp = (long) ceil(log2((double) bit_r));

    if (r > cutoff)
    {
        use_Hoeij_Novocin = 1;
        fmpz_mat_init_identity(M, r);

        fmpz_mat_mul_2exp(M, M, U_exp);
        ulong i;
    }

   if (r == 0)
   {
      printf("FLINT Exception: r == 0\n");
      nmod_poly_clear(F);
      nmod_poly_factor_clear(fac);
      fmpz_clear(lc);
      abort();
   }

   if (r == 1)
   {
      fmpz_poly_factor_insert(final_fac, f, exp);
      nmod_poly_clear(F);
      nmod_poly_factor_clear(fac);
      fmpz_clear(lc);
	  return;
   }

   /* Begin Hensel Lifting phase, make the tree in v, w, and link
      Later we'll do a check if use_Hoeij_Novocin (try for smaller a)*/
   ulong a;
   a = (long) ceil((double) M_bits / log2((double) p));
   a = (long) pow((double) 2, ceil(log2( (double) a)));

   ulong zass_a = a;
}

void fmpz_poly_factor_squarefree(fmpz_poly_factor_t fac, fmpz_poly_t F)
{
/*   fmpz_poly_factor_clear(fac);
   fmpz_poly_factor_init(fac); add to contract that this is freshly initialized*/

   fmpz_poly_content(&fac->c, F);

   fmpz_poly_t f;
   fmpz_poly_init(f);

   fmpz_poly_scalar_divexact_fmpz(f, F, &fac->c);

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

void fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t final_fac, fmpz_poly_t G)
{
   if (G->length == 0)
   {
      fmpz_set_ui(&final_fac->c, 0UL);
      
	  return;
   }
   
   if (G->length == 1)
   {
      fmpz_set(&final_fac->c, G->coeffs);
      
	  return;
   }
   
   fmpz_poly_t g;
   fmpz_poly_init(g);

   if (G->length == 2)
   {
      fmpz_poly_content(&final_fac->c, G);
      fmpz_poly_scalar_divexact_fmpz(g, G, &final_fac->c);
      fmpz_poly_factor_insert(final_fac, g, 1UL);
      fmpz_poly_clear(g);
      
	  return;
   }

  /* Does a presearch for a factor of form x^(x_pow) */
   ulong x_pow = 0;
   while (fmpz_is_zero(G->coeffs + x_pow))
      x_pow++; 
 
   if (x_pow != 0)
   {
      fmpz_poly_t temp_x;
      fmpz_poly_init(temp_x);
      fmpz_poly_set_coeff_ui(temp_x, 1, 1);
      fmpz_poly_factor_insert(final_fac, temp_x, x_pow);
      fmpz_poly_clear(temp_x);
   }
   
   fmpz_poly_shift_right(g, G, x_pow);
   fmpz_poly_factor_t sq_fr_fac;

   /* Could make other tests for x-1 or simple things 
    maybe take advantage of the composition algorithm*/
   fmpz_poly_factor_init(sq_fr_fac);
   fmpz_poly_factor_squarefree(sq_fr_fac, g);

   /* Now we can go through and factor each square free one and add it to final factors.*/
   long j;
   for (j = 0; j < sq_fr_fac->num; j++)
      fmpz_poly_factor_sq_fr_prim_internal(final_fac, sq_fr_fac->exp[j], sq_fr_fac->p + j, 10);

   fmpz_poly_factor_clear(sq_fr_fac);
   fmpz_poly_clear(g);
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
