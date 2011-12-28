#include <stdio.h>
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
    ulong M_bits;
    M_bits = fmpz_bits(lc) + FLINT_ABS(fmpz_poly_max_bits(f)) + len + FLINT_CLOG2(len);
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

            fmpz_poly_get_nmod_poly(F_tmp, f);
            if (F_tmp->length < f->length)
		    {
		        printf("shrank for p=%ld\n", p);
                p = n_nextprime( p, 0);
                nmod_poly_clear(F_d);
                nmod_poly_clear(F_sbo);
                nmod_poly_clear(F_tmp);
                continue;
            }

            nmod_poly_derivative(F_d, F_tmp);
            nmod_poly_gcd(F_sbo, F_tmp, F_d);


            if (nmod_poly_is_one(F_sbo))
                tryme = 0;
            else
		    {
                p = n_nextprime(p, 0);
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

        r = temp_fac->num;

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

        p = n_nextprime(p, 0);
        nmod_poly_clear(F_tmp);
        nmod_poly_factor_clear(temp_fac);
    }

    p = min_p;
    r = fac->num;
    lead_coeff = min_lead_coeff;
   
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
   a = (long) M_bits / FLINT_CLOG2(p);
   ulong zass_a = a;
   
    fmpz_poly_factor_t lifted_fac;
    fmpz_poly_factor_init(lifted_fac);
    long *link;
    fmpz_poly_t *v, *w;
        
    link = malloc((2*r-2) * sizeof(long));
    v = malloc((2*r-2) * sizeof(fmpz_poly_t));
    w = malloc((2*r-2) * sizeof(fmpz_poly_t));
    
    if (r > cutoff){
        printf("not ready for this\n");
        abort();
    }

    for(i = 0; i < 2*r - 2; i++)
    {
        fmpz_poly_init(v[i]);
        fmpz_poly_init(w[i]);
    }

    printf("going with p = %ld to the a = %ld, r = %ld\n", p, zass_a, r);
    ulong prev_exp;
    prev_exp = _fmpz_poly_hensel_start_lift(lifted_fac, link, v, w, f, fac, zass_a);

    nmod_poly_factor_clear(fac);
    nmod_poly_clear(F);

    fmpz_t P;
    fmpz_init(P);
    fmpz_set_ui(P, p);
    fmpz_pow_ui(P, P, zass_a);

    fmpz_poly_factor_zassenhaus_recombination(final_fac, lifted_fac, f, P, exp);
    
    for (i = 0; i < 2*r - 2; i++)
    {
        fmpz_poly_clear(v[i]);
        fmpz_poly_clear(w[i]);
    }
    free(link);
    free(v);
    free(w);
    fmpz_clear(lc);
    fmpz_clear(P);
    fmpz_poly_factor_clear(lifted_fac);   
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

/*        FILE *polyfile;
        polyfile = fopen("examples/fmpz_poly_hensel_P1", "r");

        if (!polyfile)
        {
            printf("Error.  Could not read P1 from file.\n");
            abort();
        }

        fmpz_poly_fread(polyfile, f);*/

        fmpz_poly_set_str(f, "63  1 1 1 -4 -7 -2 -6 -3 -7 18 7 25 -11 95 36 21 16 69 56 35 36 32 33 26 -26 -15 -14 -53 -96 67 72 -67 40 -79 -116 -452 -312 -260 -29 -1393 327 69 -28 -241 230 -54 -309 -125 -74 -450 -69 -3 66 -27 73 68 50 -63 -1290 372 31 -16 2");

        fmpz_poly_factor_zassenhaus(facs, f);
        
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
