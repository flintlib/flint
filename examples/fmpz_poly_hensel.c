
#include <stdio.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"

void fmpz_poly_rec_tree_hensel_check(long * link, fmpz_poly_t * v, fmpz_poly_t * w, 
	             fmpz_t P, fmpz_poly_t f, long j, fmpz_t Q, 
				       fmpz_t PQ)
{
    if (j < 0) return;

    fmpz_poly_t temp, temp1;
    fmpz_poly_init(temp);
    fmpz_poly_init(temp1);

    //f, v[j], v[j + 1], w[j], w[j + 1] == f, g, h, a, b

    fmpz_poly_mul(temp, v[j], v[j+1]);
    fmpz_poly_sub(temp, f, temp);
    fmpz_poly_scalar_mod_fmpz(temp, temp, PQ);
    long result = fmpz_poly_is_zero(temp);
    if (result <= 0)
        printf("hensel_lift part failed******\n");

    fmpz_poly_mul(temp, v[j], w[j]);
    fmpz_poly_mul(temp1, v[j+1], w[j+1]);
    fmpz_poly_add(temp, temp, temp1);
    fmpz_poly_scalar_mod_fmpz(temp, temp, PQ);
    result = fmpz_poly_is_one(temp);
    if (result <= 0)
        printf("hensel_lift inverse part failed*****\n");

    fmpz_poly_clear(temp);
    fmpz_poly_clear(temp1);

       
    fmpz_poly_rec_tree_hensel_check(link, v, w, P, v[j],   link[j], Q, PQ);
    fmpz_poly_rec_tree_hensel_check(link, v, w, P, v[j + 1], link[j + 1], Q, PQ);
}

int main(void)
{

    /*
    3 89  15 64 1 ^ 1
    3 89  15 25 1 ^ 1
    number of factors found = 2
    3 89  15 64 1
    3 89  15 25 1
    3  190713877039 0 354843
    p = 89
    3  6868 7273 1
    3  6868 648 1

    */
    {
    fmpz_poly_t f;
    fmpz_poly_init(f);

    nmod_poly_factor_t facs;
    nmod_poly_factor_init(facs);

    fmpz_poly_set_str(f, "5  190713877264 0 354248 0 1");
    long p = 89;
    fmpz_t P, Q, PQ;

    nmod_poly_t mod_pol;

    *P = p;
    *Q = p;
    *PQ = p*p;

    nmod_poly_init(mod_pol, p);

    fmpz_poly_get_nmod_poly(mod_pol, f);
    nmod_poly_factor(facs, mod_pol);

/*    nmod_poly_factor_print(facs); */
    long r = facs->num_factors;
/*    printf("number of factors found = %ld\n", r); */

    nmod_poly_t gp, hp, ap, bp, dp;
    fmpz_poly_t g, h, a, b;

    fmpz_poly_init(g);
    fmpz_poly_init(h);
    fmpz_poly_init(a);
    fmpz_poly_init(b);

    nmod_poly_init(gp, p);
    nmod_poly_init(hp, p);
    nmod_poly_init(ap, p);
    nmod_poly_init(bp, p);
    nmod_poly_init(dp, p);

    fmpz_poly_t temp, temp1;
    fmpz_poly_init(temp);
    fmpz_poly_init(temp1);

    nmod_poly_set(gp, facs->factors[0]);
    nmod_poly_set(hp, facs->factors[1]);
    nmod_poly_xgcd(dp, ap, bp, gp, hp);

/*    nmod_poly_print(gp);
      printf("\n");
      nmod_poly_print(hp);
      printf("\n"); */

    fmpz_poly_set_nmod_poly(g, gp);
    fmpz_poly_set_nmod_poly(h, hp);
    fmpz_poly_set_nmod_poly(a, ap);
    fmpz_poly_set_nmod_poly(b, bp);

    fmpz_poly_hensel_lift(g, h, a, b, f, g, h, a, b, P, Q, PQ);
    
    fmpz_poly_mul(temp, g, h);
    fmpz_poly_sub(temp, f, temp);
    fmpz_poly_scalar_mod_fmpz(temp, temp, PQ);
    long result = fmpz_poly_is_zero(temp);
    if (result > 0)
        printf("hensel_lift part worked\n");
    else
        printf("hensel_lift part failed\n");

    fmpz_poly_mul(temp, g, a);
    fmpz_poly_mul(temp1, h, b);
    fmpz_poly_add(temp, temp, temp1);
    fmpz_poly_scalar_mod_fmpz(temp, temp, PQ);
    result = fmpz_poly_is_one(temp);
    if (result > 0)
        printf("hensel_lift inverse part worked\n");
    else
        printf("hensel_lift inverse part failed\n");

/* Now lifing to 89^3*/

    fmpz_set(P, PQ);
    fmpz_mul(PQ, PQ, Q);
    
    fmpz_poly_hensel_lift_without_inverse(g, h, f, g, h, a, b, P, Q, PQ);
    fmpz_poly_mul(temp, g, h);
    fmpz_poly_sub(temp, f, temp);
    fmpz_poly_scalar_mod_fmpz(temp, temp, PQ);
    result = fmpz_poly_is_zero(temp);
    if (result > 0)
        printf("hensel_lift_without_inverse worked\n");
    else
        printf("hensel_lift_without_inverse failed\n");

    fmpz_poly_hensel_lift_only_inverse(a, b, f, g, h, a, b, P, Q, PQ);

    fmpz_poly_mul(temp, g, a);
    fmpz_poly_mul(temp1, h, b);
    fmpz_poly_add(temp, temp, temp1);
    fmpz_poly_scalar_mod_fmpz(temp, temp, PQ);
    result = fmpz_poly_is_one(temp);
    if (result > 0)
        printf("hensel_lift_only_inverse worked\n");
    else
        printf("hensel_lift_only_inverse failed\n");

    fmpz_poly_clear(temp);
    fmpz_poly_clear(temp1);

    nmod_poly_factor_clear(facs);
    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_poly_clear(h);
    fmpz_poly_clear(a);
    fmpz_poly_clear(b);
    nmod_poly_clear(mod_pol);
    nmod_poly_clear(gp);
    nmod_poly_clear(hp);
    nmod_poly_clear(ap);
    nmod_poly_clear(bp);
    nmod_poly_clear(dp);
    }
    
    {
    fmpz_poly_t f;
    fmpz_poly_init(f);

    nmod_poly_factor_t facs;
    nmod_poly_factor_init(facs);

    FILE * polyfile = fopen("P1_flint", "r");
    fmpz_poly_fread(polyfile, f);
    long p = 89;
    fmpz_t P, Q, PQ;

    nmod_poly_t mod_pol;

    *P = p;
    *Q = p;
    *PQ = p*p;

    nmod_poly_init(mod_pol, p);

    fmpz_poly_get_nmod_poly(mod_pol, f);
    nmod_poly_factor(facs, mod_pol);

/*    nmod_poly_factor_print(facs); */
    long r = facs->num_factors;
/*    printf("number of factors found = %ld\n", r); */

/* fmpz_poly_build_hensel_tree(long * link, fmpz_poly_t *v, fmpz_poly_t *w, 
                                 const nmod_poly_factor_t fac);*/
    long link[2*r - 2];
    fmpz_poly_t v[2*r - 2];
    fmpz_poly_t w[2*r - 2];

    for(long i = 0; i < 2*r - 2; i++)
    {
        fmpz_poly_init(v[i]);
        fmpz_poly_init(w[i]);
    }

    fmpz_poly_rec_tree_hensel_lift(link, v, w, P, f, 2*r - 4, 1, Q, PQ);
   
    for (long i = 0; i < 2*r - 2; i++)
    {
        fmpz_poly_clear(v[i]);
        fmpz_poly_clear(w[i]);
    }

    printf("start rec check\n");
    fmpz_poly_rec_tree_hensel_check(link, v, w, P, f, 2*r - 4, Q, PQ);
    printf("finish rec check\n");

    nmod_poly_factor_clear(facs);
    nmod_poly_clear(mod_pol);
    fmpz_poly_clear(f);
        
    }
    
    return EXIT_SUCCESS;
}

