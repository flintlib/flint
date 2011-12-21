
#include <stdio.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"

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

    fmpz_poly_t pol;
    fmpz_poly_init(pol);

    nmod_poly_factor_t facs;
    nmod_poly_factor_init(facs);

    fmpz_poly_set_str(pol, "5  190713877264 0 354248 0 1");
    long p = 89;
    fmpz_t P, Q, PQ;

    nmod_poly_t mod_pol;

    *P = p;
    *Q = p;
    *PQ = p*p;

    nmod_poly_init(mod_pol, p);

    fmpz_poly_get_nmod_poly(mod_pol, pol);
    nmod_poly_factor(facs, mod_pol);

    nmod_poly_factor_print(facs); 
    long r = facs->num_factors;
    printf("number of factors found = %ld\n", r);

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

    fmpz_poly_t temp, temp2;
    fmpz_poly_init(temp);
    fmpz_poly_init(temp2);

    nmod_poly_set(gp, facs->factors[0]);
    nmod_poly_set(hp, facs->factors[1]);
    nmod_poly_xgcd(dp, ap, bp, gp, hp);

    nmod_poly_print(gp);
    printf("\n");
    nmod_poly_print(hp);
    printf("\n");

    fmpz_poly_set_nmod_poly(g, gp);
    fmpz_poly_set_nmod_poly(h, hp);
    fmpz_poly_set_nmod_poly(a, ap);
    fmpz_poly_set_nmod_poly(b, bp);

    fmpz_poly_hensel_lift(g, h, a, b, pol, g, h, a, b, P, Q, PQ);
    
    fmpz_poly_mul(temp, g, h);
    fmpz_poly_sub(temp, f, temp);
    fmpz_poly_scalar_mod_fmpz(temp, PQ);
    long result = fmpz_poly_is_zero(temp);
    if (result > 0)
        printf("hensel_lift part worked\n");
    else
        printf("hensel_lift part failed\n");

    fmpz_poly_mul(temp, g, a);
    fmpz_poly_mul(temp1, h, b);
    fmpz_poly_add(temp, temp, temp1);
    fmpz_poly_scalar_mod_fmpz(temp, PQ);
    long result = fmpz_poly_is_one(temp);
    if (result > 0)
        printf("hensel_lift inverse part worked\n");
    else
        printf("hensel_lift inverse part failed\n");

/* Now lifing to 89^3*/
/*
    fmpz_set(P, PQ);
    fmpz_mul(PQ, PQ, Q);
    
    fmpz_poly_hensel_lift_without_inverse(g, h, f, g, h, a, b, P, Q, PQ);
    fmpz_poly_mul(temp, g, h);
    fmpz_poly_sub(temp, f, temp);
    fmpz_poly_scalar_mod_fmpz(temp, PQ);
    long result = fmpz_poly_is_zero(temp);
    if (result > 0)
        printf("hensel_lift_without_inverse worked\n");
    else
        printf("hensel_lift_without_inverse failed\n");

    fmpz_poly_hensel_lift_only_inverse(a, b, f, g, h, a, b, P, Q, PQ);

    fmpz_poly_mul(temp, g, a);
    fmpz_poly_mul(temp1, h, b);
    fmpz_poly_add(temp, temp, temp1);
    fmpz_poly_scalar_mod_fmpz(temp, PQ);
    long result = fmpz_poly_is_one(temp);
    if (result > 0)
        printf("hensel_lift_only_inverse worked\n");
    else
        printf("hensel_lift_only_inverse failed\n");
*/

    fmpz_poly_clear(temp);
    fmpz_poly_clear(temp1);

    nmod_poly_factor_clear(facs);
    fmpz_poly_clear(pol);
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

    return EXIT_SUCCESS;
}

