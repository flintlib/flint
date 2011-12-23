#include <stdio.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

static void 
fmpz_poly_tree_hensel_recursive_check(long *link, 
    fmpz_poly_t *v, fmpz_poly_t *w, fmpz_t P, fmpz_poly_t f, 
    long j, const fmpz_t Q, const fmpz_t PQ)
{
    int result;
    fmpz_poly_t temp, temp1;

    if (j < 0) return;

    fmpz_poly_init(temp);
    fmpz_poly_init(temp1);

    /*f, v[j], v[j + 1], w[j], w[j + 1] == f, g, h, a, b*/

    fmpz_poly_mul(temp, v[j], v[j+1]);
    fmpz_poly_sub(temp, f, temp);
    fmpz_poly_scalar_mod_fmpz(temp, temp, PQ);
    
    result = fmpz_poly_is_zero(temp);
    if (result <= 0)
        printf("hensel_lift part failed******\n");

    /*
    fmpz_poly_mul(temp, v[j], w[j]);
    fmpz_poly_mul(temp1, v[j+1], w[j+1]);
    fmpz_poly_add(temp, temp, temp1);
    fmpz_poly_scalar_mod_fmpz(temp, temp, PQ);

    result = fmpz_poly_is_one(temp);
    if (result <= 0)
        printf("hensel_lift inverse part failed*****\n");
     */

    fmpz_poly_clear(temp);
    fmpz_poly_clear(temp1);

    fmpz_poly_tree_hensel_recursive_check(link, v, w, P, v[j],   link[j], Q, PQ);
    fmpz_poly_tree_hensel_recursive_check(link, v, w, P, v[j + 1], link[j + 1], Q, PQ);
}

int main(void)
{
    /* Part I ****************************************************************/
    {
        int result;
        fmpz_poly_t f;
        nmod_poly_factor_t facs;
        nmod_poly_t mod_pol;
        nmod_poly_t gp, hp, ap, bp, dp;
        fmpz_poly_t g, h, a, b;
        fmpz_poly_t temp, temp1;

        const mp_limb_t p = 89;
        fmpz_t P, Q, PQ;

        fmpz_poly_init(f);
        nmod_poly_factor_init(facs);

        fmpz_poly_set_str(f, "5  190713877264 0 354248 0 1");

        *P  = p;
        *Q  = p;
        *PQ = p*p;

        nmod_poly_init(mod_pol, p);

        fmpz_poly_get_nmod_poly(mod_pol, f);
        nmod_poly_factor(facs, mod_pol);

     /* nmod_poly_factor_print(facs);
        long r = facs->num;
        printf("number of factors found = %ld\n", r); */

        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_init(a);
        fmpz_poly_init(b);

        nmod_poly_init(gp, p);
        nmod_poly_init(hp, p);
        nmod_poly_init(ap, p);
        nmod_poly_init(bp, p);
        nmod_poly_init(dp, p);

        fmpz_poly_init(temp);
        fmpz_poly_init(temp1);

        nmod_poly_set(gp, facs->p[0]);
        nmod_poly_set(hp, facs->p[1]);
        nmod_poly_xgcd(dp, ap, bp, gp, hp);

     /* nmod_poly_print(gp);
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

        result = fmpz_poly_is_zero(temp);
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

        /* Now lifing to 89^3 */

        fmpz_set(P, PQ);
        fmpz_mul(PQ, PQ, Q);

        fmpz_poly_hensel_lift_without_inverse(g, h, f, g, h, a, b, P, Q);
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
        nmod_poly_factor_clear(facs);
    }

    /* Part II ***************************************************************/
    {
        FILE *polyfile;
        fmpz_poly_t f;
        nmod_poly_factor_t facs;
        const mp_limb_t p = 89;
        fmpz_t P, Q, PQ;
        long i, r;
        nmod_poly_t mod_pol;
        long *link;
        fmpz_poly_t *v, *w;
        fmpz_poly_factor_t lifted_fac;
        ulong texp;

        fmpz_poly_init(f);
        nmod_poly_factor_init(facs);

        polyfile = fopen("examples/fmpz_poly_hensel_P1", "r");

        if (!polyfile)
        {
            printf("Error.  Could not read P1 from file.\n");
            abort();
        }

        fmpz_poly_fread(polyfile, f);

        *P  = p;
        *Q  = p;
        *PQ = p*p;

        nmod_poly_init(mod_pol, p);

        fmpz_poly_get_nmod_poly(mod_pol, f);

        nmod_poly_factor(facs, mod_pol);

        r = facs->num;

        fmpz_poly_factor_init(lifted_fac);

        link = malloc((2*r-2) * sizeof(long));
        v = malloc((2*r-2) * sizeof(fmpz_poly_t));
        w = malloc((2*r-2) * sizeof(fmpz_poly_t));

        for(i = 0; i < 2*r - 2; i++)
        {
            fmpz_poly_init(v[i]);
            fmpz_poly_init(w[i]);
        }

        texp = _fmpz_poly_start_hensel_lift(lifted_fac, link, v, w, f, facs, 2);

     /* fmpz_poly_build_hensel_tree(link, v, w, facs);
        fmpz_poly_tree_hensel_lift_recursive(link, v, w, f, 2*r - 4, 1, P, Q, PQ); */

        printf("start recursive check\n");
        fmpz_poly_tree_hensel_recursive_check(link, v, w, P, f, 2*r - 4, Q, PQ);
        printf("finish recursive check\n");

        for (i = 0; i < 2*r - 2; i++)
        {
            fmpz_poly_clear(v[i]);
            fmpz_poly_clear(w[i]);
        }
        free(link);
        free(v);
        free(w);
        fmpz_poly_factor_clear(lifted_fac);
        nmod_poly_factor_clear(facs);
        nmod_poly_clear(mod_pol);
        fmpz_poly_clear(f);
    }

    /* Part III **************************************************************/
    if (0)
    {
        long i, target_exp;

        /* Code from flint-1 */
        for (target_exp = 1; target_exp < 17; target_exp++)
        {
            long num_steps = FLINT_MIN(target_exp, 5 + FLINT_BIT_COUNT(target_exp) - 1);
            long *exponents;
            long *copy_exponents;
            long pow = target_exp;
            long max_steps;

            exponents = malloc(num_steps * sizeof(long));
            copy_exponents = malloc(num_steps * sizeof(long));

            for (i = 0; i < num_steps && pow > 1; i++)
            { 
                copy_exponents[i] = pow;
                pow = (pow + 1)/2;
            }

            max_steps = i;
            exponents[0] = 1;
            for (i = 1; i <= max_steps; i++)
                exponents[i] = copy_exponents[max_steps - i];

            num_steps = max_steps + 1;

            for (i = 0; i < num_steps - 2; i++)
                printf("-(%ld, %ld), ", exponents[i], exponents[i+1]);
            if (num_steps > 1)
                printf("+(%ld, %ld), ", exponents[i], exponents[i+1]);
            printf("\n");

            printf("Return value %ld\n", exponents[i]);

            free(exponents);
            free(copy_exponents);
        }

        /* Code from flint-2 */
        for (target_exp = 1; target_exp < 17; target_exp++)
        {
            long n, *e;

            n = FLINT_BIT_COUNT(target_exp - 1) + 1;
            e = malloc(n * sizeof(long));
            for (e[i = 0] = target_exp; e[i] > 1; i++)
                e[i + 1] = (e[i] + 1) / 2;

            for (i--; i >= 1; i--)
            {
                printf("-(%ld, %ld), ", e[i+1], e[i]);
            }
            if (target_exp > 1)
                printf("+(%ld, %ld), ", e[i+1], e[i]);
            printf("\n");

            printf("Return value %ld\n", e[i+1]);

            free(e);
        }
    }

    /* Part IV ***************************************************************/
    if (0)
    {
        long i, prev, current, target;

        /* Code from flint-1 
        for(prev = 1; prev < 10; prev++)
            for (current = prev + 10; current < prev + 30; current++)
                for (target = current + 40; target < current + 60; target++)
        {
            long num_steps = 2 + FLINT_FLOG2(target - prev);

            long *exponents = malloc(num_steps * sizeof(long));
            long pow;

            pow = current;
            exponents[0] = prev;

            for (i = 1; (i < num_steps) && (pow < target); i++)
            { 
                exponents[i] = pow;
                pow = pow * 2;
            }

            num_steps = i;
            if (exponents[num_steps - 1] != target)
            {
                exponents[num_steps] = target;
                num_steps++;
            }
   
            printf("{%ld %ld %ld}\n", prev, current, target);
            printf("0[%ld %ld]\n", exponents[0], exponents[1]);

            for (i = 1; i < num_steps - 2; i++)
                printf("1[%ld %ld]\n", exponents[i], exponents[i+1]);

            printf("2[%ld %ld]\n", exponents[i], exponents[i+1]);

            printf("{%ld}\n", exponents[i]);
            printf("--------\n");
            free(exponents);
        } */

        /* Code from flint-2 */
        for(prev = 1; prev < 10; prev++)
            for (current = prev + 10; current < prev + 30; current++)
                for (target = current + 40; target < current + 60; target++)
        {

            long n = 2 + FLINT_FLOG2(target - prev);

            long *e = malloc(n * sizeof(long));

            for (e[i = 0] = target; e[i] > current; i++)
                e[i + 1] = (e[i] + 1) / 2;
            e[i]   = current;
            e[i+1] = prev;

            printf("{%ld %ld %ld}\n", prev, current, target); fflush(stdout);
            printf("0[%ld %ld]\n", e[i+1], e[i]);

            for (i--; i > 0; i--)
                printf("1[%ld %ld]\n", e[i+1], e[i]);

            printf("2[%ld %ld]\n", e[i+1], e[i]);

            printf("{%ld}\n", e[i+1]);
            printf("--------\n");
            free(e);
        }
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
