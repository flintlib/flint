#include <stdio.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

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
