#include <stdio.h>
#include <stdlib.h>
#include <flint.h>
#include <ulong_extras.h>

mp_limb_t n_primitive_root_general_prefactor(mp_limb_t n, n_factor_t * factors)
{
    slong i;
    int found;
    mp_limb_t result, a, pm1;
    double pinv;

    pm1 = n_euler_phi(n);
    pinv = n_precompute_inverse(n);

    for (a = 2; a < n; a++)
    {
        if(n_gcd_full(a, n) == 1)
        {
            found = 1;
            for (i = 0; i < factors->num; i++)
            {
                result = n_powmod_precomp(a, pm1 / factors->p[i], n, pinv);
                if (result == 1)
                {
                    found = 0;
                    break;
                }
            }
            if (found)
            {
                return a;
            }
        }
    }
    flint_printf("Exception (n_primitive_root_prime_prefactor).  root not found.\n");
    abort();
}

mp_limb_t n_primitive_root_general(mp_limb_t n)
{
    if(n_is_prime(n)) return n_primitive_root_prime(n);
    else
    {
        if (n == 4) return 3;
        else if(n % 4 == 0)
        {
            flint_printf("Exception (n_primitive_root_prime_prefactor).  root not found.\n");
            return;
        }
        else
        {
            mp_limb_t a, m;
            n_factor_t factors;
            n_factor_init(&factors);
            if(n % 2 == 0)
            {
                m = 2;
                n /= 2;
            }
            else
            {
                m = 1;
            }
            n_factor(&factors, n, 1);
            if(factors.num >= 2)
            {
                flint_printf("Exception (n_primitive_root_prime_prefactor).  root not found.\n");
                return;
            }
            else
            {
                n *= m;
                n_factor_init(&factors);
                n_factor(&factors, n_euler_phi(n), 1);
                a = n_primitive_root_general_prefactor(n, &factors);
                return a;
            }
        }
    }
}
