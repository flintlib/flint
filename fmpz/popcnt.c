#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

#ifdef POPCNT_INTRINSICS
static __inline__ mp_bitcnt_t shortCount(slong val)
{
#if defined(_WIN64)
   return __builtin_popcountll(val);  
#else
   return __builtin_popcountl(val);
#endif
}
#else
/* A naive implementation if neither your processor nor your compiler want to
 * do the work. */
static __inline__ mp_bitcnt_t shortCount(slong val)
{
        mp_bitcnt_t cnt;
        for(cnt=0; val; val >>= 1) {
                cnt += val & 1;
        }
        return cnt;
}
#endif

mp_bitcnt_t fmpz_popcnt(const fmpz_t c)
{
        fmpz c1;
        c1 = *c;
        if(!COEFF_IS_MPZ(c1))
        {
                if(*c < 0) return 0;
                else return shortCount(*c);
        } else
        {
                __mpz_struct *t = COEFF_TO_PTR(c1);
		if(flint_mpz_cmp_si(t,0) < 0) return 0;
                else return mpz_popcount(t);
        }
}
