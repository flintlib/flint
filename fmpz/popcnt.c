#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

#ifdef POPCNT_INTRINSICS
static __inline__ mp_bitcnt_t shortCount(long val)
{
        return __builtin_popcountl(val);
}
#else
/* A naive implementation if neither your processor nor your compiler want to
 * do the work. */
static __inline__ mp_bitcnt_t shortCount(long val)
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
		if(mpz_cmp_si(t,0) < 0) return 0;
                else return mpz_popcount(t);
        }
}
