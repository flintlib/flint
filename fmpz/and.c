#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void fmpz_and(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
        fmpz_t tmp;
        fmpz c1,c2;
        c1 = *g;
        c2 = *h;
        if(!COEFF_IS_MPZ(c1))
        {
            if(!COEFF_IS_MPZ(c2)) /* both inputs are small */
            {
                fmpz_set_si(tmp, c1 & c2);
            } else /* g is small, h is large */
            {
                mpz_t g2;
                __mpz_struct * mpz3 = _fmpz_promote(tmp);
                __mpz_struct * mpz2 = COEFF_TO_PTR(c2);
                mpz_init_set_si(g2, *g);
                mpz_and(mpz3, mpz2, g2);
                mpz_clear(g2);
                _fmpz_demote_val(tmp);
            }
        } else
        {
            if(!COEFF_IS_MPZ(c2)) /* g is large, h is small */
            {
                mpz_t h2;
                __mpz_struct *mpz3 = _fmpz_promote(tmp);
                __mpz_struct *mpz1 = COEFF_TO_PTR(c1);
                mpz_init_set_si(h2, *h);
                mpz_and(mpz3, mpz1, h2);
                mpz_clear(h2);
                _fmpz_demote_val(tmp);
            } else /* g and h are large */
            {
                __mpz_struct * mpz3 = _fmpz_promote(tmp);
                __mpz_struct * mpz1 = COEFF_TO_PTR(c1);
                __mpz_struct * mpz2 = COEFF_TO_PTR(c2);
                mpz_and(mpz3, mpz1, mpz2);
                _fmpz_demote_val(tmp);
            }
        }
        fmpz_set(f,tmp);
        fmpz_clear(tmp);
}
