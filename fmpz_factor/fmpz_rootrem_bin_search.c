/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/


#define ulong ulongxx /* interferes with system includes */
#include <string.h>
#include <math.h>
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

int
fmpz_rootrem_bsearch(fmpz_t remainder, fmpz_t base, fmpz_t n, fmpz_t root) /* root >= 2 */
{
    if (fmpz_cmp_ui(n, 1)<=0)
        return 0;

    if (fmpz_cmp_ui(root, 1)<=0)
        return 0;

    fmpz_t guess, step, power, sum;
    int comp;
    unsigned long int sum_ui;

    fmpz_init(power);
    fmpz_init_set_ui(step, 1);
    fmpz_init_set_ui(guess, 1);
    fmpz_init(sum);

    while (1)
    {            
        fmpz_add(sum, step, guess);                 /* sum = step + guess */ 
        fmpz_pow_ui(power, sum, fmpz_get_ui(root)); /* power = sum ** n */
        comp = fmpz_cmp(power, n);

        if (!comp)
        {
            fmpz_set_ui(remainder, 0);
            fmpz_set(base, sum);
            goto cleanup;
        }
        else if (comp < 0)
            fmpz_mul_ui(step, step, 2);

        else if (!fmpz_cmp_ui(step, 1))
        {
            fmpz_pow_ui(power, guess, fmpz_get_ui(root));
            fmpz_set(base, guess);
            fmpz_sub(remainder, n, power);
            goto cleanup;
        }
        else
        {
            fmpz_fdiv_q_ui(step, step, 2);
            fmpz_add(guess, guess, step);
            fmpz_set_ui(step, 1);
        }
    }
    cleanup:

    fmpz_clear(guess);
    fmpz_clear(step);
    fmpz_clear(power);
    fmpz_clear(sum);

    return 1;
}
