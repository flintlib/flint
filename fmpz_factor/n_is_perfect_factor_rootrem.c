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
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"


int
n_is_perfect_power(fmpz_t n, fmpz_t base, fmpz_t exponent)
{
    fmpz_t upper_power_limit, root, remainder, currbase;
    fmpz_init_set_ui(upper_power_limit, fmpz_flog_ui(n, 2) + 1);   /* upper limit on root */
    fmpz_init(root);
    fmpz_init(remainder);
    fmpz_init(currbase);
    fmpz_set(root, upper_power_limit);

    while (fmpz_cmp_ui(root, 1))    /*root from log2n to 2 */
    {
        int i = fmpz_rootrem_bsearch(remainder, currbase, n, root); /*checking if root rem = 1 */
        if (!fmpz_cmp_ui(remainder, 0)) 
        {
            fmpz_set(base, currbase);
            fmpz_set(exponent, root);
            return 1;
        }
        fmpz_sub_ui(root, root, 1);
    }

    fmpz_clear(upper_power_limit);
    fmpz_clear(root);
    fmpz_clear(remainder);
    fmpz_clear(currbase);

    return 0;

}