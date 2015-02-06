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
fmpz_is_perfect_power(fmpz_t n, fmpz_t base, fmpz_t exponent)
{
	if (fmpz_cmp_ui(n, 1)<=0)
	{
		fmpz_set_ui(base, -1);
		fmpz_set_ui(exponent, -1);
		return 0;
	}

	fmpz_t a, b, upper, lower, middle, a_upper_limit, power;
	unsigned long currpower = 0;

	fmpz_init(a);
	fmpz_init(b);
	fmpz_init(upper);
	fmpz_init(middle);
	fmpz_init(lower);
	fmpz_init(a_upper_limit);	
	fmpz_init(power);

	fmpz_sqrt(a_upper_limit, n);				/* Set upper limit on a = n^1/2 */
	long power_upper_limit = fmpz_flog_ui(n, 2) + 1;	/* Set upper limit on b */
	fmpz_init_set_ui(b, power_upper_limit);
		
	/* checking whether an a exists for a given b */

	while (fmpz_cmp_ui(b, 1) != 0)				/* loop b from log2(a) to 2 */ 					
	{
		fmpz_set_ui(lower, 2);
		*upper = *a_upper_limit;
		currpower = fmpz_get_ui(b);
		while (fmpz_cmp(lower, upper) <= 0)		/* implementing binary serch */	
		{
			fmpz_add(middle, lower, upper);		/* find middle */
			fmpz_fdiv_q_ui(middle, middle, 2);
			fmpz_pow_ui(power, middle, currpower);				
			
			if (!fmpz_cmp(power, n))
			{
				*exponent = *b;
				*base = *middle;
				return 1;
			}	
			else if (fmpz_cmp(power, n)<0)
				fmpz_add_ui(lower, middle, 1);
			else
				fmpz_sub_ui(upper, middle, 1);	
		}
		fmpz_sub_ui(b, b, 1);				/* decrement b */
	}
	fmpz_clear(b);
	fmpz_clear(a);
	fmpz_clear(upper);
	fmpz_clear(middle);
	fmpz_clear(lower);
	fmpz_clear(a_upper_limit);
	fmpz_clear(power);

	fmpz_set_ui(base, -1);
	fmpz_set_ui(exponent, -1);
	return 0;
}