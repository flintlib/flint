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

inline double absolute(double x) { return x >= 0 ? x : -x; }

int
fmpz_rootrem_newton_iteration(fmpz_t remainder, fmpz_t base, fmpz_t n, fmpz_t root)
{
    if (fmpz_cmp_ui(n, 1)<=0)
        return 0;

    if (fmpz_cmp_ui(root, 1)<=0)
        return 0;

    double d, x, a, r;
    d = 0;
    x = 1;
    a = *n;
    r = *root;

    if (r == 2)
    {
        do {
            d = ((a / x) - x) / 2;
            x += d;
        } while (absolute(d)>=(absolute(x)*(small_float)));

    }
    else
    {
        do {
            d = (a / pow(x, r-1) - x) / r;
            x += d;
        } while (absolute(d)>=(absolute(x)*(small_float)));
    }

    *remainder = x;
    fmpz_set(base, remainder);

    fmpz_pow_ui(remainder, remainder, fmpz_get_ui(root));
    fmpz_sub(remainder, n, remainder);
    return 1;
}