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
/******************************************************************************

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "arith.h"
#include "ulong_extras.h"

void check_value(len_t n, char *ans)
{
    fmpz_t x, y;
    fmpz_init(x);
    fmpz_init(y);
    fmpz_set_si(y, n);
    arith_ramanujan_tau(x, y);
    fmpz_set_str(y, ans, 10);
    if (!fmpz_equal(x,y))
    {
          printf("FAIL:\n");
          printf("tau(%ld) gave ", n);
          fmpz_print(x);
          printf(", expected %s\n", ans); 
          abort();
    }
    fmpz_clear(x);
    fmpz_clear(y);
}

void consistency_check(len_t n)
{
    fmpz_poly_t p;
    fmpz_t x, y;
    len_t k;

    fmpz_poly_init(p);
    fmpz_init(x);
    fmpz_init(y);

    arith_ramanujan_tau_series(p, n);
    if (p->length != n && !(n == 1 && p->length == 0))
    {
        printf("FAIL:\n");
        printf("wrong length of polynomial %ld\n", n);
        abort();
    }

    for (k=0; k<n; k++)
    {
        fmpz_set_si(y, k);
        arith_ramanujan_tau(x, y);
        fmpz_poly_get_coeff_fmpz(y, p, k);
        if (!fmpz_equal(x,y))
        {
            printf("FAIL:\n");
            printf("different tau n=%ld, k=%ld\n", n, k);
            fmpz_print(x);
            printf("\n");
            fmpz_print(y);
            printf("\n");
            abort();
        }
    }
    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_poly_clear(p);
}

int main(void)
{
    printf("ramanujan_tau....");
    fflush(stdout);

    check_value(0, "0");
    check_value(1, "1");
    check_value(2, "-24");
    check_value(3, "252");
    check_value(4, "-1472");
    check_value(5, "4830");
    check_value(6, "-6048");
    check_value(7, "-16744");
    check_value(8, "84480");
    check_value(9, "-113643");
    check_value(10, "-115920");
    check_value(11, "534612");
    check_value(12, "-370944");
    check_value(13, "-577738");
    check_value(14, "401856");
    check_value(15, "1217160");
    check_value(16, "987136");
    check_value(17, "-6905934");
    check_value(18, "2727432");
    check_value(19, "10661420");
    check_value(20, "-7109760");
    check_value(21, "-4219488");
    check_value(22, "-12830688");
    check_value(23, "18643272");
    check_value(24, "21288960");
    check_value(25, "-25499225");
    check_value(26, "13865712");
    check_value(27, "-73279080");
    check_value(28, "24647168");
    check_value(29, "128406630");
    check_value(30, "-29211840");
    check_value(31, "-52843168");
    check_value(32, "-196706304");
    check_value(33, "134722224");
    check_value(34, "165742416");
    check_value(35, "-80873520");
    check_value(36, "167282496");
    check_value(37, "-182213314");
    check_value(38, "-255874080");
    check_value(39, "-145589976");
    check_value(40, "408038400");
    check_value(41, "308120442");
    check_value(42, "101267712");
    check_value(43, "-17125708");
    check_value(44, "-786948864");
    check_value(45, "-548895690");
    check_value(46, "-447438528");
    check_value(47, "2687348496");
    check_value(48, "248758272");
    check_value(49, "-1696965207");

    check_value(1000, "-30328412970240000");
    check_value(10000, "-482606811957501440000");
    check_value(100000, "-2983637890141033828147200000");
    check_value(5040, "9072480147209256960");
    check_value(25401600, "-982963272212951631424865761586105548800");
    check_value(100003, "1194906306375914517502892252");
    check_value(15251, "-67392761749743476612496");
    check_value(16777216, "5141538908507386166920374725609506471936");
    check_value(43046721, "447670851294004737003138291024309833342241");
    check_value(462594208,
        "-313042078739616847874899392539635327193629261824");

    consistency_check(0);
    consistency_check(1);
    consistency_check(2);
    consistency_check(3);
    consistency_check(10);
    consistency_check(11);
    consistency_check(100);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
