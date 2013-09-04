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

    Copyright (C) 2013 Tom Bachmann

******************************************************************************/
#include <iostream>

#include "permxx.h"
#include "flintxx/test/helpers.h"

using namespace flint;

int
main()
{
    std::cout << "permxxxx....";

    permxx p1(10), p2(10);
    p1[0] = 1;p1[1] = 0;
    tassert(p1 != p2);
    p1 = p2;
    tassert(p1 == p2);
    permxx p3(p2);
    p2[1] = 0;
    p2[0] = 1;
    tassert(p1 != p2);
    tassert(p3 != p2);

    tassert(parity(p2) == 1);

    frandxx state;
    p1 = permxx::randtest(10, state);
    tassert(p1*inv(p1) == permxx::one(10));
    p2 = permxx::randtest(10, state);
    tassert(p1*p2 == compose(p1, p2));
    p3 = p1*p2;
    p1 *= p2;
    tassert(p1 == p3);

    tassert(maybe_perm_data(&p1) == p1._data());
    tassert(maybe_perm_data(0) == 0);

    if(0)
        print(p1); // make sure this compiles

    std::cout << "PASS" << std::endl;

    return 0;
}
