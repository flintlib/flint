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

#include "cxx/test/helpers.h"
#include "cxx/mp.h"

using namespace flint;
using namespace mp;

void
test_equal_types()
{
    tassert((equal_types<int, int>::val));
    tassert((!equal_types<int, int&>::val));
}

void
test_logic()
{
    tassert(not_<true_>::val == false);
    tassert(not_<false_>::val == true);

    tassert((and_<true_, true_>::val == true));
    tassert((and_<false_, true_>::val == false));
    tassert((and_<false_, false_>::val == false));

    tassert((or_<true_, true_>::val == true));
    tassert((or_<false_, true_>::val == true));
    tassert((or_<false_, false_>::val == false));
}

template<class T>
typename enable_if<equal_types<T, int>, int>::type test_enable_if_1(T)
{
    return 0;
}
template<class T>
typename disable_if<equal_types<T, int>, int>::type test_enable_if_1(T)
{
    return 1;
}
template<class T>
int test_enable_if_2(T, typename enable_if<equal_types<T, int> >::type* = 0)
{
    return 0;
}
template<class T>
int test_enable_if_2(T, typename disable_if<equal_types<T, int> >::type* = 0)
{
    return 1;
}

void
test_enable_if()
{
    tassert(test_enable_if_1(int(1)) == 0);
    tassert(test_enable_if_1(unsigned(1)) == 1);
    tassert(test_enable_if_2(int(1)) == 0);
    tassert(test_enable_if_2(unsigned(1)) == 1);
}

int
main()
{
    std::cout << "mp....";

    test_equal_types();
    test_logic();
    test_enable_if();

    std::cout << "PASS" << std::endl;
    return 0;
}
