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
#include "cxx/traits.h"
#include "cxx/mp.h"

using namespace flint;
using namespace traits;

struct newtype { };

void
test_integrality()
{
    tassert(is_signed_integer<signed char>::val == true);
    tassert(is_signed_integer<signed int>::val == true);
    tassert(is_signed_integer<signed long>::val == true);
    tassert(is_signed_integer<unsigned char>::val == false);
    tassert(is_signed_integer<unsigned int>::val == false);
    tassert(is_signed_integer<unsigned long>::val == false);
    tassert(is_signed_integer<newtype>::val == false);

    tassert(is_unsigned_integer<unsigned char>::val == true);
    tassert(is_unsigned_integer<unsigned int>::val == true);
    tassert(is_unsigned_integer<unsigned long>::val == true);
    tassert(is_unsigned_integer<signed char>::val == false);
    tassert(is_unsigned_integer<signed int>::val == false);
    tassert(is_unsigned_integer<signed long>::val == false);
    tassert(is_unsigned_integer<newtype>::val == false);
}

void
test_manipulation()
{
    using mp::equal_types;

    tassert((equal_types<forwarding<int>::type, const int&>::val));
    tassert((equal_types<forwarding<newtype&>::type, const newtype&>::val));
    tassert((equal_types<forwarding<const long>::type, const long&>::val));
    tassert((equal_types<forwarding<const int&>::type, const int&>::val));

    tassert((equal_types<reference<int>::type, int&>::val));
    tassert((equal_types<reference<newtype&>::type, newtype&>::val));
    tassert((equal_types<reference<const long>::type, const long&>::val));
    tassert((equal_types<reference<const int&>::type, const int&>::val));

    tassert((equal_types<make_const<int>::type, const int>::val));
    tassert((equal_types<make_const<newtype&>::type, const newtype&>::val));
    tassert((equal_types<make_const<const long>::type, const long>::val));
    tassert((equal_types<make_const<const int&>::type, const int&>::val));

    tassert((equal_types<basetype<int>::type, int>::val));
    tassert((equal_types<basetype<newtype&>::type, newtype>::val));
    tassert((equal_types<basetype<const long>::type, long>::val));
    tassert((equal_types<basetype<const int&>::type, int>::val));
    tassert((equal_types<basetype<int*>::type, int*>::val));
}

class super { };
class sub : public super { };

void
test_convertibility()
{
    tassert((_is_convertible<super, sub>::val == true));
    tassert((_is_convertible<sub, super>::val == false));
    tassert((_is_convertible<sub, newtype>::val == false));
    tassert((_is_convertible<sub, sub>::val == true));
    tassert((_is_convertible<int, int>::val == true));
}

int
main()
{
    std::cout << "traits....";

    test_integrality();
    test_manipulation();
    test_convertibility();

    std::cout << "PASS" << std::endl;
    return 0;
}
