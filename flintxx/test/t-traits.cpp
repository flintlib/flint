/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flintxx/test/helpers.h"
#include "flintxx/traits.h"
#include "flintxx/mp.h"

using namespace flint;
using namespace traits;

struct newtype { };

void
test_integrality()
{
    tassert(is_signed_integer<signed char>::val == true);
    tassert(is_signed_integer<signed int>::val == true);
    tassert(is_signed_integer<slong>::val == true);
    tassert(is_signed_integer<unsigned char>::val == false);
    tassert(is_signed_integer<unsigned int>::val == false);
    tassert(is_signed_integer<ulong>::val == false);
    tassert(is_signed_integer<newtype>::val == false);

    tassert(is_unsigned_integer<unsigned char>::val == true);
    tassert(is_unsigned_integer<unsigned int>::val == true);
    tassert(is_unsigned_integer<ulong>::val == true);
    tassert(is_unsigned_integer<signed char>::val == false);
    tassert(is_unsigned_integer<signed int>::val == false);
    tassert(is_unsigned_integer<slong>::val == false);
    tassert(is_unsigned_integer<newtype>::val == false);

    tassert(is_integer<unsigned char>::val == true);
    tassert(is_integer<unsigned int>::val == true);
    tassert(is_integer<ulong>::val == true);
    tassert(is_integer<signed char>::val == true);
    tassert(is_integer<signed int>::val == true);
    tassert(is_integer<slong>::val == true);
    tassert(is_integer<newtype>::val == false);

    tassert(is_string<int>::val == false);
    tassert(is_string<const char*>::val == true);
    tassert(is_string<char*>::val == true);
    tassert(is_string<char[5]>::val == true);

    tassert(fits_into_slong<unsigned char>::val == true);
    tassert(fits_into_slong<unsigned short>::val == true);
    tassert(fits_into_slong<ulong>::val == false);
    tassert(fits_into_slong<signed char>::val == true);
    tassert(fits_into_slong<signed short>::val == true);
    tassert(fits_into_slong<signed int>::val == true);
    tassert(fits_into_slong<slong>::val == true);
    tassert(fits_into_slong<newtype>::val == false);

    typedef void function_type(int, int);
    tassert(fits_into_slong<function_type>::val == false);
}

void
test_manipulation()
{
    using mp::equal_types;

    tassert((equal_types<forwarding<int>::type, int>::val));
    tassert((equal_types<forwarding<newtype&>::type, newtype&>::val));
    tassert((equal_types<forwarding<const newtype>::type,
                const newtype&>::val));
    tassert((equal_types<forwarding<const newtype&>::type,
                const newtype&>::val));

    tassert((equal_types<reference<int>::type, int&>::val));
    tassert((equal_types<reference<newtype&>::type, newtype&>::val));
    tassert((equal_types<reference<const slong>::type, const slong&>::val));
    tassert((equal_types<reference<const int&>::type, const int&>::val));

    tassert((equal_types<make_const<int>::type, const int>::val));
    tassert((equal_types<make_const<newtype&>::type, const newtype&>::val));
    tassert((equal_types<make_const<const slong>::type, const slong>::val));
    tassert((equal_types<make_const<const int&>::type, const int&>::val));

    tassert((equal_types<basetype<int>::type, int>::val));
    tassert((equal_types<basetype<newtype&>::type, newtype>::val));
    tassert((equal_types<basetype<const slong>::type, slong>::val));
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

    // Test the HACK.
    tassert((_is_convertible<int, void(int)>::val == false));
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
