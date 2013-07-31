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

#include "flintxx/ltuple.h"
#include "fmpzxx.h"

#include "flintxx/test/helpers.h"

using namespace flint;

void
test_ltdetail()
{
    using namespace mp;

    tassert((mp::equal_types<int&, ltdetail::to_ref<int>::type>::val));
    tassert((mp::equal_types<fmpzxx_ref, ltdetail::to_ref<fmpzxx>::type>::val));
    tassert((mp::equal_types<const int&, ltdetail::to_srcref<int>::type>::val));
    tassert((mp::equal_types<fmpzxx_srcref,
                ltdetail::to_srcref<fmpzxx>::type>::val));

    tassert((mp::equal_types<make_tuple<int&, fmpzxx_ref>::type,
                ltdetail::ref_tuple<make_tuple<int, fmpzxx>::type>::type>::val));
    tassert((mp::equal_types<make_tuple<const int&, fmpzxx_srcref>::type,
                ltdetail::srcref_tuple<make_tuple<int, fmpzxx>::type>::type>::val));
}

void
test_traits()
{
    typedef mp::make_tuple<fmpzxx, int>::type tuple_t;
    typedef make_ltuple<tuple_t>::type ltuple_t;
    typedef make_ltuple<tuple_t>::ref_type ltuple_ref_t;
    typedef make_ltuple<tuple_t>::srcref_type ltuple_srcref_t;

    tassert(traits::is_ltuple_source<ltuple_t>::val);
    tassert(traits::is_ltuple_source<ltuple_ref_t>::val);
    tassert(traits::is_ltuple_source<ltuple_srcref_t>::val);
    tassert(traits::is_ltuple_target<ltuple_t>::val);
    tassert(traits::is_ltuple_target<ltuple_ref_t>::val);
    tassert(!traits::is_ltuple_target<ltuple_srcref_t>::val);
    tassert(!traits::is_ltuple_source<int>::val);
    tassert(!traits::is_ltuple_expr<int>::val);
    tassert(!traits::is_ltuple_expr<fmpzxx>::val);
}

namespace flint {
typedef make_ltuple<mp::make_tuple<fmpzxx, int>::type>::type fmpzxxint_pair;
FLINT_DEFINE_BINOP(make_lazy_test)
namespace rules {
FLINT_DEFINE_BINARY_EXPR_COND2(make_lazy_test_op, fmpzxxint_pair,
        FMPZXX_COND_S, traits::is_signed_integer,
        to._data().inner.head = e1;to._data().inner.tail.head = e2)
}
}

void
test_equals()
{
    typedef mp::make_tuple<fmpzxx, int> maker;
    typedef mp::make_tuple<fmpzxx_ref, int&> refmaker;
    typedef mp::make_tuple<fmpzxx_srcref, const int&> srcrefmaker;
    typedef make_ltuple<maker::type> lmaker;

    fmpzxx f;
    int a;
    lmaker::type ltuple(detail::INSTANTIATE_FROM_TUPLE(),
            maker::make(fmpzxx(1), 2));
    lmaker::ref_type lref(lmaker::ref_type::_make(refmaker::make(f, a)));
    lmaker::srcref_type lsrcref(lmaker::srcref_type::_make(srcrefmaker::make(
                    ltuple._data().inner.head, ltuple._data().inner.tail.head)));

    tassert(ltuple == ltuple);
    tassert(ltuple != lref);
    tassert(ltuple == lsrcref);
    tassert(lref != lsrcref);
    f = 1;
    a = 2;
    tassert(ltuple == lref);
    tassert(lref == lsrcref);

    tassert(ltuple == make_lazy_test(fmpzxx(1), 2));
}

void
test_assignment()
{
    typedef mp::make_tuple<fmpzxx, int> maker;
    typedef mp::make_tuple<fmpzxx_ref, int&> refmaker;
    typedef mp::make_tuple<fmpzxx_srcref, const int&> srcrefmaker;
    typedef make_ltuple<maker::type> lmaker;

    fmpzxx f;
    int a;
    lmaker::type ltuple(detail::INSTANTIATE_FROM_TUPLE(),
            maker::make(fmpzxx(1), 2));
    lmaker::ref_type lref(lmaker::ref_type::_make(refmaker::make(f, a)));
    lref = ltuple;
    tassert(f == 1 && a == 2);

    f = 0;
    a = 0;
    lmaker::srcref_type lsrcref(lmaker::srcref_type::_make(srcrefmaker::make(
                    ltuple._data().inner.head, ltuple._data().inner.tail.head)));
    ltuple._data().inner.head = 17;
    lref = lsrcref;
    tassert(f == 17 && a == 2);

    f = 3;
    a = 4;
    ltuple = lref;
    tassert(ltuple._data().inner.head == 3
            && ltuple._data().inner.tail.head == 4);

    lref = make_lazy_test(fmpzxx(17), 18);
    tassert(f == 17 && a == 18);
}

void
test_ltupleref()
{
    fmpzxx a, b;
    int c;
    ltupleref(c) = ltuple(2);
    tassert(c == 2);
    ltupleref(a, c) = ltuple(fmpzxx(3), 4);
    tassert(a == 3 && c == 4);
}

int
main()
{
    std::cout << "ltuple....";

    test_ltdetail();
    test_traits();
    test_equals();
    test_assignment();
    test_ltupleref();

    std::cout << "PASS" << std::endl;

    return 0;
}
