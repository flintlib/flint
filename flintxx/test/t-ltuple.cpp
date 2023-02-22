/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flintxx/ltuple.h"
#include "fmpzxx.h"

#include "flintxx/test/helpers.h"
#include "flintxx/test/myint.h"

using namespace flint;

void
test_traits()
{
    tassert(!traits::is_ltuple_expr<int>::val);
    tassert(!traits::is_ltuple_expr<fmpzxx>::val);
}

namespace flint {
typedef make_ltuple<mp::make_tuple<fmpzxx, int>::type>::type fmpzxxint_pair;
typedef make_ltuple<mp::make_tuple<mylong, int>::type>::type mylongint_pair;
FLINT_DEFINE_BINOP(make_lazy_test)
namespace rules {
FLINT_DEFINE_BINARY_EXPR_COND2(make_lazy_test_op, fmpzxxint_pair,
        FMPZXX_COND_S, traits::is_signed_integer,
        to.template get<0>() = e1;to.template get<1>() = e2)
FLINT_DEFINE_BINARY_EXPR2(make_lazy_test_op, mylongint_pair,
        mylong, int,
        to.template get<0>() = e1;to.template get<1>() = e2)
}
}

void
test_equals()
{
    typedef mp::make_tuple<fmpzxx, int> maker;
    typedef mp::make_tuple<fmpzxx&, int&> refmaker;
    typedef mp::make_tuple<const fmpzxx&, const int&> srcrefmaker;
    typedef make_ltuple<maker::type> lmaker;
    typedef make_ltuple<refmaker::type> lrefmaker;
    typedef make_ltuple<srcrefmaker::type> lsrcrefmaker;

    fmpzxx f;
    int a = 12345;
    lmaker::type ltuple(detail::INSTANTIATE_FROM_TUPLE(),
            maker::make(fmpzxx(1), 2));
    lrefmaker::ref_type lref(detail::INSTANTIATE_FROM_TUPLE(), refmaker::make(f, a));
    lsrcrefmaker::srcref_type lsrcref(detail::INSTANTIATE_FROM_TUPLE(), srcrefmaker::make(
                    ltuple._data().inner.head, ltuple._data().inner.tail.head));

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
    typedef mp::make_tuple<fmpzxx&, int&> refmaker;
    typedef mp::make_tuple<const fmpzxx&, const int&> srcrefmaker;
    typedef make_ltuple<maker::type> lmaker;

    fmpzxx f;
    int a;
    lmaker::type ltuple(detail::INSTANTIATE_FROM_TUPLE(),
            maker::make(fmpzxx(1), 2));
    lmaker::ref_type lref(detail::INSTANTIATE_FROM_TUPLE(), refmaker::make(f, a));
    lref = ltuple;
    tassert(f == 1 && a == 2);

    f = 0;
    a = 0;
    lmaker::srcref_type lsrcref(detail::INSTANTIATE_FROM_TUPLE(), srcrefmaker::make(
                    ltuple._data().inner.head, ltuple._data().inner.tail.head));
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

    // test assignment with type conversion
    ltupleref(a, b, c) = ltuple(1, 2, 4u);
    tassert(a == 1 && b == 2 && c == 4);

    // test assignment with c-style references
    ltuple(fmpzxx_ref(a)) = ltuple(b);
    tassert(a == 2);
}

template<class T>
bool is_lazy(const T&)
{
    return traits::is_lazy_expr<T>::val;
}
void
test_get()
{
    fmpzxx a(1); int b = 2;
    tassert(ltuple(a, b).get<0>() == a && ltuple(a, b).get<1>() == b);
    ltupleref(a, b).get<0>() = 17;
    ltupleref(a, b).get<1>() = 15;
    tassert(a == 17 && b == 15);
    tassert(make_lazy_test(a, 3).get<0>() == 17);
    tassert(is_lazy(make_lazy_test(a, 3).get<0>()));
    tassert(make_lazy_test(a, 3).get<1>() == 3);
}

void
test_placeholder()
{
    fmpzxx a; int b;
    ltupleref(a, _) = make_lazy_test(fmpzxx(17), 5);
    tassert(a == 17);
    ltupleref(_, b) = make_lazy_test(fmpzxx(18), 6);
    tassert(b == 6 && a == 17);
    ltupleref(_, _) = make_lazy_test(fmpzxx(1), 2);
    tassert(b == 6 && a == 17);
}

void
test_create_temporaries()
{
    tassert(make_lazy_test(mylong(1), 2).get<0>() == mylong(1));
}

int
main()
{
    std::cout << "ltuple....";

    test_traits();
    test_equals();
    test_assignment();
    test_ltupleref();
    test_get();
    test_placeholder();
    test_create_temporaries();

    std::cout << "PASS" << std::endl;

    return 0;
}
