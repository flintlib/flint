/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

// XXX NOTE: the forwarding code was never completed, and neither was this
//           test file.

#include <iostream>
#include <sstream>

#include "flintxx/forwarding.h"

#include "flintxx/test/myint.h"
#include "flintxx/test/helpers.h"

using namespace flint;

template<class Operation, class Data>
class forwarded_expression
    : public expression<derived_wrapper<forwarded_expression>, Operation, Data>
{
public:
    forwarded_expression() {};
    template<class T>
    explicit forwarded_expression(const T& t)
        : expression<derived_wrapper< ::forwarded_expression>,
              Operation, Data>(t) {}

    template<class T>
    forwarded_expression& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

protected:
    explicit forwarded_expression(const Data& d)
        : expression<derived_wrapper< ::forwarded_expression>,
              Operation, Data>(d) {}

    template<class D, class O, class Da>
    friend class ::flint::expression;
};

struct hide_myint
{
    myint i;
    template<class T>
    hide_myint(const T& t) : i(t) {}
    hide_myint() {}
};

typedef forwarded_expression<operations::immediate, hide_myint> fwint;

namespace flint {
namespace forwarding {
template<>
struct enable<fwint>
    : mp::true_
{
    typedef myint underlying_t;
    static const myint& get_underlying(const fwint& fwd)
    {
        return fwd._data().i;
    }
    static myint& get_underlying(fwint& fwd)
    {
        return fwd._data().i;
    }
};

}
}

void
test_print()
{
    fwint fwd(4);
    std::ostringstream oss;
    oss << fwd;
    tassert(oss.str() == "4");
}

void
test_assignment()
{
    fwint f1, f2(4);

    f1 = WORD(3); // TODO understand why "3" seems to lead to circular dependency
    tassert(f1 == 3);

    f1 = f2;
    tassert(f1 == 4);
}

void
test_equals()
{
    fwint f1(4), f2(5);
    myint m1(4);

    tassert(f1 == f1);
    tassert(f1 != f2);
    tassert(f1 == 4);
    tassert(5 != f1);
    tassert(f1 == m1);
    tassert(m1 != f2);
}

void
test_conversion()
{
    fwint a(4);

    tassert(typed_equals(a.to<int>(), 4));
}

void
test_evaluation()
{
    fwint a(4);
    myint b(5);

    tassert(typed_equals((a + a).evaluate(), myint(8)));
    tassert(a + b == 9);
    tassert(a + 6 == 10);
    tassert((a + 6) + b == 15);
}

int
main()
{
    std::cout << "forwarding....";

    test_print();
    test_assignment();
    // TODO test cmp
    test_equals();
    test_conversion();
    test_evaluation();

    std::cout << "PASS" << std::endl;
    return 0;
}
