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

#ifndef CXX_TEST_MYINT_H
#define CXX_TEST_MYINT_H CXX_TEST_MYINT_H

#include "cxx/expression.h"

namespace flint {
template<class Operation, class Data>
class my_expression
    : public expression<derived_wrapper<my_expression>, Operation, Data>
{
public:
    my_expression() {};
    template<class T>
    explicit my_expression(const T& t)
        : expression<derived_wrapper< ::flint::my_expression>,
              Operation, Data>(t) {}

    template<class T>
    my_expression& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

protected:
    explicit my_expression(const Data& d)
        : expression<derived_wrapper< ::flint::my_expression>,
              Operation, Data>(d) {}

    template<class D, class O, class Da>
    friend class flint::expression;
};

struct data
{
    int payload;
    bool* destroyed;
    int extra;

    data() : payload(-1), destroyed(0), extra(42) {}

    ~data()
    {
        if(destroyed)
            *destroyed = true;
    }

    data(const data& d)
        : payload(d.payload), destroyed(0), extra(1) {}
    data(int i)
        : payload(i), destroyed(0), extra(2) {}
    data(char i)
        : payload(i), destroyed(0), extra(3) {}
};

typedef my_expression<operations::immediate, data> myint;

template<class Operation, class Data>
class my_expression2
    : public expression<derived_wrapper<my_expression2>, Operation, Data>
{
public:
    // cannot have a default constructor

    typedef expression<derived_wrapper< ::flint::my_expression2>,
              Operation, Data> base_t;
    typedef typename base_t::evaluated_t evaluated_t;

    template<class T>
    explicit my_expression2(const T& t)
        : base_t(t) {}

    template<class T>
    explicit my_expression2(T& t)
        : base_t(t) {}

    template<class T>
    my_expression2& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

    evaluated_t create_temporary() const
    {
        return evaluated_t(0l);
    }

protected:
    explicit my_expression2(const Data& d)
        : base_t(d) {}

    template<class D, class O, class Da>
    friend class flint::expression;
};

struct longref_data;
struct longcref_data;
struct long_data;
typedef my_expression2<operations::immediate, long_data> mylong;
typedef my_expression2<operations::immediate, longref_data> mylong_ref;
typedef my_expression2<operations::immediate, longcref_data> mylong_cref;

struct long_data
{
    long payload;
    // no default constructor
    long_data(long d) : payload(d) {}
    long_data(const myint& m) : payload(m._data().payload) {}

    long_data(const mylong_ref&);
    long_data(const mylong_cref&);
};

struct longref_data
{
    long& payload;

    longref_data(mylong& l) : payload(l._data().payload) {}
};

struct longcref_data
{
    const long& payload;

    longcref_data(const mylong& l) : payload(l._data().payload) {}
    longcref_data(mylong_ref lr) : payload(lr._data().payload) {}
};

inline long_data::long_data(const mylong_ref& mlr) : payload(mlr._data().payload) {}
inline long_data::long_data(const mylong_cref& mlr) : payload(mlr._data().payload) {}

namespace mylong_traits {
template<class T> struct is_source : mp::false_ { };
template<class T> struct is_target : mp::false_ { };

template<> struct is_source<mylong> : mp::true_ { };
template<> struct is_source<mylong_ref> : mp::true_ { };
template<> struct is_source<mylong_cref> : mp::true_ { };
template<> struct is_target<mylong> : mp::true_ { };
template<> struct is_target<mylong_ref> : mp::true_ { };
}

namespace rules {

template<>
struct print<myint>
{
    static void doit(const myint& v, std::ostream& o)
    {
        o << v._data().payload;
    }
};

template<>
struct assignment<myint, myint>
{
    static void doit(myint& to, const myint& from)
    {
        to._data().payload = from._data().payload;
        to._data().extra = 4;
    }
};

template<>
struct assignment<myint, long>
{
    static void doit(myint& to, long from)
    {
        to._data().payload = from;
        to._data().extra = 5;
    }
};

template<>
struct equals<myint, myint>
{
    static bool get(const myint& i1, const myint& i2)
    {
        return i1._data().payload == i2._data().payload;
    }
};

template<>
struct equals<myint, int>
{
    static bool get(const myint& i1, int i2)
    {
        return i1._data().payload == i2;
    }
};

template<>
struct conversion<int, myint>
{
    static int get(const myint& from)
    {
        return from._data().payload;
    }
};

template<>
struct commutative_binary_expression<myint, operations::plus, myint>
{
    typedef myint return_t;
    static void doit(myint& to, const myint& a1, const myint& a2)
    {
        to._data().payload = a1._data().payload + a2._data().payload;
    }
};

template<class T>
struct commutative_binary_expression<myint,
    typename mp::enable_if<traits::is_integer<T>, operations::plus>::type,
    T>
{
    typedef myint return_t;
    static void doit(myint& to, const myint& a1, T a2)
    {
        to._data().payload = a1._data().payload + a2;
    }
};

template<>
struct commutative_binary_expression<myint, operations::times, myint>
{
    typedef myint return_t;
    static void doit(myint& to, const myint& a1, const myint& a2)
    {
        to._data().payload = a1._data().payload * a2._data().payload;
    }
};

template<>
struct binary_expression<myint, operations::minus, myint>
{
    typedef myint return_t;
    static void doit(myint& to, const myint& a1, const myint& a2)
    {
        to._data().payload = a1._data().payload - a2._data().payload;
    }
};

template<>
struct binary_expression<myint, operations::divided_by, myint>
{
    typedef myint return_t;
    static void doit(myint& to, const myint& a1, const myint& a2)
    {
        to._data().payload = a1._data().payload / a2._data().payload;
    }
};

template<>
struct binary_expression<myint, operations::modulo, myint>
{
    typedef myint return_t;
    static void doit(myint& to, const myint& a1, const myint& a2)
    {
        to._data().payload = a1._data().payload % a2._data().payload;
    }
};

template<>
struct binary_expression<myint, operations::shift, int>
{
    typedef myint return_t;
    static void doit(myint& to, const myint& a1, int a2)
    {
        if (a2 >= 0)
            to._data().payload = a1._data().payload << a2;
        else
            to._data().payload = a1._data().payload >> (-a2);
    }
};

template<>
struct unary_expression<operations::negate, myint>
{
    typedef myint return_t;
    static void doit(myint& to, const myint& from)
    {
        to._data().payload = - from._data().payload;
    }
};


/////////////////////////////////////////////////////////////////////////////
// Minimal rules for mylong
/////////////////////////////////////////////////////////////////////////////

template<class T, class U>
struct equals<T, U, typename mp::enable_if<mp::and_<
    mylong_traits::is_source<T>, mylong_traits::is_source<U> > >::type>
{
    static bool get(const T& i1, const U& i2)
    {
        return i1._data().payload == i2._data().payload;
    }
};

template<class T>
struct equals<T, long,
    typename mp::enable_if<mylong_traits::is_source<T> >::type>
{
    static bool get(const T& i1, long i2)
    {
        return i1._data().payload == i2;
    }
};

#if 0
template<bool c, class Op, class Data>
struct evaluation<
    my_expression<
        operations::plus,
        make_tuple<const mylong&, const mylong&>::type>,
    Op, Data,
    c, 0>
{
    typedef mylong return_t;
    typedef empty_tuple temporaries_t;
    typedef my_expression<
        operations::plus,
        make_tuple<const mylong&, const mylong&>::type> expr_t;
    static void doit(const expr_t& input, temporaries_t temps, return_t* output)
    {
        output->_data().payload = input._data().first()._data().payload
                                + input._data().second()._data().payload;
    }
};
#endif

template<class T, class U>
struct commutative_binary_expression<T, typename mp::enable_if<mp::and_<
        mylong_traits::is_source<T>,
        mylong_traits::is_source<U> >,
    operations::plus>::type, U>
{
    typedef mylong return_t;

    template<class V>
    static void doit(V& to, const T& a1, const U& a2)
    {
        to._data().payload = a1._data().payload + a2._data().payload;
    }
};

template<class U>
struct commutative_binary_expression<typename mp::enable_if<
    mylong_traits::is_source<U>, myint>::type, operations::plus, U>
{
    typedef mylong return_t;

    template<class V>
    static void doit(V& to, const myint& a1, const U& a2)
    {
        to._data().payload = a1._data().payload + a2._data().payload;
    }
};

template<class T, class U>
struct assignment<T, U, typename mp::enable_if<mp::and_<
    mylong_traits::is_target<T>, mylong_traits::is_source<U> > >::type>
{
    static void doit(T& to, const U& from)
    {
        to._data().payload = from._data().payload;
    }
};

template<class T>
struct assignment<T, long,
    typename mp::enable_if<mylong_traits::is_target<T> >::type>
{
    static void doit(T& to, long from)
    {
        to._data().payload = from;
    }
};
} // rules
} // flint

#endif
