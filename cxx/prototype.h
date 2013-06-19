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

// TODO
// * formatting
// * static asserts
// * contexts
// * tuple merging

#ifndef CXX_PROTOTYPE_H
#define CXX_PROTOTYPE_H

#include <iosfwd>
#include <cstdlib>

namespace flint {

namespace traits {
struct true_
{
    static const bool val = true;
};

struct false_
{
    static const bool val = false;
};

template<class T> struct is_signed_integer : false_ { };
template<> struct is_signed_integer<signed char> : true_ { };
template<> struct is_signed_integer<signed short> : true_ { };
template<> struct is_signed_integer<signed int> : true_ { };
template<> struct is_signed_integer<signed long> : true_ { };

template<class T> struct is_unsigned_integer : false_ { };
template<> struct is_unsigned_integer<unsigned char> : true_ { };
template<> struct is_unsigned_integer<unsigned short> : true_ { };
template<> struct is_unsigned_integer<unsigned int> : true_ { };
template<> struct is_unsigned_integer<unsigned long> : true_ { };

template<class T> struct forwarding {typedef const T& type;};
template<class T> struct forwarding<T&> {typedef const T& type;};
template<class T> struct forwarding<const T&> {typedef const T& type;};
template<class T> struct reference {typedef T& type;};
template<class T> struct reference<T&> {typedef T& type;};
template<class T> struct reference<const T&> {typedef const T& type;};
template<class T> struct make_const {typedef const T type;};
template<class T> struct make_const<T&> {typedef const T& type;};
template<class T> struct basetype {typedef T type;};
template<class T> struct basetype<T&> {typedef T type;};
template<class T> struct basetype<const T&> {typedef T type;};
}

template<class Head, class Tail>
struct tuple
{
    Head head;
    Tail tail;

    typedef Head head_t;

    typename traits::reference<head_t>::type first() {return head;}
    typename traits::forwarding<head_t>::type first() const {return head;}
    typename traits::reference<typename Tail::head_t>::type second() {return tail.head;}
    typename traits::forwarding<typename Tail::head_t>::type second() const {return tail.head;}

    tuple() {};
    tuple(typename traits::forwarding<Head>::type h,
          typename traits::forwarding<Tail>::type t)
        : head(h), tail(t)
    {
    }
};
struct empty {};

template<>
struct tuple<void, void>
{
    typedef empty head_t;
    typedef empty tail_t;
    empty head;
    empty tail;
};
typedef tuple<void, void> empty_tuple;

namespace mp {
template<class T1, class T2 = void, class T3 = void>
struct make_tuple
{
    typedef tuple<T1, tuple<T2, tuple<T3, empty_tuple> > > type;
    static type make(typename traits::forwarding<T1>::type t1,
              typename traits::forwarding<T2>::type t2,
              typename traits::forwarding<T3>::type t3)
    {
        return type(t1, make_tuple<T2, T3>::make(t2, t3));
    }
};
template<class T1, class T2>
struct make_tuple<T1, T2, void>
{
    typedef tuple<T1, tuple<T2, empty_tuple> > type;
    static type make(typename traits::forwarding<T1>::type t1,
              typename traits::forwarding<T2>::type t2)
    {
        return type(t1, make_tuple<T2>::make(t2));
    }
};
template<class T1>
struct make_tuple<T1, void, void>
{
    typedef tuple<T1, empty_tuple> type;
    static type make(typename traits::forwarding<T1>::type t1)
    {
        return type(t1, empty_tuple());
    }
};

// Create a tuple backing a tuple of points.
// If Tuple::head_t is the same as Return*, do not back the head
// and instead feed it in separately.
template<class Tuple, class Return = void>
struct back_tuple;

template<class T, class Return>
struct back_tuple<tuple<T*, empty_tuple>, Return>
{
    typedef tuple<T, empty_tuple> type;
    static void init(tuple<T*, empty_tuple>& to, type& from, Return* ret /* unused */)
    {
        to.head = &from.head;
    }
};
template<class Head, class Tail, class Return>
struct back_tuple<tuple<Head*, Tail>, Return>
{
    typedef tuple<Head, typename back_tuple<Tail, empty_tuple>::type> type;
    static void init(tuple<Head*, Tail>& to, type& from, Return* ret /* unused */)
    {
        back_tuple<Tail, Return>::init(to.tail, from.tail, ret);
        to.head = &from.head;
    }
};

template<class T>
struct back_tuple<tuple<T*, empty_tuple>, T>
{
    typedef void type;
    static void init(tuple<T*, empty_tuple>& to, empty& from, T* ret)
    {
        to.head = ret;
    }
};
template<class Head, class Tail>
struct back_tuple<tuple<Head*, Tail>, Head>
{
    typedef typename back_tuple<Tail, void>::type type;
    static void init(tuple<Head*, Tail>& to, type& from, Head* ret)
    {
        to.head = ret;
        back_tuple<Tail, void>::init(to.tail, from, ret /* unused */ );
    }
};
template<class Return>
struct back_tuple<empty_tuple, Return>
{
    typedef empty_tuple type;
    static void init(empty_tuple& to, type& from, Return* ret)
    {
    }
};

template<class Tuple1, class Tuple2>
struct concat_tuple;

template<class Tuple2>
struct concat_tuple<empty_tuple, Tuple2>
{
    typedef Tuple2 type;
    static Tuple2& get_second(Tuple2& t) {return t;}
    static void init_first(empty_tuple&, type&) {};
};
template<class Head, class Tail, class Tuple2>
struct concat_tuple<tuple<Head, Tail>, Tuple2>
{
    typedef tuple<Head, typename concat_tuple<Tail, Tuple2>::type> type;
    static Tuple2& get_second(type& t)
    {
        return concat_tuple<Tail, Tuple2>::get_second(t.tail);
    }
    static void init_first(tuple<Head, Tail>& t, type& o)
    {
        t.head = o.head;
        concat_tuple<Tail, Tuple2>::init_first(t.tail, o.tail);
    }
};
}

namespace operations {
struct immediate {};
struct plus {};
} // operations

namespace rules {
template<class T>
struct no_op
{
    template<class U>
    static void doit(const U&) {}
};

struct UNIMPLEMENTED
{
    static const bool unimplemented_marker = true;
};

template<class T, class From, class Enable = void>
struct initialization
    : UNIMPLEMENTED
{ };

template<class T, class Enable = void>
struct print
    : UNIMPLEMENTED
{ };

template<class T, class Enable = void>
struct destruction
    : public no_op<T>
{};

template<class T, class Enable = void>
struct empty_initialization
    : UNIMPLEMENTED
{ };

template<class T, class U, class Enable = void>
struct assignment
    : UNIMPLEMENTED
{ };

// TODO priorities
// If result_is_temporary is true, then the result coincides with the
// first temporary (provided these have the same type)
template<class T, bool result_is_temporary, class Enable = void>
struct evaluation;
//{
//    typedef X return_t;
//    typedef Y temporaries_t; // a tuple of *pointers*
//    static void doit(const T& input, temporaries_t temps, return_t* output);
//};
} // rules

namespace mp {

template<class T>
struct not_
{
    static const bool val = !T::val;
};

template<class T, class U>
struct and_
{
    static const bool val = T::val && U::val;
};

template<bool, class U = void>
struct enable_if_v
{
    typedef U type;
    static const int val = 0;
};
template<class U>
struct enable_if_v<false, U>
{ };

template<bool v, class U = void>
struct disable_if_v : enable_if_v<!v, U> { };

template<class T, class U = void>
struct enable_if
    : public enable_if_v<T::val, U>
{ };

template<class T, class U = void>
struct disable_if
    : public disable_if_v<T::val, U>
{ };
}

namespace detail{
template<class T>
struct wrap
{
    T t;
};
}

namespace traits {
struct no { int data[2]; };
typedef int yes;
template<class T> detail::wrap<T> fakeinstance();

// use with care
template<class To, class From>
struct _is_convertible
{
private:
    static yes test(...) {return yes();}
    static no test(To) {return no();}
public:
    static const bool val = (sizeof(test(fakeinstance<From>().t)) != sizeof(yes));
};

template<class T>
struct is_implemented
    : public mp::not_<_is_convertible<rules::UNIMPLEMENTED, T> >
{ };

template<class T>
struct is_immediate
    : public _is_convertible<typename basetype<T>::type::operation_t, operations::immediate>
{ };
} // traits

namespace detail {
struct EXPRESSION
{ };

template<class Operation, class Expr>
struct evaluation_traits
{
    typedef typename Expr::derived_t derived_t;
    typedef rules::evaluation<derived_t, false> rule_t;
    typedef rules::evaluation<derived_t, true> temp_rule_t;
    typedef typename rule_t::return_t evaluation_return_t;
    typedef evaluation_return_t evaluated_t;

    static evaluation_return_t evaluate(const derived_t& from)
    {
        evaluation_return_t res;
        evaluate_into_fresh(res, from);
        return res;
    }

    static void evaluate_into(evaluation_return_t& to, const derived_t& from)
    {
        typedef mp::back_tuple<typename rule_t::temporaries_t> back_t;
        typename back_t::type temps_backing;
        typename rule_t::temporaries_t temps;
        back_t::init(temps, temps_backing, 0);
        rule_t::doit(from, temps, &to);
    }

    static void evaluate_into_fresh(evaluation_return_t& to, const derived_t& from)
    {
        typedef mp::back_tuple<typename temp_rule_t::temporaries_t, evaluation_return_t> back_t;
        typename back_t::type temps_backing;
        typename temp_rule_t::temporaries_t temps;
        back_t::init(temps, temps_backing, &to);
        temp_rule_t::doit(from, temps, &to);
    }
};

template<class Expr>
struct evaluation_traits<operations::immediate, Expr>
{
    typedef typename Expr::derived_t derived_t;
    typedef typename Expr::derived_t evaluated_t;
    typedef evaluated_t& evaluation_return_t;

    static evaluated_t& evaluate(derived_t& d) {return d;}
    static const evaluated_t& evaluate(const derived_t& d) {return d;}

    static void evaluate_into(derived_t& to, const derived_t& from)
    {
        rules::assignment<derived_t, derived_t>::doit(to, from);
    }

    static void evaluate_into_fresh(derived_t& to, const derived_t& from)
    {
        evaluate_into(to, from);
    }
};

template<class Operation, class Expr>
struct storage_traits
{
    typedef Expr type;
};

template<class Expr>
struct storage_traits<operations::immediate, Expr>
{
    typedef const Expr& type;
};
} // detail

namespace traits {
template<class T>
struct is_expression
    : public _is_convertible<detail::EXPRESSION, T>
{ };
} // traits

template<class Derived, class Operation, class Data>
class expression
    : public detail::EXPRESSION
{
protected:
    Data data;

    explicit expression(const Data & d) : data(d) {}

public:
    typedef detail::evaluation_traits<Operation, expression> ev_traits_t;
    typedef typename Derived::template type<Operation, Data>::result derived_t;
    typedef typename ev_traits_t::evaluated_t evaluated_t;
    typedef typename ev_traits_t::evaluation_return_t evaluation_return_t;
    typedef Data data_t;
    typedef Operation operation_t;

    // XXX private??
    derived_t& downcast() {return *static_cast<derived_t*>(this);}
    const derived_t& downcast() const {return *static_cast<const derived_t*>(this);}

public:
    // TODO strip qualifiers?
    template<class T>
    explicit expression(const T& t, typename mp::enable_if<traits::is_implemented<rules::initialization<derived_t, T> > >::type* = 0)
    {
        rules::initialization<derived_t, T>::doit(downcast(), t);
    }

    template<class T>
    explicit expression(const T& t,
            typename mp::disable_if<traits::is_implemented<rules::initialization<derived_t, T> > >::type* = 0,
            typename mp::enable_if<traits::is_expression<T> >::type* = 0)
    {
        T::ev_traits_t::evaluate_into_fresh(downcast(), t.downcast());
    }

    // XXX is this ok if non-implemented (i.e. may the compiler eagerly instantiate)?
    expression() {rules::empty_initialization<derived_t>::doit(downcast());}

    ~expression() {rules::destruction<derived_t>::doit(downcast());}

    Data& _data() {return data;}
    const Data& _data() const {return data;}

    // XXX move this stuff below the rules?
    void print(std::ostream& o) const
    {
        rules::print<evaluated_t>::doit(evaluate(), o);
    }

    typename traits::make_const<evaluation_return_t>::type evaluate() const {return ev_traits_t::evaluate(downcast());}
    evaluation_return_t evaluate() {return ev_traits_t::evaluate(downcast());}

    template<class T>
    void set(const T& t, typename mp::enable_if<traits::is_expression<T> >::type* = 0)
    {
        T::ev_traits_t::evaluate_into(downcast(), t.downcast());
    }
    template<class T>
    void set(const T& t, typename mp::disable_if<traits::is_expression<T> >::type* = 0)
    {
        rules::assignment<derived_t, T>::doit(downcast(), t);
    }

protected:
    template<class T, class Op>
    struct helper
    {
        typedef mp::make_tuple<
            typename detail::storage_traits<Operation, derived_t>::type,
            typename detail::storage_traits<typename T::operation_t, T>::type
          > maker;
        typedef typename maker::type type;
        typedef typename Derived::template type<Op, typename maker::type>::result new_derived_t;

        static new_derived_t
        make(const expression& self, const T& other)
        {
            return new_derived_t(maker::make(self.downcast(), other.downcast()));
        }
    };

public:
    template<class T>
    typename helper<T, operations::plus>::new_derived_t plus(const T& t) const
    {
        return helper<T, operations::plus>::make(*this, t);
    }

    template<class T>
    struct return_types
    {
        typedef typename helper<T, operations::plus>::new_derived_t plus_return_t;
    };
};

template<template<class O, class D> class Derived>
struct derived_wrapper
{
    template<class Operation, class Data>
    struct type
    {
        typedef Derived<Operation, Data> result;
    };
};


// operators

template<class Expr>
inline typename mp::enable_if<traits::is_expression<Expr>, std::ostream&>::type
operator<<(std::ostream& o, const Expr& e)
{
    e.print(o);
    return o;
}

// TODO other order
template<class Expr1, class Expr2>
inline typename mp::enable_if<traits::is_expression<Expr1>, typename Expr1::template return_types<Expr2>::plus_return_t>::type
operator+(const Expr1& e1, const Expr2& e2)
{
    return e1.plus(e2);
}


// "concrete" expression classes

template<class Operation, class Data>
class mpz_expression
    : public expression<derived_wrapper<mpz_expression>, Operation, Data>
{
public:
    mpz_expression() {}
    template<class T>
    explicit mpz_expression(const T& t) : mpz_expression::expression(t) {}

    template<class T>
    mpz_expression& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

//protected: // XXX
    explicit mpz_expression(const Data& d) : mpz_expression::expression(d) {}
};

} // flint

#include "fmpz.h"

namespace flint {
typedef mpz_expression<operations::immediate, fmpz_t> mpz;

namespace rules {
template<>
struct initialization<mpz, mpz>
{
    static void doit(mpz& target, const mpz& source)
    {
        fmpz_init_set(target._data(), source._data());
    }
};

template<class T>
struct initialization<mpz, T, typename mp::enable_if<traits::is_unsigned_integer<T> >::type>
{
    static void doit(mpz& target, T source)
    {
        fmpz_init_set_ui(target._data(), source);
    }
};

template<class T>
struct initialization<mpz, T, typename mp::enable_if<traits::is_signed_integer<T> >::type>
{
    static void doit(mpz& target, T source)
    {
        fmpz_init(target._data());
        fmpz_set_si(target._data(), source);
    }
};

template<>
struct assignment<mpz, mpz>
{
    static void doit(mpz& target, const mpz& source)
    {
        fmpz_set(target._data(), source._data());
    }
};

template<class T>
struct assignment<mpz, T, typename mp::enable_if<traits::is_unsigned_integer<T> >::type>
{
    static void doit(mpz& target, T source)
    {
        fmpz_set_ui(target._data(), source);
    }
};

template<class T>
struct assignment<mpz, T, typename mp::enable_if<traits::is_signed_integer<T> >::type>
{
    static void doit(mpz& target, T source)
    {
        fmpz_set_si(target._data(), source);
    }
};

template<int n>
struct assignment<mpz, char[n]>
{
    static void doit(mpz& target, const char* source)
    {
        fmpz_set_str(target._data(), const_cast<char*>(source), 10);
    }
};

template<>
struct empty_initialization<mpz>
{
    static void doit(mpz& v)
    {
        fmpz_init(v._data());
    }
};

template<>
struct destruction<mpz>
{
    static void doit(mpz& v)
    {
        fmpz_clear(v._data());
    }
};

template<>
struct print<mpz>
{
    static void doit(const mpz& v, std::ostream& o)
    {
        char* str = fmpz_get_str(0, 10, v._data());
        o << str;
        std::free(str);
    }
};

template<bool result_is_temporary>
struct evaluation<mpz_expression<operations::plus, tuple<const mpz&, tuple<const mpz&, empty_tuple> > >, result_is_temporary>
{
    typedef mpz return_t;
    typedef empty_tuple temporaries_t;
    static void doit(const mpz_expression<operations::plus, tuple<const mpz&, tuple<const mpz&, empty_tuple> > >& input, temporaries_t temps, return_t* output)
    {
        fmpz_add(output->_data(), input._data().first()._data(), input._data().second()._data());
    }
};

template<bool result_is_temporary, class Op, class Data1, class Data2>
struct evaluation<
    mpz_expression<Op, tuple<Data1, tuple<Data2, empty_tuple> > >,
    result_is_temporary,
    typename mp::enable_if<mp::and_<mp::not_<traits::is_immediate<Data1> >,
                                    mp::not_<traits::is_immediate<Data1> > > >::type>
{
    typedef mpz return_t;
    typedef evaluation<Data1, true> ev1_t;
    typedef evaluation<Data2, true> ev2_t;
    typedef mp::concat_tuple<tuple<mpz*, typename ev1_t::temporaries_t>, tuple<mpz*, typename ev2_t::temporaries_t> > concater;
    typedef typename concater::type temporaries_t;

    static void doit(const mpz_expression<Op, tuple<Data1, tuple<Data2, empty_tuple> > >& input, temporaries_t temps, return_t* output)
    {
        tuple<mpz*, typename ev1_t::temporaries_t> temps1;
        concater::init_first(temps1, temps);
        tuple<mpz*, typename ev2_t::temporaries_t> temps2 = concater::get_second(temps);
        ev1_t::doit(input._data().first(), temps1.tail, temps1.head);
        ev2_t::doit(input._data().second(), temps2.tail, temps2.head);
        *output = *temps1.head + *temps2.head;
    }
};
} // rules
} // flint

#endif
