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

#ifndef CXX_EXPRESSION_H
#define CXX_EXPRESSION_H

// TODO
// * static asserts

#include <iosfwd>

#include "cxx/mp.h"
#include "cxx/traits.h"
#include "cxx/tuple.h"

namespace flint {

namespace operations {
// TODO document these
struct immediate { };
struct plus { };
} // operations

namespace rules {
// TODO document these
struct no_op
{
    template<class U>
    static void doit(const U&) {}
};

struct UNIMPLEMENTED
{
    static const bool unimplemented_marker = true;
};

// NB: this is mostly internal
template<class T, class Enable = void>
struct copy_initialization : UNIMPLEMENTED { };

template<class T, class From, class Enable = void>
struct initialization : UNIMPLEMENTED { };

template<class T, class Enable = void>
struct print : UNIMPLEMENTED { };

template<class T, class Enable = void>
struct destruction : no_op { };

template<class T, class Enable = void>
struct empty_initialization : UNIMPLEMENTED { };

template<class T, class U, class Enable = void>
struct assignment : UNIMPLEMENTED { };

template<class T, class U, class Enable = void>
struct equals : UNIMPLEMENTED { };

// If result_is_temporary is true, then the result coincides with the
// first temporary (provided these have the same type)
// Priorities 2, 1, 0 can be used to resolve conflicts.
template<
    class T, class Op, class Data,
    bool result_is_temporary,
    unsigned priority,
    class Enable = void>
struct evaluation : UNIMPLEMENTED { };
//{
//    typedef X return_t;
//    typedef Y temporaries_t; // a tuple of *pointers*
//    static void doit(const T& input, temporaries_t temps, return_t* output);
//};

// Convenience helpers, instantiate by evaluation if necessary
template<class T, class Op, class U>
struct binary_expression : UNIMPLEMENTED { };
// typedef X return_t;
// static void doit(return_t& to, const T&, const U&);
template<class T, class Op, class U>
struct commutative_binary_expression : UNIMPLEMENTED { };
// similarly
} // rules

namespace traits {
template<class T>
struct is_implemented : mp::not_<_is_convertible<rules::UNIMPLEMENTED, T> > { };
} // traits

namespace mp {
template<class T, class Op, class Data,
    bool result_is_temporary, unsigned min_prio = 0>
struct find_evaluation
{
private:
    typedef rules::evaluation<T, Op, Data, result_is_temporary, 2> r2;
    typedef rules::evaluation<T, Op, Data, result_is_temporary, 1> r1;
    typedef rules::evaluation<T, Op, Data, result_is_temporary, 0> r0;

    typedef traits::is_implemented<r2> i2;
    typedef traits::is_implemented<r1> i1;
    typedef traits::is_implemented<r0> i0;

public:
    typedef typename mp::select<rules::UNIMPLEMENTED, // TODO
        mp::and_c<i2, min_prio <= 2>, r2,
        mp::and_c<i1, min_prio <= 1>, r1,
        mp::and_c<i0, min_prio <= 0>, r0
      >::type type;
};
} // mp

namespace detail {
struct EXPRESSION { };
}

namespace traits {
template<class T>
struct is_expression : _is_convertible<detail::EXPRESSION, T> { };

template<class T>
struct _is_immediate_expr
    : _is_convertible<
          typename basetype<T>::type::operation_t,
          operations::immediate
        >
{ };

// Compute if T is an expression, with operation "immediate"
template<class T, class Enable = void>
struct is_immediate_expr : _is_immediate_expr<T> { };
template<class T>
struct is_immediate_expr<T,
  typename mp::enable_if<mp::not_<is_expression<T> > >::type>
    : mp::false_ { };

// Compute if T is an immediate expression, *or not an expression at all*
template<class T>
struct is_immediate
    : mp::or_<mp::not_<is_expression<T> >, is_immediate_expr<T> > { };

// Compute if T is a non-immediate expression
template<class T>
struct is_lazy_expr
    : mp::and_<is_expression<T>, mp::not_<is_immediate_expr<T> > > { };
} // traits

namespace detail {
template<class Operation, class Expr, class Data>
struct evaluation_traits
{
    typedef typename Expr::derived_t derived_t;
    typedef typename mp::find_evaluation<derived_t, Operation, Data, false>::type rule_t;
    typedef typename mp::find_evaluation<derived_t, Operation, Data, true>::type temp_rule_t;
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
        typedef mp::back_tuple<
            typename temp_rule_t::temporaries_t,
            evaluation_return_t
          > back_t;
        typename back_t::type temps_backing;
        typename temp_rule_t::temporaries_t temps;
        back_t::init(temps, temps_backing, &to);
        temp_rule_t::doit(from, temps, &to);
    }
};

template<class Expr, class Data>
struct evaluation_traits<operations::immediate, Expr, Data>
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
} // detail

template<class Derived, class Operation, class Data>
class expression
    : public detail::EXPRESSION
{
private:
    Data data;

protected:
    explicit expression(const Data & d) : data(d) {}

public:
    typedef detail::evaluation_traits<Operation, expression, Data> ev_traits_t;
    typedef typename Derived::template type<Operation, Data>::result derived_t;
    typedef typename ev_traits_t::evaluated_t evaluated_t;
    typedef typename ev_traits_t::evaluation_return_t evaluation_return_t;
    typedef Data data_t;
    typedef Operation operation_t;

private:
    derived_t& downcast() {return *static_cast<derived_t*>(this);}
    const derived_t& downcast() const
    {
        return *static_cast<const derived_t*>(this);
    }

public:
    // TODO strip qualifiers?
    template<class T>
    explicit expression(const T& t,
            typename mp::enable_if<traits::is_implemented<
                rules::initialization<derived_t, T> > >::type* = 0)
    {
        rules::initialization<derived_t, T>::doit(downcast(), t);
    }

    template<class T>
    explicit expression(const T& t,
            typename mp::disable_if<traits::is_implemented<
                rules::initialization<derived_t, T> > >::type* = 0,
            typename mp::enable_if<traits::is_expression<T> >::type* = 0)
    {
        T::ev_traits_t::evaluate_into_fresh(downcast(), t);
    }

    expression(const expression& e)
        : data(rules::copy_initialization<derived_t>::get(e.downcast()))
    {
    }

    // NB: the compiler is not allowed to eagerly instantiate!
    expression() {rules::empty_initialization<derived_t>::doit(downcast());}

    ~expression() {rules::destruction<derived_t>::doit(downcast());}

    Data& _data() {return data;}
    const Data& _data() const {return data;}

    void print(std::ostream& o) const
    {
        rules::print<evaluated_t>::doit(evaluate(), o);
    }

    typename traits::make_const<evaluation_return_t>::type evaluate() const
    {
        return ev_traits_t::evaluate(downcast());
    }
    evaluation_return_t evaluate() {return ev_traits_t::evaluate(downcast());}

    template<class T>
    void set(const T& t,
            typename mp::enable_if<traits::is_expression<T> >::type* = 0)
    {
        T::ev_traits_t::evaluate_into(downcast(), t);
    }
    template<class T>
    void set(const T& t,
            typename mp::disable_if<traits::is_expression<T> >::type* = 0)
    {
        rules::assignment<derived_t, T>::doit(downcast(), t);
    }

    template<class T>
    bool equals(const T& t,
            typename mp::enable_if<traits::is_lazy_expr<T> >::type* = 0) const
    {
        return equals(t.evaluate());
    }
    template<class T>
    bool equals(const T& t,
            typename mp::disable_if<traits::is_lazy_expr<T> >::type* = 0) const
    {
        return rules::equals<evaluated_t, T>::get(evaluate(), t);
    }

    template<class Op, class NData>
    struct make_helper
    {
        typedef typename Derived::template type<Op, NData>::result type;
        static type make(const NData& ndata)
        {
            return type(ndata);
        }
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

namespace detail {
template<class Expr>
struct storage_traits
    : mp::if_<
          traits::is_immediate<Expr>,
          typename traits::forwarding<Expr>::type,
          Expr
        > { };

template<class Expr1, class Op, class Expr2>
struct binary_op_helper
{
    typedef mp::make_tuple<
        typename storage_traits<Expr1>::type,
        typename storage_traits<Expr2>::type
      > maker;
    typedef typename maker::type type;
    typedef typename mp::if_<traits::is_expression<Expr1>, Expr1, Expr2>::type Expr;
    typedef typename Expr::template make_helper<Op, type> make_helper;
    typedef typename make_helper::type return_t;

    static return_t make(const Expr1& left, const Expr2& right)
    {
        return make_helper::make(maker::make(left, right));
    }
};
}

template<class Expr>
inline typename mp::enable_if<traits::is_expression<Expr>, std::ostream&>::type
operator<<(std::ostream& o, const Expr& e)
{
    e.print(o);
    return o;
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<traits::is_expression<Expr1>, bool>::type
operator==(const Expr1& e1, const Expr2& e2)
{
  return e1.equals(e2);
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<mp::and_<
        mp::not_<traits::is_expression<Expr1> >,
        traits::is_expression<Expr2> >,
    bool>::type
operator==(const Expr1& e1, const Expr2& e2)
{
  return e2.equals(e1);
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<mp::or_<
        traits::is_expression<Expr1>,
        traits::is_expression<Expr2> >,
    bool>::type
operator!=(const Expr1& e1, const Expr2& e2)
{
  return !(e1 == e2);
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<mp::or_<
        traits::is_expression<Expr1>,
        traits::is_expression<Expr2> >,
    typename detail::binary_op_helper<Expr1, operations::plus, Expr2>::return_t>::type
operator+(const Expr1& e1, const Expr2& e2)
{
    return detail::binary_op_helper<Expr1, operations::plus, Expr2>::make(e1, e2);
}

// default rules

namespace rules {
// Copy initialization is special in that data does not necessarily have
// a default constructor.
// TODO could store pointers instead of references.
template<class T>
struct copy_initialization<T,
    typename mp::disable_if<traits::is_immediate<T> >::type>
{
    static typename T::data_t get(const T& from)
    {
        return typename T::data_t(from._data());
    }
};

template<class T>
struct copy_initialization<T,
    typename mp::enable_if<traits::is_immediate<T> >::type>
{
    static typename T::data_t get(const T& from)
    {
        T tmp;
        initialization<T, T>::doit(tmp, from);
        return tmp._data();
    }
};

// Composite binary operators
// These rules implement binary operators by implementing both arguments
// separately, then performing the operation on the evaluated types by
// instantiating the appropriate rule again.
//
// Hence to evaluate expressions like a + (b + c), it suffices to write
// rules for composition of two immediates.

// TODO make these more efficient
template<class T, bool result_is_temporary, class Op, class Data1, class Data2>
struct evaluation<
    T, Op, tuple<Data1, tuple<Data2, empty_tuple> >, result_is_temporary, 0,
    typename mp::enable_if<
        mp::and_<
            traits::is_lazy_expr<Data1>,
            traits::is_lazy_expr<Data2> > >::type
  >
{
    typedef typename mp::find_evaluation<
        Data1,
        typename Data1::operation_t,
        typename Data1::data_t,
        true
      >::type ev1_t;
    typedef typename mp::find_evaluation<
        Data2,
        typename Data2::operation_t,
        typename Data2::data_t,
        true
      >::type ev2_t;
    typedef typename ev1_t::return_t return1_t;
    typedef typename ev1_t::temporaries_t temporaries1_t;
    typedef typename ev2_t::return_t return2_t;
    typedef typename ev2_t::temporaries_t temporaries2_t;

    typedef detail::binary_op_helper<return1_t, Op, return2_t> binop_helper;
    typedef typename binop_helper::return_t::evaluated_t return_t;

    typedef mp::concat_tuple<
        tuple<return1_t*, typename ev1_t::temporaries_t>,
        tuple<return2_t*, typename ev2_t::temporaries_t> > concater;
    typedef typename concater::type temporaries_t;

    static void doit(const T& input, temporaries_t temps, return_t* output)
    {
        tuple<return1_t*, typename ev1_t::temporaries_t> temps1 =
            concater::get_first(temps);
        tuple<return2_t*, typename ev2_t::temporaries_t> temps2 =
            concater::get_second(temps);
        ev1_t::doit(input._data().first(), temps1.tail, temps1.head);
        ev2_t::doit(input._data().second(), temps2.tail, temps2.head);
        *output = binop_helper::make(*temps1.head, *temps2.head);
    }
};

template<class T, bool result_is_temporary, class Op, class Data1, class Data2>
struct evaluation<
    T, Op, tuple<Data1, tuple<Data2, empty_tuple> >, result_is_temporary, 0,
    typename mp::enable_if<
        mp::and_<
            traits::is_immediate<typename traits::basetype<Data1>::type>,
            traits::is_lazy_expr<Data2> > >::type
  >
{
    typedef typename mp::find_evaluation<
        Data2,
        typename Data2::operation_t,
        typename Data2::data_t,
        true
      >::type ev2_t;
    typedef typename ev2_t::return_t return2_t;
    typedef typename ev2_t::temporaries_t temporaries2_t;
    typedef typename traits::basetype<Data1>::type return1_t;

    typedef detail::binary_op_helper<return1_t, Op, return2_t> binop_helper;
    typedef typename binop_helper::return_t::evaluated_t return_t;

    typedef tuple<return2_t*, temporaries2_t> temporaries_t;

    static void doit(const T& input, temporaries_t temps, return_t* output)
    {
        ev2_t::doit(input._data().second(), temps.tail, temps.head);
        *output = binop_helper::make(input._data().first(), *temps.head);
    }
};

template<class T, bool result_is_temporary, class Op, class Data1, class Data2>
struct evaluation<
    T, Op, tuple<Data1, tuple<Data2, empty_tuple> >, result_is_temporary, 0,
    typename mp::enable_if<
        mp::and_<
            traits::is_immediate<typename traits::basetype<Data2>::type>,
            traits::is_lazy_expr<Data1> > >::type
  >
{
    typedef typename mp::find_evaluation<
        Data1,
        typename Data1::operation_t,
        typename Data1::data_t,
        true
      >::type ev1_t;
    typedef typename ev1_t::return_t return1_t;
    typedef typename ev1_t::temporaries_t temporaries1_t;
    typedef typename traits::basetype<Data2>::type return2_t;

    typedef detail::binary_op_helper<return1_t, Op, return2_t> binop_helper;
    typedef typename binop_helper::return_t::evaluated_t return_t;

    typedef tuple<return1_t*, temporaries1_t> temporaries_t;

    static void doit(const T& input, temporaries_t temps, return_t* output)
    {
        ev1_t::doit(input._data().first(), temps.tail, temps.head);
        *output = binop_helper::make(*temps.head, input._data().second());
    }
};

// Automatically invoke binary_expression or commutative_binary_expression
namespace rdetail {
template<class Expr1, class Op, class Expr2>
struct inverted_binary_expression
{
  typedef commutative_binary_expression<Expr2, Op, Expr1> wrapped_t;
  typedef typename wrapped_t::return_t return_t;
  static void doit(return_t& to, const Expr1& e1, const Expr2& e2)
  {
    return wrapped_t::doit(to, e2, e1);
  }
};

template<template<class E1, class O, class E2> class BE, class T, class Data1, class Op, class Data2>
struct binary_expr_helper
{
    typedef typename traits::basetype<Data1>::type data1_t;
    typedef typename traits::basetype<Data2>::type data2_t;
    typedef BE<data1_t, Op, data2_t> wrapped_t;
    typedef typename wrapped_t::return_t return_t;
    typedef empty_tuple temporaries_t;
    static void doit(const T& input, temporaries_t temps, return_t* output)
    {
        wrapped_t::doit(*output, input._data().first(), input._data().second());
    }
};
}

template<class T, bool result_is_temporary, class Op, class Data1, class Data2>
struct evaluation<
    T, Op, tuple<Data1, tuple<Data2, empty_tuple> >, result_is_temporary, 0,
    typename mp::enable_if<
        mp::and_<
            traits::is_immediate<typename traits::basetype<Data1>::type>,
            mp::and_<
                traits::is_immediate<typename traits::basetype<Data2>::type>,
                mp::or_<
                    traits::is_implemented<binary_expression<
                        typename traits::basetype<Data1>::type,
                        Op,
                        typename traits::basetype<Data2>::type
                      > >,
                    mp::or_<
                        traits::is_implemented<commutative_binary_expression<
                            typename traits::basetype<Data1>::type,
                            Op,
                            typename traits::basetype<Data2>::type
                          > >,
                        traits::is_implemented<commutative_binary_expression<
                            typename traits::basetype<Data2>::type,
                            Op,
                            typename traits::basetype<Data1>::type
                          > >
                      >
                  >
              >
          >
      >::type>
    : mp::if_<
            traits::is_implemented<binary_expression<
                typename traits::basetype<Data1>::type,
                Op,
                typename traits::basetype<Data2>::type
              > >,
            rdetail::binary_expr_helper<binary_expression, T, Data1, Op, Data2>,
            typename mp::if_<
                traits::is_implemented<commutative_binary_expression<
                    typename traits::basetype<Data1>::type,
                    Op,
                    typename traits::basetype<Data2>::type
                  > >,
                rdetail::binary_expr_helper<
                    commutative_binary_expression, T, Data1, Op, Data2>,
                rdetail::binary_expr_helper<
                    rdetail::inverted_binary_expression, T, Data1, Op, Data2>
              >::type
          >::type
{ };
} // rules
} // flint

#endif
