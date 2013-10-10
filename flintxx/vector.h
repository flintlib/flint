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

// Sketch of a generic vector class.

#ifndef CXX_VECTOR_H
#define CXX_VECTOR_H

#include <string>
#include <sstream>

#include "expression.h"
#include "evaluation_tools.h"
#include "ltuple.h"
#include "mp.h"

namespace flint {
FLINT_DEFINE_BINOP(vector_at)

template<class Underlying_traits, class Operation, class Data>
class vector_expression;

namespace detail {
template<class Traits>
struct vector_wrapper : derived_wrapper2<vector_expression, Traits> { };

template<class Idx, class Operation, class Expr, class Traits>
struct vector_at_traits
{
    typedef FLINT_BINOP_ENABLE_RETTYPE(vector_at, Expr, Idx) ref_t;
    typedef ref_t cref_t;
    static ref_t at(const Expr& v, Idx i)
        {return vector_at(v, i);}
};
template<class Idx, class Expr, class Traits>
struct vector_at_traits<Idx, operations::immediate, Expr, Traits>
    : Traits { };
}

template<class Underlying_traits, class Operation, class Data>
class vector_expression
    : public expression<detail::vector_wrapper<Underlying_traits>, Operation, Data>
{
public:
    typedef expression<detail::vector_wrapper<Underlying_traits>,
                Operation, Data> base_t;
    typedef typename Underlying_traits::ref_t ref_t;
    typedef typename Underlying_traits::cref_t cref_t;
    typedef typename Underlying_traits::idx_t idx_t;
    typedef typename Underlying_traits::underlying_t underlying_t;
    typedef typename Underlying_traits::arrayref_t arrayref_t;
    typedef typename Underlying_traits::arraysrcref_t arraysrcref_t;

    vector_expression() {}

    template<class T>
    explicit vector_expression(const T& t) : base_t(t) {}
    template<class T, class U>
    vector_expression(const T& t, const U& u) : base_t(t, u) {}
    template<class T, class U, class V>
    vector_expression(const T& t, const U& u, const V& v)
        : base_t(t, u, v) {}

    template<class T>
    vector_expression& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

    template<class Idx>
    typename detail::vector_at_traits<Idx, Operation, vector_expression,
                 Underlying_traits>::ref_t operator[](Idx idx)
    {
        return detail::vector_at_traits<Idx, Operation, vector_expression,
                 Underlying_traits>::at(*this, idx);
    }
    template<class Idx>
    typename detail::vector_at_traits<Idx, Operation, vector_expression,
                 Underlying_traits>::cref_t operator[](Idx idx) const
    {
        return detail::vector_at_traits<Idx, Operation, vector_expression,
                 Underlying_traits>::at(*this, idx);
    }

    idx_t size() const {return Underlying_traits::size(*this);}

    arrayref_t _array() {return Underlying_traits::array(*this);}
    arraysrcref_t _array() const {return Underlying_traits::array(*this);}

    typename base_t::evaluated_t create_temporary() const
    {
        return Underlying_traits::create_temporary(*this);
    }

protected:
    explicit vector_expression(const Data& d) : base_t(d) {}

    template<class D, class O, class Da>
    friend class expression;
};

namespace vectors {
// Similar to matrices, the size of a vector expression has to be known in
// order to allocate temporary objects. In this case, the generic
// implementation looks for any vector immediate subexpression and returs its
// size. This makes sense since mixing vectors of differing sizes usually makes
// no sense.
// Thus specialisation is usually only necessary in constructor-like operations,
// which do not involve vector immediates.
template<class Operation>
struct outsize
{
    template<class Expr>
    static unsigned get(const Expr& e)
    {
        return tools::find_subexpr_T<typename Expr::evaluated_t>(e)._data().size;
    }
};

// Hack for ltuple_get, similar to the matrices case.
template<unsigned n>
struct outsize<operations::ltuple_get_op<n> >
{
    template<class Expr>
    static unsigned get(const Expr& e)
    {
        return outsize<typename Expr::data_t::head_t::operation_t>::get(
                e._data().head);
    }
};
}

namespace detail {
template<class T, class Ref, class Cref, class ArrayT>
struct basic_vector_traits
{
    typedef unsigned idx_t;
    typedef Ref ref_t;
    typedef const Cref cref_t;
    typedef T underlying_t;
    typedef ArrayT* arrayref_t;
    typedef const ArrayT* arraysrcref_t;

    template<class Expr>
    static ref_t at(Expr& e, unsigned i)
    {
        return e.evaluate()._data().array[i];
    }

    template<class Expr>
    static cref_t at(const Expr& e, unsigned i)
    {
        return e.evaluate()._data().array[i];
    }

    template<class Expr>
    static arrayref_t array(Expr& e)
    {
        return e.evaluate()._data().array;
    }

    template<class Expr>
    static arraysrcref_t array(const Expr& e)
    {
        return e.evaluate()._data().array;
    }
};
template<class T, class Ref = T&, class Cref = const T&, class ArrayT = T>
struct rtfixed_size_traits
    : basic_vector_traits<T, Ref, Cref, ArrayT>
{
    template<class Expr>
    static unsigned size(const Expr& e)
    {
        return vectors::outsize<typename Expr::operation_t>::get(e);
    }

    template<class Expr>
    static typename Expr::evaluated_t create_temporary(const Expr& e)
    {
        return typename Expr::evaluated_t(e.size());
    }
};
template<class T, class Ref = T&, class Cref = const T&, class ArrayT = T>
struct fixed_size_traits
    : basic_vector_traits<T, Ref, Cref, ArrayT>
{
    template<class Expr>
    static unsigned size(const Expr& e)
    {
        return Expr::evaluated_t::data_t::size;
    }

    template<class Expr>
    static typename Expr::evaluated_t create_temporary(const Expr& e)
    {
        return typename Expr::evaluated_t();
    }
};

template<class T, class Size, class Ref, class Cref, class ArrayT>
struct wrapped_vector_traits
    : rtfixed_size_traits<T, Ref, Cref, ArrayT>
{
    typedef Size idx_t;

    template<class Expr>
    static Ref at(Expr& e, idx_t i)
    {
        return e.evaluate()._data().at(i);
    }

    template<class Expr>
    static Cref at(const Expr& e, idx_t i)
    {
        return e.evaluate()._data().at(i);
    }
};

template<class T>
struct rtfixed_size_data
{
    const unsigned size;
    T* array;

    rtfixed_size_data(unsigned n)
        : size(n), array(new T[n]) {}
    ~rtfixed_size_data() {delete[] array;}

    rtfixed_size_data(const rtfixed_size_data& o)
        : size(o.size)
    {
        // TODO this is very non-optimal ... (?)
        array = new T[size];
        for(unsigned i = 0;i < size;++i)
        {
            array[i] = o.array[i];
        }
    }
};
template<class T, unsigned n>
struct fixed_size_data
{
    static const unsigned size = n;
    T array[n];
};
} // detail

template<class T>
struct make_vector
{
    typedef vector_expression<detail::rtfixed_size_traits<T>,
                operations::immediate, detail::rtfixed_size_data<T> > type;
};
template<class T, unsigned n>
struct make_vector_n
{
    typedef vector_expression<detail::fixed_size_traits<T>,
                operations::immediate, detail::fixed_size_data<T, n> > type;
};

template<class Expr>
struct enable_vector_rules : mp::false_ { };

template<class Traits, class Data>
struct enable_vector_rules<
    vector_expression<Traits, operations::immediate, Data> >
    : mp::true_ { };

namespace rules {
// temporary allocation inside ltuples
template<class Operation, class Data, class U,
    class Traits, class Op, class Da>
struct instantiate_temporaries<ltuple_expression<U, Operation, Data>,
    vector_expression<Traits, Op, Da> >
{
    typedef ltuple_expression<U, Operation, Data> Expr;
    typedef vector_expression<Traits, Op, Da> T;
    static T get(const Expr& e)
    {
        return T(vectors::outsize<Operation>::get(e));
    }
};

template<class Traits, class Data, class T>
struct binary_expression<vector_expression<Traits, operations::immediate, Data>,
    operations::vector_at_op, T>
{
    typedef typename Traits::underlying_t return_t;
    template<class V>
    static void doit(V& to,
            const vector_expression<Traits, operations::immediate, Data>& v,
            T i)
    {
        to = Traits::at(v, i);
    }
};


template<class Expr>
struct to_string<Expr, typename mp::enable_if<mp::and_<
    enable_vector_rules<Expr>,
    traits::is_implemented<to_string<typename Expr::underlying_t> > > >::type>
{
    static std::string get(const Expr& e, int base)
    {
        // TODO inefficient
        std::string res = "(";
        for(typename Expr::idx_t i = 0;i < e.size();++i)
        {
            res += e[i].to_string();
            if(i != e.size() - 1)
                res += ", ";
        }
        res += ")";
        return res;
    }
};

template<class Expr>
struct equals<Expr, Expr,
    typename mp::enable_if<enable_vector_rules<Expr> >::type>
{
    static bool get(const Expr& e1, const Expr& e2)
    {
        if(e1.size() != e2.size())
            return false;
        for(typename Expr::idx_t i = 0;i < e1.size();++i)
            if(e1[i] != e2[i])
                return false;
        return true;
    }
};

namespace rvdetail {
template<class Tuple>
struct translate_data;

template<class Expr, class enable = void>
struct translate_expr
{
    typedef translate_data<typename Expr::data_t> trdata_t;
    typedef typename Expr::underlying_t ul_t;
    typedef typename ul_t::template make_helper<
        typename Expr::operation_t, typename trdata_t::type> make_helper;
    typedef typename make_helper::type type;

    template<class Idx>
    static type make(const Expr& e, Idx idx)
    {
        return make_helper::make(trdata_t::make(e._data(), idx));
    }
};

template<class Expr>
struct translate_expr<Expr,
    typename mp::enable_if<traits::is_immediate<Expr> >::type>
{
    typedef typename Expr::cref_t type;

    template<class Idx>
    static type make(const Expr& e, Idx idx)
    {
        return e[idx];
    }
};

template<class Head, class Tail>
struct translate_data<tuple<Head, Tail> >
{
    typedef translate_expr<typename traits::basetype<Head>::type> trexpr;
    typedef translate_data<Tail> trtail;
    typedef tuple<typename trexpr::type, typename trtail::type> type;

    template<class Idx>
    static type make(const tuple<Head, Tail>& e, Idx idx)
    {
        return type(trexpr::make(e.head, idx), trtail::make(e.tail, idx));
    }
};
template<>
struct translate_data<empty_tuple>
{
    typedef empty_tuple type;
    template<class Idx>
    static type make(empty_tuple, Idx) {return empty_tuple();}
};

template<class Data, class Enable = void>
struct enable_evaluation : mp::false_ {typedef void vector_t;};

template<class Data>
struct enable_evaluation<Data,
    typename mp::enable_if<traits::is_expression<
        typename traits::basetype<Data>::type> >::type>
    : enable_vector_rules<typename traits::basetype<Data>::type::evaluated_t>
{
    typedef typename traits::basetype<Data>::type::evaluated_t vector_t;
};
template<class Head, class Tail>
struct enable_evaluation<tuple<Head, Tail> >
    : mp::and_<enable_evaluation<Head>, enable_evaluation<Tail> >
{
    typedef typename enable_evaluation<Head>::vector_t vector_t;
};
template<>
struct enable_evaluation<empty_tuple>
    : mp::true_ { };
} //rvdetail

// TODO this is a bit greedy ..
template<class Op, class Data, bool result_is_temporary>
struct evaluation<Op, Data, result_is_temporary, 1,
    typename mp::enable_if<rvdetail::enable_evaluation<Data> >::type>
{
    typedef rvdetail::translate_data<Data> translator;
    typedef typename translator::type trdata_t;
    typedef typename mp::find_evaluation<
        Op, trdata_t, result_is_temporary>::type rule_t;
    typedef typename rvdetail::enable_evaluation<Data>::vector_t vector_t;
    typedef typename vector_t::evaluated_t return_t; // TODO
    typedef typename rule_t::temporaries_t temporaries_t;
    typedef typename rule_t::return_t trreturn_t;

    template<class Return>
    static void doit(const Data& input, temporaries_t temps, Return* output)
    {
        for(typename return_t::idx_t i = 0;i < output->size();++i)
        {
            rule_t::doit(translator::make(input, i), temps, &((*output)[i]));
        }
    }
};

// TODO scalar multiplication etc
} // rules
} // flint
#endif
