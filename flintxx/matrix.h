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

// Common code shared among matrix classes

#ifndef FLINTXX_MATRIX_H
#define FLINTXX_MATRIX_H

#include "flint_classes.h"
#include "mp.h"
#include "rules.h"
#include "stdmath.h"
#include "ltuple.h"
#include "traits.h"
#include "tuple.h"
#include "../permxx.h"

namespace flint {
FLINT_DEFINE_BINOP(solve)
FLINT_DEFINE_BINOP(solve_fflu)
FLINT_DEFINE_THREEARY(mat_at)
FLINT_DEFINE_THREEARY(solve_fflu_precomp)
FLINT_DEFINE_UNOP(charpoly)
FLINT_DEFINE_UNOP(det)
FLINT_DEFINE_UNOP(det_fflu)
FLINT_DEFINE_UNOP(det_interpolate)
FLINT_DEFINE_UNOP(nullspace)
FLINT_DEFINE_UNOP(rref)
FLINT_DEFINE_UNOP(trace)
FLINT_DEFINE_UNOP(transpose)

FLINT_DEFINE_THREEARY(fflu)
FLINT_DEFINE_THREEARY_HERE_2DEFAULT(fflu, permxx*, 0, bool, false)

template<class T> struct matrix_traits : rules::UNIMPLEMENTED { };
// override this for the non-reference type of your choice
// {
//     template<class M> static slong rows(const M&);
//     template<class M> static slong cols(const M&);
//     template<class M> static ??? at(const M&, slong, slong);
//     template<class M> static ??? at(M&, slong, slong);
// };
namespace traits {
template <class T, class Enable = void> struct is_mat : is_implemented<
    matrix_traits<typename flint_classes::to_nonref<T>::type> > { };
template<class T> struct is_mat<T,
    typename mp::disable_if<flint_classes::is_flint_class<T> >::type>
        : mp::false_ { };
} // traits
namespace matrices {
namespace mdetail {
// Some helper traits used below.
template<class Data> struct second_is_mat_data : mp::false_ { };
template<class Data1, class Data2>
struct second_is_mat_data<tuple<Data1, tuple<Data2, empty_tuple> > >
    : traits::is_mat<typename traits::basetype<Data2>::type> { };
template<class Expr> struct second_is_mat
    : second_is_mat_data<typename Expr::data_t> { };

template<class Mat>
struct both_mat : mp::and_<
    traits::is_mat<
        typename traits::basetype<typename Mat::data_t::head_t>::type>,
    traits::is_mat<
        typename traits::basetype<typename Mat::data_t::tail_t::head_t>::type>
  > { };

// A more convenient way to obtain the traits associated to a non-immediate
// or non-nonref expression.
template<class Expr> struct immediate_traits
    : matrix_traits<typename flint_classes::to_nonref<Expr>::type> { };
} // mdetail

// For matrix expressions to create temporaries, it is necessary to know the
// dimensions of the result of a computation. This is a generic implementation,
// which assumes that the output dimensions are the same as the dimensions of
// the first argument, which is assumed to be a matrix, except if there are
// precisely two arguments only the second of which is a matrix, in which case
// we assume its the dimension of that.
// This implementation works correctly in many cases, e.g. matrix-addition or
// matrix-scalar multiplication.
template<class Operation>
struct outsize_generic
{
    template<class Mat>
    static slong rows(const Mat& m,
            typename mp::disable_if<mdetail::second_is_mat<Mat> >::type* = 0)
    {
        return m._data().head.rows();
    }
    template<class Mat>
    static slong cols(const Mat& m,
            typename mp::disable_if<mdetail::second_is_mat<Mat> >::type* = 0)
    {
        return m._data().head.cols();
    }

    template<class Mat>
    static slong rows(const Mat& m,
            typename mp::enable_if<mdetail::second_is_mat<Mat> >::type* = 0)
    {
        return m._data().tail.head.rows();
    }
    template<class Mat>
    static slong cols(const Mat& m,
            typename mp::enable_if<mdetail::second_is_mat<Mat> >::type* = 0)
    {
        return m._data().tail.head.cols();
    }
};

// This is the expression template used for computing the dimensions of an
// operation. Without further specialisation, it is just the generic
// implementation described above.
// If you introduce a new operation where the generic implementation is
// incorrect, you must specialise this template.
template<class Operation>
struct outsize : outsize_generic<Operation> { };

// Specialise immediates, where the dimensions are stored with the object.
template<>
struct outsize<operations::immediate>
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return mdetail::immediate_traits<Mat>::rows(m);
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return mdetail::immediate_traits<Mat>::cols(m);
    }
};

// Specialise multiplication. For matrix-matrix multiplication, use
// the usual formula. For matrix-scalar multiplication, use the generic
// implementation.
template<>
struct outsize<operations::times>
{
    template<class Mat>
    static slong rows(const Mat& m,
            typename mp::enable_if<mdetail::both_mat<Mat> >::type* = 0)
    {
        return m._data().head.rows();
    }
    template<class Mat>
    static slong cols(const Mat& m,
            typename mp::enable_if<mdetail::both_mat<Mat> >::type* = 0)
    {
        return m._data().tail.head.cols();
    }

    template<class Mat>
    static slong rows(const Mat& m,
            typename mp::disable_if<mdetail::both_mat<Mat> >::type* = 0)
    {
        return outsize_generic<operations::times>::rows(m);
    }
    template<class Mat>
    static slong cols(const Mat& m,
            typename mp::disable_if<mdetail::both_mat<Mat> >::type* = 0)
    {
        return outsize_generic<operations::times>::cols(m);
    }
};
// Any particular multipication algorithm also has to be specialised.
template<>
struct outsize<operations::mul_classical_op>
    : outsize<operations::times> { };

// Specialise transpose.
template<>
struct outsize<operations::transpose_op>
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return m._data().head.cols();
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return m._data().head.rows();
    }
};

// Specialise nullspace. Note that the nullspace computation functions in
// flint return a matrix the columns of which span the nullspace. Since the
// nullity is not known in advance in general, we have to allocate a square
// matrix.
template<>
struct outsize<operations::nullspace_op>
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return m._data().head.cols();
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return m._data().head.cols();
    }
};

// This is a bit of a hack. Matrix operations returning a tuple typically
// only return one matrix. We key outsize on the inner operation to find out
// the dimensions. So e.g. solve(A, X).get<1>() (say) will invoke outsize
// with ltuple_get_op<1> as argument, which then invokes outsize with solve_op
// and (A, X) as argument.
template<unsigned n>
struct outsize<operations::ltuple_get_op<n> >
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return outsize<
            typename Mat::data_t::head_t::operation_t>::rows(m._data().head);
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return outsize<
            typename Mat::data_t::head_t::operation_t>::cols(m._data().head);
    }
};

// This is not actually a matrix expression, but called by the above ...
template<>
struct outsize<operations::solve_op>
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return m._data().second().rows();
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return m._data().second().cols();
    }
};
template<> struct outsize<operations::solve_fflu_op>
    : outsize<operations::solve_op> { };
template<>
struct outsize<operations::solve_fflu_precomp_op>
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return m._data().tail.second().rows();
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return m._data().tail.second().cols();
    }
};

namespace mdetail {
struct base_traits
{
    template<class M>
    static slong rows(const M& m)
    {
        return matrices::outsize<typename M::operation_t>::rows(m);
    }
    template<class M>
    static slong cols(const M& m)
    {
        return matrices::outsize<typename M::operation_t>::cols(m);
    }
};
} // mdetail

// These traits classes are useful for implementing unified coefficient access.
// See fmpz_matxx etc for example usage.
template<class Mat>
struct generic_traits : mdetail::base_traits
{
    template<class T, class U>
    struct at
    {
        typedef FLINT_THREEARY_ENABLE_RETTYPE(mat_at, Mat, T, U)
            entry_ref_t;
        typedef entry_ref_t entry_srcref_t;

        static entry_srcref_t get(const Mat& m, T i, U j)
        {
            return mat_at(m, i, j);
        }
    };
};

template<class Srcref>
struct generic_traits_srcref : mdetail::base_traits
{
    template<class T, class U>
    struct at
    {
        typedef Srcref entry_ref_t;
        typedef Srcref entry_srcref_t;

        template<class M>
        static Srcref get(M m, T i, U j)
        {
            return mdetail::immediate_traits<M>::at(m, i, j);
        }
    };
};

template<class Ref>
struct generic_traits_ref : mdetail::base_traits
{
    template<class T, class U>
    struct at
    {
        typedef Ref entry_ref_t;
        typedef Ref entry_srcref_t;

        template<class M>
        static Ref get(M m, T i, U j)
        {
            return mdetail::immediate_traits<M>::at(m, i, j);
        }
    };
};

template<class Ref, class Srcref>
struct generic_traits_nonref : mdetail::base_traits
{
    template<class T, class U>
    struct at
    {
        typedef Ref entry_ref_t;
        typedef Srcref entry_srcref_t;

        template<class M>
        static Ref get(M& m, T i, U j)
        {
            return mdetail::immediate_traits<M>::at(m, i, j);
        }

        template<class M>
        static Srcref get(const M& m, T i, U j)
        {
            return mdetail::immediate_traits<M>::at(m, i, j);
        }
    };
};
} // matrices

// immediate functions
template<class Mat>
inline typename mp::enable_if<traits::is_mat<Mat>, slong>::type
rank(const Mat& m)
{
    return m.rank();
}

template<class Mat>
inline typename mp::enable_if<traits::is_mat<Mat>, slong>::type
find_pivot_any(const Mat& m, slong start, slong end, slong c)
{
    return m.find_pivot_any(start, end, c);
}

template<class Mat>
inline typename mp::enable_if<traits::is_mat<Mat>, slong>::type
find_pivot_partial(const Mat& m, slong start, slong end, slong c)
{
    return m.find_pivot_partial(start, end, c);
}
} // flint

// Define rows(), cols(), create_temporary() and at() methods.
// For this to work, Traits must behave like the above generic traits.
// Also your matrix class must have a static create_temporary_rowscols function.
// See fmpz_mat for an example.
#define FLINTXX_DEFINE_MATRIX_METHODS(Traits) \
template<class T, class U> \
typename Traits::template at<T, U>::entry_ref_t at(T i, U j) \
    {return Traits::template at<T, U>::get(*this, i, j);} \
template<class T, class U> \
typename Traits::template at<T, U>::entry_srcref_t at(T i, U j) const \
    {return Traits::template at<T, U>::get(*this, i, j);} \
\
slong rows() const {return Traits::rows(*this);} \
slong cols() const {return Traits::cols(*this);} \
evaluated_t create_temporary() const \
{ \
    return create_temporary_rowscols(*this, rows(), cols()); \
}

// Disable temporary merging. Requires create_temporary_rowscols.
// TODO do we really need the ltuple code everywhere?
#define FLINTXX_DEFINE_TEMPORARY_RULES(Matrix) \
namespace traits { \
template<> struct use_temporary_merging<Matrix> : mp::false_ { }; \
} /* traits */ \
namespace rules { \
template<class Expr> \
struct use_default_temporary_instantiation<Expr, Matrix> : mp::false_ { }; \
template<class Expr> \
struct instantiate_temporaries<Expr, Matrix, \
    typename mp::disable_if<flint_classes::is_Base<Matrix, Expr> >::type> \
{ \
    /* The only case where this should ever happen is if Expr is an ltuple */ \
    /* TODO static assert this */ \
    static Matrix get(const Expr& e) \
    { \
        typedef typename Expr::operation_t op_t; \
        return Matrix::create_temporary_rowscols(e, \
                matrices::outsize<op_t>::rows(e), \
                matrices::outsize<op_t>::cols(e)); \
    } \
}; \
template<class Expr> \
struct instantiate_temporaries<Expr, Matrix, \
    typename mp::enable_if<flint_classes::is_Base<Matrix, Expr> >::type> \
{ \
    static Matrix get(const Expr& e) \
    { \
        return e.create_temporary(); \
    } \
}; \
} /* rules */

// Add a fflu() member function to the matrix class.
#define FLINTXX_DEFINE_MEMBER_FFLU \
template<class T> typename detail::nary_op_helper2<operations::fflu_op, \
    typename base_t::derived_t, permxx*, T>::enable::type \
fflu(permxx* p, const T& t) const \
{ \
    return flint::fflu(*this, p, t); \
}

#endif
