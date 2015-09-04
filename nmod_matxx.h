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

#ifndef NMOD_MATXX_H
#define NMOD_MATXX_H

#include <string>
#include <vector>

#include "nmod_mat.h"

#include "nmod_vecxx.h"
#include "fmpz_matxx.h" // for modular reduction
#include "permxx.h"

#include "flintxx/flint_exception.h"
#include "flintxx/ltuple.h"
#include "flintxx/matrix.h"

// TODO addmul
// TODO default argument for mat_solve_triu etc?
// TODO nullspace member
// TODO unnecessary perm copies in set_lu*

namespace flint {
FLINT_DEFINE_BINOP(solve_vec)
FLINT_DEFINE_BINOP(mul_strassen)
FLINT_DEFINE_THREEARY(solve_tril)
FLINT_DEFINE_THREEARY(solve_tril_classical)
FLINT_DEFINE_THREEARY(solve_tril_recursive)
FLINT_DEFINE_THREEARY(solve_triu)
FLINT_DEFINE_THREEARY(solve_triu_classical)
FLINT_DEFINE_THREEARY(solve_triu_recursive)

FLINT_DEFINE_THREEARY(multi_CRT_precomp)

namespace detail {
template<class Mat>
struct nmod_matxx_traits : matrices::generic_traits<Mat> { };
} // detail

template<class Operation, class Data>
class nmod_matxx_expression
    : public expression<derived_wrapper<nmod_matxx_expression>, Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::nmod_matxx_expression>,
              Operation, Data> base_t;
    typedef detail::nmod_matxx_traits<nmod_matxx_expression> traits_t;

    FLINTXX_DEFINE_BASICS(nmod_matxx_expression)
    FLINTXX_DEFINE_CTORS(nmod_matxx_expression)
    FLINTXX_DEFINE_C_REF(nmod_matxx_expression, nmod_mat_struct, _mat)

    // These only make sense with immediates
    nmodxx_ctx_srcref _ctx() const
        {return nmodxx_ctx_srcref::make(_mat()->mod);}

    // These work on any expression without evaluation
    nmodxx_ctx_srcref estimate_ctx() const
    {
        return tools::find_nmodxx_ctx(*this);
    }
    mp_limb_t modulus() const {return estimate_ctx().n();}

    template<class Expr>
    static evaluated_t create_temporary_rowscols(
            const Expr& e, slong rows, slong cols)
    {
        return evaluated_t(rows, cols, tools::find_nmodxx_ctx(e).n());
    }
    FLINTXX_DEFINE_MATRIX_METHODS(traits_t)

    template<class Fmpz_mat>
    static nmod_matxx_expression reduce(const Fmpz_mat& mat,
            mp_limb_t modulus,
            typename mp::enable_if<traits::is_fmpz_matxx<Fmpz_mat> >::type* = 0)
    {
        nmod_matxx_expression res(mat.rows(), mat.cols(), modulus);
        fmpz_mat_get_nmod_mat(res._mat(), mat.evaluate()._mat());
        return res;
    }

    static nmod_matxx_expression randtest(slong rows, slong cols, mp_limb_t n,
            frandxx& state)
    {
        nmod_matxx_expression res(rows, cols, n);
        res.set_randtest(state);
        return res;
    }
    static nmod_matxx_expression randfull(slong rows, slong cols, mp_limb_t n,
            frandxx& state)
    {
        nmod_matxx_expression res(rows, cols, n);
        res.set_randfull(state);
        return res;
    }
    static nmod_matxx_expression randrank(slong rows, slong cols, mp_limb_t n,
            frandxx& state, slong rank)
    {
        nmod_matxx_expression res(rows, cols, n);
        res.set_randrank(state, rank);
        return res;
    }
    static nmod_matxx_expression randtril(slong rows, slong cols, mp_limb_t n,
            frandxx& state, bool unit)
    {
        nmod_matxx_expression res(rows, cols, n);
        res.set_randtril(state, unit);
        return res;
    }
    static nmod_matxx_expression randtriu(slong rows, slong cols, mp_limb_t n,
            frandxx& state,
            bool unit)
    {
        nmod_matxx_expression res(rows, cols, n);
        res.set_randtriu(state, unit);
        return res;
    }

    template<class Vec>
    static nmod_matxx_expression randpermdiag(slong rows, slong cols, mp_limb_t n,
            frandxx& state, const Vec& v)
    {
        nmod_matxx_expression res(rows, cols, n);
        res.set_randpermdiag(state, v);
        return res;
    }

    static nmod_matxx_expression zero(slong rows, slong cols, mp_limb_t n)
        {return nmod_matxx_expression(rows, cols, n);}

    // these only make sense with targets
    void set_randtest(frandxx& state)
        {nmod_mat_randtest(_mat(), state._data());}
    void set_randfull(frandxx& state)
        {nmod_mat_randfull(_mat(), state._data());}
    void set_randrank(frandxx& state, slong rank)
        {nmod_mat_randrank(_mat(), state._data(), rank);}
    void set_randtril(frandxx& state, bool unit)
        {nmod_mat_randtril(_mat(), state._data(), unit);}
    void set_randtriu(frandxx& state, bool unit)
        {nmod_mat_randtriu(_mat(), state._data(), unit);}

    template<class Vec>
    int set_randpermdiag(frandxx& state, const Vec& v)
    {
        return nmod_mat_randpermdiag(_mat(), state._data(), v._array(), v.size());
    }

    void apply_randops(frandxx& state, slong count)
        {nmod_mat_randops(_mat(), count, state._data());}

    slong set_rref() {return nmod_mat_rref(_mat());}
    void set_zero() {nmod_mat_zero(_mat());}

    typedef mp::make_tuple<slong, permxx>::type lu_rt;
    lu_rt set_lu(bool rank_check = false)
    {
        lu_rt res = mp::make_tuple<slong, permxx>::make(0, permxx(rows()));
        res.first() = nmod_mat_lu(res.second()._data(), _mat(), rank_check);
        return res;
    }
    lu_rt set_lu_classical(bool rank_check = false)
    {
        lu_rt res = mp::make_tuple<slong, permxx>::make(0, permxx(rows()));
        res.first() = nmod_mat_lu_classical(
                res.second()._data(), _mat(), rank_check);
        return res;
    }
    lu_rt set_lu_recursive(bool rank_check = false)
    {
        lu_rt res = mp::make_tuple<slong, permxx>::make(0, permxx(rows()));
        res.first() = nmod_mat_lu_recursive(
                res.second()._data(), _mat(), rank_check);
        return res;
    }


    // these cause evaluation
    slong rank() const {return nmod_mat_rank(this->evaluate()._mat());}
    bool is_zero() const {return nmod_mat_is_zero(this->evaluate()._mat());}
    bool is_empty() const {return nmod_mat_is_empty(this->evaluate()._mat());}
    bool is_square() const {return nmod_mat_is_square(this->evaluate()._mat());}

    // lazy members
    FLINTXX_DEFINE_MEMBER_BINOP(solve)
    FLINTXX_DEFINE_MEMBER_BINOP(mul_classical)
    FLINTXX_DEFINE_MEMBER_BINOP(mul_strassen)
    FLINTXX_DEFINE_MEMBER_UNOP(inv)
    FLINTXX_DEFINE_MEMBER_UNOP(transpose)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(nmodxx, trace)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(nmodxx, det)
    //FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(???, nullspace) // TODO
    FLINTXX_DEFINE_MEMBER_3OP(solve_tril)
    FLINTXX_DEFINE_MEMBER_3OP(solve_tril_recursive)
    FLINTXX_DEFINE_MEMBER_3OP(solve_tril_classical)
    FLINTXX_DEFINE_MEMBER_3OP(solve_triu)
    FLINTXX_DEFINE_MEMBER_3OP(solve_triu_recursive)
    FLINTXX_DEFINE_MEMBER_3OP(solve_triu_classical)
};

namespace detail {
struct nmod_mat_data;
} // detail

typedef nmod_matxx_expression<operations::immediate, detail::nmod_mat_data> nmod_matxx;
typedef nmod_matxx_expression<operations::immediate,
            flint_classes::ref_data<nmod_matxx, nmod_mat_struct> > nmod_matxx_ref;
typedef nmod_matxx_expression<operations::immediate, flint_classes::srcref_data<
    nmod_matxx, nmod_matxx_ref, nmod_mat_struct> > nmod_matxx_srcref;

template<>
struct matrix_traits<nmod_matxx>
{
    template<class M> static slong rows(const M& m)
    {
        return nmod_mat_nrows(m._mat());
    }
    template<class M> static slong cols(const M& m)
    {
        return nmod_mat_ncols(m._mat());
    }

    template<class M> static nmodxx_srcref at(const M& m, slong i, slong j)
    {
        return nmodxx_srcref::make(nmod_mat_entry(m._mat(), i, j),
                m.estimate_ctx());
    }
    template<class M> static nmodxx_ref at(M& m, slong i, slong j)
    {
        return nmodxx_ref::make(nmod_mat_entry(m._mat(), i, j),
                m.estimate_ctx());
    }
};

namespace traits {
template<> struct has_nmodxx_ctx<nmod_matxx> : mp::true_ { };
template<> struct has_nmodxx_ctx<nmod_matxx_ref> : mp::true_ { };
template<> struct has_nmodxx_ctx<nmod_matxx_srcref> : mp::true_ { };
} // traits

namespace detail {
template<>
struct nmod_matxx_traits<nmod_matxx_srcref>
    : matrices::generic_traits_srcref<nmodxx_srcref> { };
template<>
struct nmod_matxx_traits<nmod_matxx_ref>
    : matrices::generic_traits_ref<nmodxx_ref> { };
template<> struct nmod_matxx_traits<nmod_matxx>
    : matrices::generic_traits_nonref<nmodxx_ref, nmodxx_srcref> { };

struct nmod_mat_data
{
    typedef nmod_mat_t& data_ref_t;
    typedef const nmod_mat_t& data_srcref_t;

    nmod_mat_t inner;

    nmod_mat_data(slong m, slong n, mp_limb_t modulus)
    {
        nmod_mat_init(inner, m, n, modulus);
    }

    nmod_mat_data(const nmod_mat_data& o)
    {
        nmod_mat_init_set(inner, o.inner);
    }

    nmod_mat_data(nmod_matxx_srcref o)
    {
        nmod_mat_init_set(inner, o._data().inner);
    }

    ~nmod_mat_data() {nmod_mat_clear(inner);}
};
} // detail

namespace matrices {
template<>
struct outsize<operations::mul_strassen_op>
    : outsize<operations::times> { };

template<> struct outsize<operations::solve_tril_op>
    : outsize<operations::solve_op> { };
template<> struct outsize<operations::solve_tril_classical_op>
    : outsize<operations::solve_op> { };
template<> struct outsize<operations::solve_tril_recursive_op>
    : outsize<operations::solve_op> { };
template<> struct outsize<operations::solve_triu_op>
    : outsize<operations::solve_op> { };
template<> struct outsize<operations::solve_triu_classical_op>
    : outsize<operations::solve_op> { };
template<> struct outsize<operations::solve_triu_recursive_op>
    : outsize<operations::solve_op> { };
}

// temporary instantiation stuff
FLINTXX_DEFINE_TEMPORARY_RULES(nmod_matxx)

#define NMOD_MATXX_COND_S FLINTXX_COND_S(nmod_matxx)
#define NMOD_MATXX_COND_T FLINTXX_COND_T(nmod_matxx)

namespace traits {
template<class T> struct is_nmod_matxx
    : flint_classes::is_Base<nmod_matxx, T> { };
} // traits

namespace rules {
FLINT_DEFINE_DOIT_COND2(assignment, NMOD_MATXX_COND_T, NMOD_MATXX_COND_S,
        nmod_mat_set(to._mat(), from._mat()))

FLINTXX_DEFINE_SWAP(nmod_matxx, nmod_mat_swap(e1._mat(), e2._mat()))

FLINTXX_DEFINE_EQUALS(nmod_matxx, nmod_mat_equal(e1._mat(), e2._mat()))

FLINT_DEFINE_PRINT_PRETTY_COND(NMOD_MATXX_COND_S,
        (nmod_mat_print_pretty(from._mat()), 1))

FLINT_DEFINE_THREEARY_EXPR_COND3(mat_at_op, nmodxx,
        NMOD_MATXX_COND_S, traits::fits_into_slong, traits::fits_into_slong,
        to.set_nored(nmod_mat_entry(e1._mat(), e2, e3)))

FLINT_DEFINE_BINARY_EXPR_COND2(times, nmod_matxx,
        NMOD_MATXX_COND_S, NMOD_MATXX_COND_S,
        nmod_mat_mul(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, nmod_matxx,
        NMOD_MATXX_COND_S, NMODXX_COND_S,
        nmod_mat_scalar_mul(to._mat(), e1._mat(), e2._limb()))

FLINT_DEFINE_BINARY_EXPR_COND2(plus, nmod_matxx,
        NMOD_MATXX_COND_S, NMOD_MATXX_COND_S,
        nmod_mat_add(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(minus, nmod_matxx,
        NMOD_MATXX_COND_S, NMOD_MATXX_COND_S,
        nmod_mat_sub(to._mat(), e1._mat(), e2._mat()))

FLINT_DEFINE_UNARY_EXPR_COND(negate, nmod_matxx, NMOD_MATXX_COND_S,
        nmod_mat_neg(to._mat(), from._mat()))

FLINT_DEFINE_UNARY_EXPR_COND(transpose_op, nmod_matxx, NMOD_MATXX_COND_S,
        nmod_mat_transpose(to._mat(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(trace_op, nmodxx, NMOD_MATXX_COND_S,
        to.set_nored(nmod_mat_trace(from._mat())))

FLINT_DEFINE_BINARY_EXPR_COND2(mul_classical_op, nmod_matxx,
        NMOD_MATXX_COND_S, NMOD_MATXX_COND_S,
        nmod_mat_mul(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(mul_strassen_op, nmod_matxx,
        NMOD_MATXX_COND_S, NMOD_MATXX_COND_S,
        nmod_mat_mul(to._mat(), e1._mat(), e2._mat()))

FLINT_DEFINE_UNARY_EXPR_COND(det_op, nmodxx, NMOD_MATXX_COND_S,
        to.set_nored(nmod_mat_det(from._mat())))

FLINT_DEFINE_UNARY_EXPR_COND(inv_op, nmod_matxx, NMOD_MATXX_COND_S,
        execution_check(nmod_mat_inv(to._mat(), from._mat()),
            "inv", "nmod_mat"))

#define NMOD_MATXX_DEFINE_SOLVE_TRI(name) \
FLINT_DEFINE_THREEARY_EXPR_COND3(name##_op, nmod_matxx, \
        NMOD_MATXX_COND_S, NMOD_MATXX_COND_S, tools::is_bool, \
        nmod_mat_##name(to._mat(), e1._mat(), e2._mat(), e3))
NMOD_MATXX_DEFINE_SOLVE_TRI(solve_tril)
NMOD_MATXX_DEFINE_SOLVE_TRI(solve_tril_classical)
NMOD_MATXX_DEFINE_SOLVE_TRI(solve_tril_recursive)
NMOD_MATXX_DEFINE_SOLVE_TRI(solve_triu)
NMOD_MATXX_DEFINE_SOLVE_TRI(solve_triu_classical)
NMOD_MATXX_DEFINE_SOLVE_TRI(solve_triu_recursive)

FLINT_DEFINE_BINARY_EXPR_COND2(solve_op, nmod_matxx,
        NMOD_MATXX_COND_S, NMOD_MATXX_COND_S,
        execution_check(nmod_mat_solve(to._mat(), e1._mat(), e2._mat()),
            "solve", "nmod_mat"))

FLINT_DEFINE_BINARY_EXPR_COND2(solve_op, nmod_vecxx,
        NMOD_MATXX_COND_S, NMOD_VECXX_COND_S,
        execution_check(nmod_mat_solve_vec(to._array(), e1._mat(), e2._array()),
            "solve_vec", "nmod_mat"))

namespace rdetail {
typedef make_ltuple<mp::make_tuple<slong, nmod_matxx>::type >::type
    nmod_mat_nullspace_rt;
} // rdetail
FLINT_DEFINE_UNARY_EXPR_COND(nullspace_op, rdetail::nmod_mat_nullspace_rt,
        NMOD_MATXX_COND_S, to.template get<0>() = nmod_mat_nullspace(
            to.template get<1>()._mat(), from._mat()))
} // rules

//////////////////////////////////////////////////////////////////////////////
// nmod_mat_vector class
//////////////////////////////////////////////////////////////////////////////
// This class stores a vector of nmod_matxx with differing moduli. It is *not*
// an expression template class!

class nmod_mat_vector
{
private:
    nmod_mat_t* data;
    std::size_t size_;

    void init(const nmod_mat_vector& o)
    {
        size_ = o.size_;
        data = new nmod_mat_t[size_];
        for(std::size_t i = 0;i < size_;++i)
            nmod_mat_init_set(data[i], o.data[i]);
    }

public:
    ~nmod_mat_vector() {delete[] data;}
    nmod_mat_vector(slong rows, slong cols, const std::vector<mp_limb_t>& primes)
    {
        size_ = primes.size();
        data = new nmod_mat_t[primes.size()];
        for(std::size_t i = 0;i < primes.size();++i)
            nmod_mat_init(data[i], rows, cols, primes[i]);
    }

    nmod_mat_vector(const nmod_mat_vector& o)
    {
        init(o);
    }

    nmod_mat_vector& operator=(const nmod_mat_vector& o)
    {
        delete[] data;
        init(o);
        return *this;
    }

    nmod_matxx_ref operator[](std::size_t idx)
        {return nmod_matxx_ref::make(data[idx]);}
    nmod_matxx_srcref operator[](std::size_t idx) const
        {return nmod_matxx_srcref::make(data[idx]);}

    std::size_t size() const {return size_;}

    const nmod_mat_t* _data() const {return data;}
    nmod_mat_t* _data() {return data;}

    bool operator==(const nmod_mat_vector& o)
    {
        if(size() != o.size())
            return false;
        for(std::size_t i = 0;i < size();++i)
            if((*this)[i] != o[i])
                return false;
        return true;
    }
    bool operator!=(const nmod_mat_vector& o)
    {
        return !(*this == o);
    }

    template<class Fmpz_mat>
    void set_multi_mod(const Fmpz_mat& m,
            typename mp::enable_if<traits::is_fmpz_matxx<Fmpz_mat> >::type* = 0)
    {
        fmpz_mat_multi_mod_ui(data, size(), m.evaluate()._mat());
    }
    template<class Fmpz_mat>
    void set_multi_mod_precomp(const Fmpz_mat& m,
            const fmpz_combxx& comb,
            typename mp::enable_if<traits::is_fmpz_matxx<Fmpz_mat> >::type* = 0)
    {
        fmpz_mat_multi_mod_ui_precomp(data, size(), m.evaluate()._mat(),
                comb._comb(), comb._temp());
    }
};

/////////////////////////////////////////////////////////////////////////////
// chinese remaindering
/////////////////////////////////////////////////////////////////////////////
// Note this operates on fmpz_matxx and fmpz_combxx (as well as nmod_matxx).
// We define it here to deal with the circular dependencies. fmpz_matxx.h
// includes nmod_matxx.h at the bottom.

template<class Fmpz_mat>
inline nmod_mat_vector multi_mod(const Fmpz_mat& m,
        const std::vector<mp_limb_t>& primes,
        typename mp::enable_if<traits::is_fmpz_matxx<Fmpz_mat> >::type* = 0)
{
    nmod_mat_vector res(m.rows(), m.cols(), primes);
    res.set_multi_mod(m);
    return res;
}
template<class Fmpz_mat>
inline nmod_mat_vector multi_mod_precomp(const Fmpz_mat& m,
        const std::vector<mp_limb_t>& primes,
        const fmpz_combxx& comb,
        typename mp::enable_if<traits::is_fmpz_matxx<Fmpz_mat> >::type* = 0)
{
    nmod_mat_vector res(m.rows(), m.cols(), primes);
    res.set_multi_mod_precomp(m, comb);
    return res;
}

namespace matrices {
// outsize computation for multi-CRT
struct outsize_CRT
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return m._data().first()[0].rows();
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return m._data().first()[0].cols();
    }
};

template<> struct outsize<operations::multi_CRT_op> : outsize_CRT { };
template<> struct outsize<operations::multi_CRT_precomp_op> : outsize_CRT { };
}

namespace rules {
FLINT_DEFINE_FOURARY_EXPR_COND4(CRT_op, fmpz_matxx,
        FMPZ_MATXX_COND_T, FMPZXX_COND_S, NMOD_MATXX_COND_S, tools::is_bool,
        fmpz_mat_CRT_ui(to._mat(), e1._mat(), e2._fmpz(), e3._mat(), e4))

FLINT_DEFINE_BINARY_EXPR2(multi_CRT_op, fmpz_matxx, nmod_mat_vector, bool,
        fmpz_mat_multi_CRT_ui(to._mat(), (nmod_mat_t * const) e1._data(), e1.size(), e2))
FLINT_DEFINE_THREEARY_EXPR(multi_CRT_precomp_op, fmpz_matxx,
        nmod_mat_vector, fmpz_combxx, bool,
        fmpz_mat_multi_CRT_ui_precomp(to._mat(), (nmod_mat_t * const) e1._data(), e1.size(),
            e2._comb(), e2._temp(), e3))
} // rules
} // flint

#endif
