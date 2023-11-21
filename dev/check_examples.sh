#!/usr/bin/env bash

# Copyright (C) 2023 Albin Ahlbäck
#
# This file is part of FLINT.
#
# FLINT is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License (LGPL) as published
# by the Free Software Foundation; either version 2.1 of the License, or
# (at your option) any later version.  See <https://www.gnu.org/licenses/>.

if test "$1" = "bernoulli";
then
    echo -n "bernoulli...."
    res=$($2/bernoulli 20 -threads 2)
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 1
    fi
    echo $res | perl -0ne 'if (/-174611\/330/) { $found=1; last } END { exit !$found }'
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 2
    fi
    echo "PASS"
    exit 0
elif test "$1" = "binet";
then
    echo -n "binet...."
    res=$($2/binet -limit 2 8)
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 1
    fi
    echo $res | perl -0ne 'if (/21 cpu/) { $found=1; last } END { exit !$found }'
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 2
    fi
    echo "PASS"
    exit 0
elif test "$1" = "class_poly";
then
    echo "class_poly....SKIPPED"
    exit 0
elif test "$1" = "complex_plot";
then
    echo "complex_plot....SKIPPED"
    exit 0
elif test "$1" = "crt";
then
    echo "crt....SKIPPED"
    exit 0
elif test "$1" = "delta_qexp";
then
    echo "delta_qexp....SKIPPED"
    exit 0
elif test "$1" = "dft";
then
    echo "dft....SKIPPED"
    exit 0
elif test "$1" = "elementary";
then
    echo -n "elementary...."
    res=$($2/elementary)
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 1
    fi
    # Check that the first and last expression prints as intended
    echo "$res" | perl -0ne 'if (/>>> Exp\(Pi\*I\) \+ 1\n0\n/) { $found=1; last } END { exit !$found }'
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 2
    fi
    echo "$res" | perl -0ne 'if (/>>> Erf\(2\*Log\(Sqrt\(1\/2-Sqrt\(2\)\/4\)\)\+Log\(4\)\) - Erf\(Log\(2-Sqrt\(2\)\)\)\n0\n/) { $found=1; last } END { exit !$found }'
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 2
    fi
    echo "PASS"
    exit 0
elif test "$1" = "factor_integer";
then
    echo -n "factor_integer...."
    res=$($2/factor_integer -threads 3 144)
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 1
    fi
    echo "$res" | perl -0ne 'if (/144 =\n2\^4 \* 3\^2/) { $found=1; last } END { exit !$found }'
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 2
    fi
    echo "PASS"
    exit 0
elif test "$1" = "factor_polynomial";
then
    echo -n "factor_polynomial...."
    res=$($2/factor_polynomial -threads 2 "4*x*(1+x^2+y^2*x^3+z^4)^2")
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 1
    fi
    echo "$res" | perl -0ne 'if (/4\*x\^7\*y\^4\+8\*x\^6\*y\^2\+4\*x\^5\+8\*x\^4\*y\^2\*z\^4\+8\*x\^4\*y\^2\+8\*x\^3\*z\^4\+8\*x\^3\+4\*x\*z\^8\+8\*x\*z\^4\+4\*x =\n4\*\(x\)\^1\*\(x\^3\*y\^2\+x\^2\+z\^4\+1\)\^2/) { $found=1; last } END { exit !$found }'
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 2
    fi
    echo "PASS"
    exit 0
elif test "$1" = "fmpq_poly";
then
    echo "fmpq_poly....SKIPPED"
elif test "$1" = "fmpz_mod_poly";
then
    echo "fmpz_mod_poly....SKIPPED"
elif test "$1" = "fmpz_poly_factor_zassenhaus";
then
    echo "fmpz_poly_factor_zassenhaus....SKIPPED"
elif test "$1" = "fmpz_poly_q";
then
    echo "fmpz_poly_q....SKIPPED"
elif test "$1" = "fpwrap";
then
    echo "fpwrap....SKIPPED"
elif test "$1" = "fq_poly";
then
    echo "fq_poly....SKIPPED"
elif test "$1" = "functions_benchmark";
then
    echo "functions_benchmark....SKIPPED"
elif test "$1" = "hilbert_matrix";
then
    echo -n "hilbert_matrix...."
    res=$($2/hilbert_matrix -eig 20)
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 1
    fi
    echo "$res" | perl -0ne 'if (/prec=20: nan\nprec=40: nan\nprec=80: nan\nprec=160: nan\nprec=320: \[7\.777377397e-29 \+\/- 1.44e-39\]\nsuccess!\n/) { $found=1; last } END { exit !$found }'
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 2
    fi
    echo "PASS"
    exit 0
elif test "$1" = "hilbert_matrix_ca";
then
    echo "hilbert_matrix_ca....SKIPPED"
elif test "$1" = "huge_expr";
then
    echo "huge_expr....SKIPPED"
elif test "$1" = "integrals";
then
    echo -n "integrals...."
    res=$($2/integrals -i 21 -threads 3)
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 1
    fi
    echo "$res" | perl -0ne 'if (/I21 = \[43\.3273391411077 \+\/- 4\.52e-14] \+ \[\+\/- 1\.95e-14\]\*I/) { $found=1; last } END { exit !$found }'
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 2
    fi
    echo "PASS"
    exit 0
elif test "$1" = "keiper_li";
then
    echo "keiper_li....SKIPPED"
elif test "$1" = "lcentral";
then
    echo "lcentral....SKIPPED"
elif test "$1" = "logistic";
then
    echo "logistic....SKIPPED"
elif test "$1" = "lvalue";
then
    echo "lvalue....SKIPPED"
elif test "$1" = "machin";
then
    echo "machin....SKIPPED"
elif test "$1" = "multi_crt";
then
    echo "multi_crt....SKIPPED"
elif test "$1" = "padic";
then
    echo "padic....SKIPPED"
elif test "$1" = "partitions";
then
    echo "partitions....SKIPPED"
elif test "$1" = "pi";
then
    echo -n "pi...."
    res=$($2/pi 1000 -threads 2)
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 1
    fi
    echo "$res" | perl -0ne 'if (/\[3\.14159265358979323846\{\.\.\.959 digits\.\.\.\}76611195909216420199 \+\/- 8\.09e-1001\]/) { $found=1; last } END { exit !$found }'
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 2
    fi
    echo "PASS"
    exit 0
elif test "$1" = "poly_roots";
then
    echo -n "poly_roots...."
    res=$($2/poly_roots t 10)
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 1
    fi
    echo "$res" | perl -0ne 'if (/\n10 roots with multiplicity 1\nsearching for 10 roots, 5 deflated\nprec=32: 5 isolated roots/) { $found=1; last } END { exit !$found }'
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 2
    fi
    echo "PASS"
    exit 0
elif test "$1" = "primegen";
then
    echo "primegen....SKIPPED"
elif test "$1" = "qadic";
then
    echo "qadic....SKIPPED"
elif test "$1" = "radix";
then
    echo "radix....SKIPPED"
elif test "$1" = "real_roots";
then
    echo "real_roots....SKIPPED"
elif test "$1" = "stirling_matrix";
then
    echo "stirling_matrix....SKIPPED"
elif test "$1" = "swinnerton_dyer_poly";
then
    echo "swinnerton_dyer_poly....SKIPPED"
elif test "$1" = "taylor_integrals";
then
    echo "taylor_integrals....SKIPPED"
elif test "$1" = "zeta_zeros";
then
    echo -n "zeta_zeros...."
    res=$($2/zeta_zeros -count 100 -threads 3)
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 1
    fi
    echo "$res" | perl -0ne 'if (/97\t231\.2501887004991648\n98\t231\.9872352531802486\n99\t233\.6934041789083006\n100\t236\.5242296658162058/) { $found=1; last } END { exit !$found }'
    if test "$?" != "0";
    then
        echo "FAIL"
        exit 2
    fi
    echo "PASS"
    exit 0
else
    exit 3
fi
