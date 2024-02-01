#
#   Copyright (C) 2024 Albin Ahlbäck
#
#   This file is part of FLINT.
#
#   FLINT is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Lesser General Public License (LGPL) as published
#   by the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#

###############################################################################
###############################################################################
# helper functions
###############################################################################
###############################################################################

# Example:
#
# str = "19,31,21"
# new_str, ix = get_comma_separated_item("19,31,21")
# # new_str = "19"
# # ix = 4
# # str[ix:end] = "31,21"
function get_comma_separated_item(str::String)
    ix = 2
    while str[ix] != ','
        ix += 1
    end
    return (str[1:ix - 1], ix + 1)
end

function push_coefficients(op::String, ip::Vector{UInt}, num::UInt,
        finalseparator::String = ",")
    for ix in 1:num
        op *= "$(ip[ix]),"
    end
    return op[1:end - 1] * finalseparator
end

###############################################################################
# prime < 260
###############################################################################

# Idea: Search for prime in __flint_cp_primes0, then search for corresponding
# degree in __flint_cp_degrees0, then search for non-trivial coefficients in
# __flint_ntcoeffs0 by using __flint_numntcoeffs0.

function primes_lt_260(lines::Vector{String})
    # Lines are on the form "p, deg, c_{0}, c_{1}, ..., c_{deg}"
    coeffs = Vector{UInt}(undef, 410)

    # Current line
    lx = 1

    primes = "uint8_t __flint_cp_primes0[] = {"
    degrees = "uint16_t __flint_cp_degrees0[] = {"
    numntcoeffs = "uint8_t __flint_numntcoeffs0[] = {"
    ntcoeffs = "uint8_t __flint_ntcoeffs0[] = {"

    # Get first prime
    ps, ix = get_comma_separated_item(lines[lx])
    lines[lx] = lines[lx][ix:end]
    prime = parse(UInt, ps)
    primes *= "$(prime - 2),"

    # Check that we start with the correct one
    ds, _ = get_comma_separated_item(lines[lx])
    degree = parse(UInt, ds)
    if prime != 2 || degree != 1
        error()
    end

    while prime < 260
        # Get degree
        ds, ix = get_comma_separated_item(lines[lx])
        lines[lx] = lines[lx][ix:end]
        degree = parse(UInt, ds)
        degrees *= "$ds,"

        # Get coefficients
        for jx in 0:degree
            cs, ix = get_comma_separated_item(lines[lx])
            lines[lx] = lines[lx][ix:end]
            coeffs[jx + 1] = parse(UInt, cs)
        end

        # Check that leading coefficient indeed is one
        if coeffs[degree + 1] != 1
            error()
        end

        # Non-trivial coefficients
        for jx in degree:-1:1
            if coeffs[jx] != 0
                numntcoeffs *= "$jx,"
                ntcoeffs = push_coefficients(ntcoeffs, coeffs, jx)
                break
            end
        end

        # Increment one line in file
        lx += 1

        # Get new prime
        ps, ix = get_comma_separated_item(lines[lx])
        newprime = parse(UInt, ps)
        if newprime != prime
            prime = newprime
            if newprime < 260
                primes *= "$(newprime - 2),"
            end
        end
        lines[lx] = lines[lx][ix:end]
    end

    primes = primes[1:end] * "0};"
    degrees = degrees[1:end] * "0};"
    numntcoeffs = numntcoeffs[1:end - 1] * "};"
    ntcoeffs = ntcoeffs[1:end - 1] * "};"

    return lx - 1, primes, degrees, numntcoeffs, ntcoeffs
end

###############################################################################
# 260 < prime < 300
###############################################################################

function check_deg_1_300(coeffs::Vector{UInt}, a0::UInt, b0::UInt)
    if (coeffs[1] != b0 || coeffs[2] != 1)
        error("a0 = $a0, b0 = $b0\ncoeffs = $coeffs")
    end
end

function check_deg_2_300(coeffs::Vector{UInt}, a0::UInt, b0::UInt)
    if (coeffs[1] != a0 || coeffs[2] >= 2^16 || coeffs[3] != 1)
        error("a0 = $a0, b0 = $b0\ncoeffs = $coeffs")
    end
end

function check_deg_3_300(coeffs::Vector{UInt}, a0::UInt, b0::UInt)
    if (coeffs[1] != b0 || coeffs[2] >= 2^8 || coeffs[3] != 0
        || coeffs[4] != 1)
        error("a0 = $a0, b0 = $b0\ncoeffs = $coeffs")
    end
end

function check_deg_4_300(coeffs::Vector{UInt}, a0::UInt, b0::UInt)
    if (coeffs[1] != a0 || coeffs[2] >= 2^16 || coeffs[3] >= 2^8 ||
        coeffs[4] != 0 || coeffs[5] != 1)
        error("a0 = $a0, b0 = $b0\ncoeffs = $coeffs")
    end
end

function check_deg_5_300(coeffs::Vector{UInt}, a0::UInt, b0::UInt)
    if (coeffs[1] != b0 || coeffs[2] >= 2^8 || coeffs[3] != 0 ||
        coeffs[4] != 0 || coeffs[5] != 0 || coeffs[6] != 1)
        error("a0 = $a0, b0 = $b0\ncoeffs = $coeffs")
    end
end

function check_deg_6_300(coeffs::Vector{UInt}, a0::UInt, b0::UInt)
    if (coeffs[1] != a0 || coeffs[2] >= 2^16 || coeffs[3] >= 2^8 ||
        coeffs[4] >= 2^8 || coeffs[5] >= 2^8 || coeffs[6] != 0 ||
        coeffs[7] != 1)
        error("a0 = $a0, b0 = $b0\ncoeffs = $coeffs")
    end
end

function check_deg_7_300(coeffs::Vector{UInt}, a0::UInt, b0::UInt)
    if (coeffs[1] != b0 || coeffs[2] >= 2^8 || coeffs[3] != 0 ||
        coeffs[4] != 0 || coeffs[5] != 0 || coeffs[6] != 0 ||
        coeffs[7] != 0 || coeffs[8] != 1)
        error("a0 = $a0, b0 = $b0\ncoeffs = $coeffs")
    end
end

function check_deg_8_300(coeffs::Vector{UInt}, a0::UInt, b0::UInt)
    if (coeffs[1] != a0 || coeffs[2] >= 2^8 || coeffs[3] >= 2^16 ||
        coeffs[4] >= 2^8 || coeffs[5] >= 2^8 || coeffs[6] != 0 ||
        coeffs[7] != 0 || coeffs[8] != 0 || coeffs[9] != 1)
        error("a0 = $a0, b0 = $b0\ncoeffs = $coeffs")
    end
end

function check_deg_9_300(coeffs::Vector{UInt}, a0::UInt, b0::UInt)
    if (coeffs[1] != b0 || coeffs[2] >= 2^16 || coeffs[3] >= 2^16 ||
        coeffs[4] >= 2^8 || coeffs[5] != 0 || coeffs[6] != 0 ||
        coeffs[7] != 0 || coeffs[8] != 0 || coeffs[9] != 0 ||
        coeffs[10] != 1)
        error("a0 = $a0, b0 = $b0\ncoeffs = $coeffs")
    end
end

function check_deg_10_300(coeffs::Vector{UInt}, a0::UInt, b0::UInt)
    if (coeffs[1] != a0 || coeffs[2] >= 2^16 || coeffs[3] >= 2^8 ||
        coeffs[4] >= 2^16 || coeffs[5] >= 2^8 || coeffs[6] >= 2^16 ||
        coeffs[7] >= 2^8 || coeffs[8] != 0 || coeffs[9] != 0 ||
        coeffs[10] != 0 || coeffs[11] != 1)
        error("a0 = $a0, b0 = $b0\ncoeffs = $coeffs")
    end
end

function check_deg_11_300(coeffs::Vector{UInt}, a0::UInt, b0::UInt)
    if (coeffs[1] != b0 || coeffs[2] >= 2^8 || coeffs[3] != 0 ||
        coeffs[4] != 0 || coeffs[5] != 0 || coeffs[6] != 0 ||
        coeffs[7] != 0 || coeffs[8] != 0 || coeffs[9] != 0 ||
        coeffs[10] != 0 || coeffs[11] != 0 || coeffs[12] != 1)
        error("a0 = $a0, b0 = $b0\ncoeffs = $coeffs")
    end
end

function check_deg_12_300(coeffs::Vector{UInt}, a0::UInt, b0::UInt)
    if (coeffs[1] != a0 || coeffs[2] >= 2^8 || coeffs[3] >= 2^16 ||
        coeffs[4] >= 2^8 || coeffs[5] >= 2^8 || coeffs[6] >= 2^8 ||
        coeffs[7] >= 2^8 || coeffs[8] >= 2^8 || coeffs[9] >= 2^8 ||
        coeffs[10] != 0 || coeffs[11] != 0 || coeffs[12] != 0 ||
        coeffs[13] != 1)
        error("a0 = $a0, b0 = $b0\ncoeffs = $coeffs")
    end
end

# [b0,   1],
# [a0,  b1,   1],
# [b0,  a1,   0,   1],
# [a0,  b2,  a2,   0,   1],
# [b0,  a3,   0,   0,   0,   1],
# [a0,  b3,  a4,  a5,  a6,   0,   1],
# [b0,  a7,   0,   0,   0,   0,   0,   1],
# [a0,  a8,  b4,  a9, a10,   0,   0,   0,   1],
# [b0,  b5,  b6, a11,   0,   0,   0,   0,   0, 1],
# [a0,  b7, a12,  b8, a13,  b9, a14,   0,   0, 0, 1],
# [b0, a15,   0,   0,   0,   0,   0,   0,   0, 0, 0, 1],
# [a0, a16, b10, a17, a18, a19, a20, a21, a22, 0, 0, 0, 1],

function check_coefficients_300(coeffs::Vector{Vector{UInt}})
    a0 = coeffs[2][1]
    b0 = coeffs[1][1]

    check_deg_1_300(coeffs[1], a0, b0)
    check_deg_2_300(coeffs[2], a0, b0)
    check_deg_3_300(coeffs[3], a0, b0)
    check_deg_4_300(coeffs[4], a0, b0)
    check_deg_5_300(coeffs[5], a0, b0)
    check_deg_6_300(coeffs[6], a0, b0)
    check_deg_7_300(coeffs[7], a0, b0)
    check_deg_8_300(coeffs[8], a0, b0)
    check_deg_9_300(coeffs[9], a0, b0)
    check_deg_10_300(coeffs[10], a0, b0)
    check_deg_11_300(coeffs[11], a0, b0)
    check_deg_12_300(coeffs[12], a0, b0)
end

function primes_lt_300(lines::Vector{String})
    coeffs = [Vector{UInt}(undef, ix + 1) for ix in 1:12]

    # Current line
    lx = 1

    primes = "uint16_t __flint_cp_primes1[] = {"
    sm_coeffs = "uint8_t __flint_cp_sm_coeffs1[] = {"
    md_coeffs = "uint16_t __flint_cp_md_coeffs1[] = {"

    # Get first prime
    ps, ix = get_comma_separated_item(lines[lx])
    prime = parse(UInt, ps)

    # Check that we start with the correct one
    ds, _ = get_comma_separated_item(lines[lx][ix:end])
    degree = parse(UInt, ds)
    if prime != 263 || degree != 1
        error()
    end

    while prime < 300
        primes *= "$(prime),"

        for degree in 1:12
            ps, ix = get_comma_separated_item(lines[lx])
            thisprime = parse(UInt, ps)
            lines[lx] = lines[lx][ix:end]

            if thisprime != prime
                error("lx = $lx, degree = $degree")
            end

            ds, ix = get_comma_separated_item(lines[lx])
            lines[lx] = lines[lx][ix:end]
            thisdegree = parse(UInt, ds)

            if thisdegree != degree
                error()
            end

            for jx in 0:degree
                cs, ix = get_comma_separated_item(lines[lx])
                lines[lx] = lines[lx][ix:end]
                coeffs[degree][jx + 1] = parse(UInt, cs)
            end

            lx += 1
        end

        check_coefficients_300(coeffs)

        sm_coeffs *= "$(coeffs[2][1]),"
        sm_coeffs *= "$(coeffs[3][2]),"
        sm_coeffs *= "$(coeffs[4][3]),"
        sm_coeffs *= "$(coeffs[5][2]),"
        sm_coeffs *= "$(coeffs[6][3]),$(coeffs[6][4]),$(coeffs[6][5]),"
        sm_coeffs *= "$(coeffs[7][2]),"
        sm_coeffs *= "$(coeffs[8][2]),$(coeffs[8][4]),$(coeffs[8][5]),"
        sm_coeffs *= "$(coeffs[9][4]),"
        sm_coeffs *= "$(coeffs[10][3]),$(coeffs[10][5]),$(coeffs[10][7]),"
        sm_coeffs *= "$(coeffs[11][2]),"
        sm_coeffs *= "$(coeffs[12][2]),$(coeffs[12][4]),$(coeffs[12][5]),$(coeffs[12][6]),$(coeffs[12][7]),$(coeffs[12][8]),$(coeffs[12][9]),"

        md_coeffs *= "$(coeffs[1][1]),"
        md_coeffs *= "$(coeffs[2][2]),"
        md_coeffs *= "$(coeffs[4][2]),"
        md_coeffs *= "$(coeffs[6][2]),"
        md_coeffs *= "$(coeffs[8][3]),"
        md_coeffs *= "$(coeffs[9][2]),$(coeffs[9][3]),"
        md_coeffs *= "$(coeffs[10][2]),$(coeffs[10][4]),$(coeffs[10][6]),"
        md_coeffs *= "$(coeffs[12][3]),"

        # Get new prime
        ps, _ = get_comma_separated_item(lines[lx])
        newprime = parse(UInt, ps)
        if newprime == prime
            error()
        end
        prime = newprime
    end

    primes = primes[1:end] * "0};"
    sm_coeffs = sm_coeffs[1:end - 1] * "};"
    md_coeffs = md_coeffs[1:end - 1] * "};"

    return lx - 1, primes, sm_coeffs, md_coeffs
end

###############################################################################
# 300 < prime < 2^16
###############################################################################

check_deg_1(coeffs::Vector{UInt}, a0::UInt, b0::UInt) =
    check_deg_1_300(coeffs, a0, b0)

check_deg_2(coeffs::Vector{UInt}, a0::UInt, b0::UInt) =
    check_deg_2_300(coeffs, a0, b0)

check_deg_3(coeffs::Vector{UInt}, a0::UInt, b0::UInt) =
    check_deg_3_300(coeffs, a0, b0)

check_deg_4(coeffs::Vector{UInt}, a0::UInt, b0::UInt) =
    check_deg_4_300(coeffs, a0, b0)

check_deg_5(coeffs::Vector{UInt}, a0::UInt, b0::UInt) =
    check_deg_5_300(coeffs, a0, b0)

function check_deg_6(coeffs::Vector{UInt}, a0::UInt, b0::UInt)
    if (coeffs[1] != a0 || coeffs[2] >= 2^16 || coeffs[3] >= 2^16 ||
        coeffs[4] >= 2^16 || coeffs[5] >= 2^8 || coeffs[6] != 0 ||
        coeffs[7] != 1)
        error("a0 = $a0, b0 = $b0\ncoeffs = $coeffs")
    end
end

check_deg_7(coeffs::Vector{UInt}, a0::UInt, b0::UInt) =
    check_deg_7_300(coeffs, a0, b0)

function check_deg_8(coeffs::Vector{UInt}, a0::UInt, b0::UInt)
    if (coeffs[1] != a0 || coeffs[2] >= 2^16 || coeffs[3] >= 2^16 ||
        coeffs[4] >= 2^16 || coeffs[5] >= 2^8 || coeffs[6] != 0 ||
        coeffs[7] != 0 || coeffs[8] != 0 || coeffs[9] != 1)
        error("a0 = $a0, b0 = $b0\ncoeffs = $coeffs")
    end
end

check_deg_9(coeffs::Vector{UInt}, a0::UInt, b0::UInt) =
    check_deg_9_300(coeffs, a0, b0)

# [b0,  1],
# [a0, b1,   1],
# [b0, a1,   0,   1],
# [a0, b2,  a2,   0,  1],
# [b0, a3,   0,   0,  0, 1],
# [a0, b3,  b4,  b5, a4, 0, 1],
# [b0, a5,   0,   0,  0, 0, 0, 1],
# [a0, b6,  b7,  b8, a6, 0, 0, 0, 1],
# [b0, b9, b10,  a7,  0, 0, 0, 0, 0, 1],

# Case 2: 300 < prime < 1000
# Case 3: 1000 < prime < 3370
# Case 4: 3370 < prime < 11000
# Case 5: 11000 < prime < 2^16
function check_coefficients(coeffs::Vector{Vector{UInt}}, case::Int, special::Bool = false)
    if case < 2 || case > 5
        error()
    end

    a0 = coeffs[2][1]
    b0 = coeffs[1][1]

    check_deg_1(coeffs[1], a0, b0)
    check_deg_2(coeffs[2], a0, b0)
    check_deg_3(coeffs[3], a0, b0)
    check_deg_4(coeffs[4], a0, b0)
    if case < 5
        check_deg_5(coeffs[5], a0, b0)
        check_deg_6(coeffs[6], a0, b0)
        if case < 4 && !special
            check_deg_7(coeffs[7], a0, b0)
            if case < 3
                check_deg_8(coeffs[8], a0, b0)
            end
            check_deg_9(coeffs[9], a0, b0)
        end
    end
end

function primes_lt_2p16(lines::Vector{String}, case::Int)
    if case < 2 || case > 5
        error()
    end

    specialprimes = [2689, 2797, 2833, 3019, 3163, 3209, 3331]

    if case == 2
        upperdegree = 9
        startingprime = 307
        upperprime = 1009
    elseif case == 3
        upperdegree = 9
        startingprime = 1009
        upperprime = 3371
    elseif case == 4
        upperdegree = 6
        startingprime = 3371
        upperprime = 11003
    elseif case == 5
        upperdegree = 4
        startingprime = 11003
        upperprime = 2^16
    end

    coeffs = [Vector{UInt}(undef, ix + 1) for ix in 1:upperdegree]

    # Current line
    lx = 1

    primes = "uint16_t __flint_cp_primes$case[] = {"
    sm_coeffs = "uint8_t __flint_cp_sm_coeffs$case[] = {"
    md_coeffs = "uint16_t __flint_cp_md_coeffs$case[] = {"

    # Get first prime
    ps, ix = get_comma_separated_item(lines[lx])
    prime = parse(UInt, ps)

    # Check that we start with the correct one
    ds, _ = get_comma_separated_item(lines[lx][ix:end])
    degree = parse(UInt, ds)
    if degree != 1 || prime != startingprime
        error()
    end

    while prime < upperprime
        primes *= "$(prime),"

        for degree in 1:upperdegree
            # Treat special cases
            if case == 3
                if degree == 8
                    continue
                elseif degree in 7:9 && prime in specialprimes
                    coeffs[degree] = zeros(UInt, degree + 1)
                    continue
                end
            end

            ps, ix = get_comma_separated_item(lines[lx])
            thisprime = parse(UInt, ps)
            lines[lx] = lines[lx][ix:end]

            if thisprime != prime
                error("lx = $lx, degree = $degree")
            end

            ds, ix = get_comma_separated_item(lines[lx])
            lines[lx] = lines[lx][ix:end]
            thisdegree = parse(UInt, ds)

            if thisdegree != degree
                error()
            end

            for jx in 0:degree
                cs, ix = get_comma_separated_item(lines[lx])
                lines[lx] = lines[lx][ix:end]
                coeffs[degree][jx + 1] = parse(UInt, cs)
            end

            lx += 1
        end

        check_coefficients(coeffs, case, prime in specialprimes)

        sm_coeffs *= "$(coeffs[2][1]),"
        sm_coeffs *= "$(coeffs[3][2]),"
        sm_coeffs *= "$(coeffs[4][3]),"
        if case < 5
            sm_coeffs *= "$(coeffs[5][2]),"
            sm_coeffs *= "$(coeffs[6][5]),"
            if case < 4
                sm_coeffs *= "$(coeffs[7][2]),"
                if case < 3
                    sm_coeffs *= "$(coeffs[8][5]),"
                end
                sm_coeffs *= "$(coeffs[9][4]),"
            end
        end

        md_coeffs *= "$(coeffs[1][1]),"
        md_coeffs *= "$(coeffs[2][2]),"
        md_coeffs *= "$(coeffs[4][2]),"
        if case < 5
            md_coeffs *= "$(coeffs[6][2]),$(coeffs[6][3]),$(coeffs[6][4]),"
            if case < 4
                if case < 3
                    md_coeffs *= "$(coeffs[8][2]),$(coeffs[8][3]),$(coeffs[8][4]),"
                end
                md_coeffs *= "$(coeffs[9][2]),$(coeffs[9][3]),"
            end
        end

        # Get new prime
        ps, _ = get_comma_separated_item(lines[lx])
        newprime = parse(UInt, ps)
        if newprime == prime
            error()
        end
        prime = newprime
    end

    primes = primes[1:end] * "0};"
    sm_coeffs = sm_coeffs[1:end - 1] * "};"
    md_coeffs = md_coeffs[1:end - 1] * "};"

    return lx - 1, primes, sm_coeffs, md_coeffs
end

###############################################################################
# 2^16 < prime < 110000
###############################################################################

function check_deg_4_gt_2p16(coeffs::Vector{UInt})
    if (coeffs[1] >= 2^8 || coeffs[2] >= 2^32 || coeffs[3] >= 2^8 ||
        coeffs[4] != 0 || coeffs[5] != 1)
        error("coeffs = $coeffs")
    end
end

function primes_gt_2p16(lines::Vector{String})
    coeffs = Vector{UInt}(undef, 5)

    # Current line
    lx = 1

    primes = "uint16_t __flint_cp_primes6[] = {"
    sm_coeffs = "uint8_t __flint_cp_sm_coeffs6[] = {"
    lg_coeffs = "uint32_t __flint_cp_lg_coeffs6[] = {"

    # Get first prime
    ps, ix = get_comma_separated_item(lines[lx])
    prime = parse(UInt, ps)

    # Check that we start with the correct one
    ds, _ = get_comma_separated_item(lines[lx][ix:end])
    degree = parse(UInt, ds)
    if degree != 4 || prime != 65537
        error()
    end

    while true
        primes *= "$(prime - 2^16),"

        ps, ix = get_comma_separated_item(lines[lx])
        thisprime = parse(UInt, ps)
        lines[lx] = lines[lx][ix:end]

        ds, ix = get_comma_separated_item(lines[lx])
        lines[lx] = lines[lx][ix:end]
        thisdegree = parse(UInt, ds)

        if thisdegree != 4
            error()
        end

        for jx in 0:degree
            cs, ix = get_comma_separated_item(lines[lx])
            lines[lx] = lines[lx][ix:end]
            coeffs[jx + 1] = parse(UInt, cs)
        end

        check_deg_4_gt_2p16(coeffs)

        sm_coeffs *= "$(coeffs[1]),"
        sm_coeffs *= "$(coeffs[3]),"
        lg_coeffs *= "$(coeffs[2]),"

        # If no more lines are to be found, break
        if length(lines[lx:end]) < 2
            break
        end

        lx += 1

        # Get new prime
        ps, _ = get_comma_separated_item(lines[lx])
        newprime = parse(UInt, ps)
        if newprime <= prime
            error()
        end
        prime = newprime
    end

    primes = primes[1:end] * "0};"
    sm_coeffs = sm_coeffs[1:end - 1] * "};"
    lg_coeffs = lg_coeffs[1:end - 1] * "};"

    return lx - 1, primes, sm_coeffs, lg_coeffs
end

###############################################################################
###############################################################################
# main
###############################################################################
###############################################################################

function main(filename::String)
    lines = read_lines()
    ix = 1

    lx, p0, d0, numntcoeffs0, ntcoeffs0 = primes_lt_260(lines[ix:end])
    ix += lx
    lx, p1, sm_coeffs1, md_coeffs1 = primes_lt_300(lines[ix:end])
    ix += lx
    lx, p2, sm_coeffs2, md_coeffs2 = primes_lt_2p16(lines[ix:end], 2)
    ix += lx
    lx, p3, sm_coeffs3, md_coeffs3 = primes_lt_2p16(lines[ix:end], 3)
    ix += lx
    lx, p4, sm_coeffs4, md_coeffs4 = primes_lt_2p16(lines[ix:end], 4)
    ix += lx
    lx, p5, sm_coeffs5, md_coeffs5 = primes_lt_2p16(lines[ix:end], 5)
    ix += lx
    lx, p6, sm_coeffs6, lg_coeffs6 = primes_gt_2p16(lines[ix:end])

    file = open(filename, "w")
    write(file, p0)
    write(file, "\n")
    write(file, d0)
    write(file, "\n")
    write(file, numntcoeffs0)
    write(file, "\n")
    write(file, ntcoeffs0)
    write(file, "\n\n")
    write(file, p1)
    write(file, "\n")
    write(file, sm_coeffs1)
    write(file, "\n")
    write(file, md_coeffs1)
    write(file, "\n\n")
    write(file, p2)
    write(file, "\n")
    write(file, sm_coeffs2)
    write(file, "\n")
    write(file, md_coeffs2)
    write(file, "\n\n")
    write(file, p3)
    write(file, "\n")
    write(file, sm_coeffs3)
    write(file, "\n")
    write(file, md_coeffs3)
    write(file, "\n\n")
    write(file, p4)
    write(file, "\n")
    write(file, sm_coeffs4)
    write(file, "\n")
    write(file, md_coeffs4)
    write(file, "\n\n")
    write(file, p5)
    write(file, "\n")
    write(file, sm_coeffs5)
    write(file, "\n")
    write(file, md_coeffs5)
    write(file, "\n\n")
    write(file, p6)
    write(file, "\n")
    write(file, sm_coeffs6)
    write(file, "\n")
    write(file, lg_coeffs6)
    write(file, "\n")
    close(file)
end
