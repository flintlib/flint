/*
Evaluate two algebraic numbers defined by huge expressions
and verify that they are equal.

Example from https://ask.sagemath.org/question/52653/
    equality-of-algebraic-numbers-given-by-huge-symbolic-expressions/

This file is public domain. Author: Fredrik Johansson.
*/

#include <string.h>
#include <flint/profiler.h>
#include <flint/calcium.h>
#include <flint/ca.h>
#include <flint/gr.h>

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif

const char * EXPR_N = 
"1/16*(44*(7*sqrt(2) - 10)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) + 2*(11*"
"(7*sqrt(2) - 10)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - 10*(63*sqrt(2) "
"- 89)*sqrt(sqrt(2) + 2) - (3*(3*sqrt(2) - 4)*sqrt(sqrt(2) + 2)*sqrt(-17*sqr"
"t(2) + 26) - (85*sqrt(2) - 122)*sqrt(sqrt(2) + 2))*sqrt(-sqrt(2) + 2) + 2*("
"(3*(3*sqrt(2) - 4)*sqrt(-17*sqrt(2) + 26) - 85*sqrt(2) + 122)*sqrt(-sqrt(2)"
" + 2) - 11*(7*sqrt(2) - 10)*sqrt(-17*sqrt(2) + 26) + 630*sqrt(2) - 890)*sqr"
"t(sqrt(sqrt(2) + 2) - 1))*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3) - 40"
"*(63*sqrt(2) - 89)*sqrt(sqrt(2) + 2) - 4*(3*(3*sqrt(2) - 4)*sqrt(sqrt(2) + "
"2)*sqrt(-17*sqrt(2) + 26) - (85*sqrt(2) - 122)*sqrt(sqrt(2) + 2))*sqrt(-sqr"
"t(2) + 2) + (22*(5*sqrt(2) - 7)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) + "
"(11*(5*sqrt(2) - 7)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - 5*(89*sqrt(2"
") - 126)*sqrt(sqrt(2) + 2) - (3*(2*sqrt(2) - 3)*sqrt(sqrt(2) + 2)*sqrt(-17*"
"sqrt(2) + 26) - (61*sqrt(2) - 85)*sqrt(sqrt(2) + 2))*sqrt(-sqrt(2) + 2) + 2"
"*((3*(2*sqrt(2) - 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt(2"
") + 2) - 11*(5*sqrt(2) - 7)*sqrt(-17*sqrt(2) + 26) + 445*sqrt(2) - 630)*sqr"
"t(sqrt(sqrt(2) + 2) - 1))*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3) - 10"
"*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2) - 2*(3*(2*sqrt(2) - 3)*sqrt(sqrt(2) +"
" 2)*sqrt(-17*sqrt(2) + 26) - (61*sqrt(2) - 85)*sqrt(sqrt(2) + 2))*sqrt(-sqr"
"t(2) + 2) + 4*((3*(2*sqrt(2) - 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)"
"*sqrt(-sqrt(2) + 2) - 11*(5*sqrt(2) - 7)*sqrt(-17*sqrt(2) + 26) + 445*sqrt("
"2) - 630)*sqrt(sqrt(sqrt(2) + 2) - 1))*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt(2) +"
" 2) - 2*sqrt(-17*sqrt(2) + 26) + 24) + 8*((3*(3*sqrt(2) - 4)*sqrt(-17*sqrt("
"2) + 26) - 85*sqrt(2) + 122)*sqrt(-sqrt(2) + 2) - 11*(7*sqrt(2) - 10)*sqrt("
"-17*sqrt(2) + 26) + 630*sqrt(2) - 890)*sqrt(sqrt(sqrt(2) + 2) - 1))*(((sqrt"
"(2)*sqrt(sqrt(2) + 2) - sqrt(2) - 1)*sqrt(sqrt(sqrt(2) + 2) - 1) + sqrt(sqr"
"t(2) + 2) - 1)*(8*((5*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2) - ((61*sqrt(2) -"
" 85)*sqrt(sqrt(2) + 2) - 3*((2*sqrt(2) - 3)*sqrt(sqrt(2) + 2) - 4*sqrt(2) +"
" 6)*sqrt(-17*sqrt(2) + 26) - 122*sqrt(2) + 170)*sqrt(-sqrt(2) + 2) - 11*((5"
"*sqrt(2) - 7)*sqrt(sqrt(2) + 2) - 10*sqrt(2) + 14)*sqrt(-17*sqrt(2) + 26) -"
" 890*sqrt(2) + 1260)*sqrt(3*sqrt(2) + sqrt(-sqrt(2) + 2) - 5)*sqrt(-12*sqrt"
"(2) - 2*sqrt(-sqrt(2) + 2) - 2*sqrt(-17*sqrt(2) + 26) + 24) + 2*(10*(63*sqr"
"t(2) - 89)*sqrt(sqrt(2) + 2) - ((85*sqrt(2) - 122)*sqrt(sqrt(2) + 2) - 3*(("
"3*sqrt(2) - 4)*sqrt(sqrt(2) + 2) - 6*sqrt(2) + 8)*sqrt(-17*sqrt(2) + 26) - "
"170*sqrt(2) + 244)*sqrt(-sqrt(2) + 2) - 11*((7*sqrt(2) - 10)*sqrt(sqrt(2) +"
" 2) - 14*sqrt(2) + 20)*sqrt(-17*sqrt(2) + 26) - 1260*sqrt(2) + 1780)*sqrt(3"
"*sqrt(2) + sqrt(-sqrt(2) + 2) - 5))*(1/((sqrt(sqrt(2) + 2) + 1)*sqrt(sqrt(s"
"qrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1) + ((sqrt(2) + sqrt(sqrt(2) + 2))*"
"sqrt(sqrt(sqrt(2) + 2) - 1) + 3*sqrt(2) - 5*sqrt(sqrt(2) + 2) + 8)*((sqrt(2"
") + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2) - 1) - 3*sqrt(2) + 5*sqrt(sqr"
"t(2) + 2) - 8)*(sqrt(sqrt(2) + 2) + 2)/(((sqrt(sqrt(2) + 2) + 1)*sqrt(sqrt("
"sqrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1)^2*((sqrt(sqrt(2) + 2) + 1)*sqrt("
"sqrt(sqrt(2) + 2) - 1) - sqrt(sqrt(2) + 2) - ((sqrt(2) + sqrt(sqrt(2) + 2))"
"*sqrt(sqrt(sqrt(2) + 2) - 1) + 3*sqrt(2) - 5*sqrt(sqrt(2) + 2) + 8)*((sqrt("
"2) + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2) - 1) - 3*sqrt(2) + 5*sqrt(sq"
"rt(2) + 2) - 8)*(sqrt(sqrt(2) + 2) + 2)/(((sqrt(sqrt(2) + 2) + 1)*sqrt(sqrt"
"(sqrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1)*(sqrt(sqrt(2) + 2) - 2)^3) - 1)"
"*(sqrt(sqrt(2) + 2) - 2)^3))/(2*((3*(3*sqrt(2) - 4)*sqrt(-17*sqrt(2) + 26) "
"- 85*sqrt(2) + 122)*sqrt(-sqrt(2) + 2) - 11*(7*sqrt(2) - 10)*sqrt(-17*sqrt("
"2) + 26) + 630*sqrt(2) - 890)*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3)*"
"sqrt(sqrt(sqrt(2) + 2) - 1) - 2*(4896*sqrt(2) - 6923)*sqrt(sqrt(2) + 2) + ("
"20*(79*sqrt(2) - 112)*sqrt(sqrt(2) + 2) - (7*(27*sqrt(2) - 38)*sqrt(sqrt(2)"
" + 2) - 342*sqrt(2) + 484)*sqrt(-17*sqrt(2) + 26) - 2820*sqrt(2) + 3992)*sq"
"rt(-sqrt(2) + 2) + (((3*(2*sqrt(2) - 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2)"
" + 85)*sqrt(-sqrt(2) + 2) - 11*(5*sqrt(2) - 7)*sqrt(-17*sqrt(2) + 26) + 445"
"*sqrt(2) - 630)*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3)*sqrt(sqrt(sqrt"
"(2) + 2) - 1) - 10*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2) + 2*((61*sqrt(2) - "
"85)*sqrt(sqrt(2) + 2) - 3*((2*sqrt(2) - 3)*sqrt(sqrt(2) + 2) - 2*sqrt(2) + "
"3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt(2) + 2) + 22*((5*sq"
"rt(2) - 7)*sqrt(sqrt(2) + 2) - 5*sqrt(2) + 7)*sqrt(-17*sqrt(2) + 26) + 890*"
"sqrt(2) - 1260)*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt(2) + 2) - 2*sqrt(-17*sqrt(2"
") + 26) + 24) + 4*((319*sqrt(2) - 452)*sqrt(sqrt(2) + 2) - 561*sqrt(2) + 79"
"4)*sqrt(-17*sqrt(2) + 26) + 17064*sqrt(2) - 24132) + (44*(7*sqrt(2) - 10)*s"
"qrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) + 2*(11*(7*sqrt(2) - 10)*sqrt(sqrt("
"2) + 2)*sqrt(-17*sqrt(2) + 26) - 10*(63*sqrt(2) - 89)*sqrt(sqrt(2) + 2) - ("
"3*(3*sqrt(2) - 4)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - (85*sqrt(2) - "
"122)*sqrt(sqrt(2) + 2))*sqrt(-sqrt(2) + 2) + 2*((3*(3*sqrt(2) - 4)*sqrt(-17"
"*sqrt(2) + 26) - 85*sqrt(2) + 122)*sqrt(-sqrt(2) + 2) - 11*(7*sqrt(2) - 10)"
"*sqrt(-17*sqrt(2) + 26) + 630*sqrt(2) - 890)*sqrt(sqrt(sqrt(2) + 2) - 1))*s"
"qrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3) - 40*(63*sqrt(2) - 89)*sqrt(sqr"
"t(2) + 2) - 4*(3*(3*sqrt(2) - 4)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) -"
" (85*sqrt(2) - 122)*sqrt(sqrt(2) + 2))*sqrt(-sqrt(2) + 2) + (22*(5*sqrt(2) "
"- 7)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) + (11*(5*sqrt(2) - 7)*sqrt(sq"
"rt(2) + 2)*sqrt(-17*sqrt(2) + 26) - 5*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2) "
"- (3*(2*sqrt(2) - 3)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - (61*sqrt(2)"
" - 85)*sqrt(sqrt(2) + 2))*sqrt(-sqrt(2) + 2) + 2*((3*(2*sqrt(2) - 3)*sqrt(-"
"17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt(2) + 2) - 11*(5*sqrt(2) - 7)"
"*sqrt(-17*sqrt(2) + 26) + 445*sqrt(2) - 630)*sqrt(sqrt(sqrt(2) + 2) - 1))*s"
"qrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3) - 10*(89*sqrt(2) - 126)*sqrt(sq"
"rt(2) + 2) - 2*(3*(2*sqrt(2) - 3)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) "
"- (61*sqrt(2) - 85)*sqrt(sqrt(2) + 2))*sqrt(-sqrt(2) + 2) + 4*((3*(2*sqrt(2"
") - 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt(2) + 2) - 11*(5"
"*sqrt(2) - 7)*sqrt(-17*sqrt(2) + 26) + 445*sqrt(2) - 630)*sqrt(sqrt(sqrt(2)"
" + 2) - 1))*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt(2) + 2) - 2*sqrt(-17*sqrt(2) + "
"26) + 24) + 8*((3*(3*sqrt(2) - 4)*sqrt(-17*sqrt(2) + 26) - 85*sqrt(2) + 122"
")*sqrt(-sqrt(2) + 2) - 11*(7*sqrt(2) - 10)*sqrt(-17*sqrt(2) + 26) + 630*sqr"
"t(2) - 890)*sqrt(sqrt(sqrt(2) + 2) - 1))*((sqrt(2) + sqrt(sqrt(2) + 2))*sqr"
"t(sqrt(sqrt(2) + 2) - 1) + 3*sqrt(2) - 5*sqrt(sqrt(2) + 2) + 8)*sqrt(sqrt(s"
"qrt(2) + 2) + 2)/((2*((3*(3*sqrt(2) - 4)*sqrt(-17*sqrt(2) + 26) - 85*sqrt(2"
") + 122)*sqrt(-sqrt(2) + 2) - 11*(7*sqrt(2) - 10)*sqrt(-17*sqrt(2) + 26) + "
"630*sqrt(2) - 890)*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3)*sqrt(sqrt(s"
"qrt(2) + 2) - 1) - 2*(4896*sqrt(2) - 6923)*sqrt(sqrt(2) + 2) + (20*(79*sqrt"
"(2) - 112)*sqrt(sqrt(2) + 2) - (7*(27*sqrt(2) - 38)*sqrt(sqrt(2) + 2) - 342"
"*sqrt(2) + 484)*sqrt(-17*sqrt(2) + 26) - 2820*sqrt(2) + 3992)*sqrt(-sqrt(2)"
" + 2) + (((3*(2*sqrt(2) - 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt"
"(-sqrt(2) + 2) - 11*(5*sqrt(2) - 7)*sqrt(-17*sqrt(2) + 26) + 445*sqrt(2) - "
"630)*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3)*sqrt(sqrt(sqrt(2) + 2) - "
"1) - 10*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2) + 2*((61*sqrt(2) - 85)*sqrt(sq"
"rt(2) + 2) - 3*((2*sqrt(2) - 3)*sqrt(sqrt(2) + 2) - 2*sqrt(2) + 3)*sqrt(-17"
"*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt(2) + 2) + 22*((5*sqrt(2) - 7)*"
"sqrt(sqrt(2) + 2) - 5*sqrt(2) + 7)*sqrt(-17*sqrt(2) + 26) + 890*sqrt(2) - 1"
"260)*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt(2) + 2) - 2*sqrt(-17*sqrt(2) + 26) + 2"
"4) + 4*((319*sqrt(2) - 452)*sqrt(sqrt(2) + 2) - 561*sqrt(2) + 794)*sqrt(-17"
"*sqrt(2) + 26) + 17064*sqrt(2) - 24132)*((sqrt(sqrt(2) + 2) + 1)*sqrt(sqrt("
"sqrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1)*((sqrt(sqrt(2) + 2) + 1)*sqrt(sq"
"rt(sqrt(2) + 2) - 1) - sqrt(sqrt(2) + 2) - ((sqrt(2) + sqrt(sqrt(2) + 2))*s"
"qrt(sqrt(sqrt(2) + 2) - 1) + 3*sqrt(2) - 5*sqrt(sqrt(2) + 2) + 8)*((sqrt(2)"
" + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2) - 1) - 3*sqrt(2) + 5*sqrt(sqrt"
"(2) + 2) - 8)*(sqrt(sqrt(2) + 2) + 2)/(((sqrt(sqrt(2) + 2) + 1)*sqrt(sqrt(s"
"qrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1)*(sqrt(sqrt(2) + 2) - 2)^3) - 1)*("
"-1/4*sqrt(sqrt(2) + 2) + 1/2)^(3/2)))*sqrt(sqrt(sqrt(2) + 2) + 2)/sqrt(-1/4"
"*sqrt(sqrt(2) + 2) + 1/2) + 2*(sqrt(sqrt(2) + 2)*(sqrt(2) - 1)*sqrt(sqrt(sq"
"rt(2) + 2) - 1) - sqrt(2) + sqrt(sqrt(2) + 2))*(8*(44*(7*sqrt(2) - 10)*sqrt"
"(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - 2*(11*(7*sqrt(2) - 10)*sqrt(sqrt(2) "
"+ 2)*sqrt(-17*sqrt(2) + 26) - 10*(63*sqrt(2) - 89)*sqrt(sqrt(2) + 2) - (3*("
"3*sqrt(2) - 4)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - (85*sqrt(2) - 122"
")*sqrt(sqrt(2) + 2))*sqrt(-sqrt(2) + 2) - 2*((3*(3*sqrt(2) - 4)*sqrt(-17*sq"
"rt(2) + 26) - 85*sqrt(2) + 122)*sqrt(-sqrt(2) + 2) - 11*(7*sqrt(2) - 10)*sq"
"rt(-17*sqrt(2) + 26) + 630*sqrt(2) - 890)*sqrt(sqrt(sqrt(2) + 2) - 1))*sqrt"
"(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3) - 40*(63*sqrt(2) - 89)*sqrt(sqrt(2"
") + 2) - 4*(3*(3*sqrt(2) - 4)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - (8"
"5*sqrt(2) - 122)*sqrt(sqrt(2) + 2))*sqrt(-sqrt(2) + 2) + (22*(5*sqrt(2) - 7"
")*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - (11*(5*sqrt(2) - 7)*sqrt(sqrt("
"2) + 2)*sqrt(-17*sqrt(2) + 26) - 5*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2) - ("
"3*(2*sqrt(2) - 3)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - (61*sqrt(2) - "
"85)*sqrt(sqrt(2) + 2))*sqrt(-sqrt(2) + 2) - 2*((3*(2*sqrt(2) - 3)*sqrt(-17*"
"sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt(2) + 2) - 11*(5*sqrt(2) - 7)*sq"
"rt(-17*sqrt(2) + 26) + 445*sqrt(2) - 630)*sqrt(sqrt(sqrt(2) + 2) - 1))*sqrt"
"(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3) - 10*(89*sqrt(2) - 126)*sqrt(sqrt("
"2) + 2) - 2*(3*(2*sqrt(2) - 3)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - ("
"61*sqrt(2) - 85)*sqrt(sqrt(2) + 2))*sqrt(-sqrt(2) + 2) - 4*((3*(2*sqrt(2) -"
" 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt(2) + 2) - 11*(5*sq"
"rt(2) - 7)*sqrt(-17*sqrt(2) + 26) + 445*sqrt(2) - 630)*sqrt(sqrt(sqrt(2) + "
"2) - 1))*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt(2) + 2) - 2*sqrt(-17*sqrt(2) + 26)"
" + 24) - 8*((3*(3*sqrt(2) - 4)*sqrt(-17*sqrt(2) + 26) - 85*sqrt(2) + 122)*s"
"qrt(-sqrt(2) + 2) - 11*(7*sqrt(2) - 10)*sqrt(-17*sqrt(2) + 26) + 630*sqrt(2"
") - 890)*sqrt(sqrt(sqrt(2) + 2) - 1))*(1/((sqrt(sqrt(2) + 2) + 1)*sqrt(sqrt"
"(sqrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1) + ((sqrt(2) + sqrt(sqrt(2) + 2)"
")*sqrt(sqrt(sqrt(2) + 2) - 1) + 3*sqrt(2) - 5*sqrt(sqrt(2) + 2) + 8)*((sqrt"
"(2) + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2) - 1) - 3*sqrt(2) + 5*sqrt(s"
"qrt(2) + 2) - 8)*(sqrt(sqrt(2) + 2) + 2)/(((sqrt(sqrt(2) + 2) + 1)*sqrt(sqr"
"t(sqrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1)^2*((sqrt(sqrt(2) + 2) + 1)*sqr"
"t(sqrt(sqrt(2) + 2) - 1) - sqrt(sqrt(2) + 2) - ((sqrt(2) + sqrt(sqrt(2) + 2"
"))*sqrt(sqrt(sqrt(2) + 2) - 1) + 3*sqrt(2) - 5*sqrt(sqrt(2) + 2) + 8)*((sqr"
"t(2) + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2) - 1) - 3*sqrt(2) + 5*sqrt("
"sqrt(2) + 2) - 8)*(sqrt(sqrt(2) + 2) + 2)/(((sqrt(sqrt(2) + 2) + 1)*sqrt(sq"
"rt(sqrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1)*(sqrt(sqrt(2) + 2) - 2)^3) - "
"1)*(sqrt(sqrt(2) + 2) - 2)^3))/(2*((3*(3*sqrt(2) - 4)*sqrt(-17*sqrt(2) + 26"
") - 85*sqrt(2) + 122)*sqrt(-sqrt(2) + 2) - 11*(7*sqrt(2) - 10)*sqrt(-17*sqr"
"t(2) + 26) + 630*sqrt(2) - 890)*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3"
")*sqrt(sqrt(sqrt(2) + 2) - 1) - 2*(4896*sqrt(2) - 6923)*sqrt(sqrt(2) + 2) +"
" (20*(79*sqrt(2) - 112)*sqrt(sqrt(2) + 2) - (7*(27*sqrt(2) - 38)*sqrt(sqrt("
"2) + 2) - 342*sqrt(2) + 484)*sqrt(-17*sqrt(2) + 26) - 2820*sqrt(2) + 3992)*"
"sqrt(-sqrt(2) + 2) + (((3*(2*sqrt(2) - 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt("
"2) + 85)*sqrt(-sqrt(2) + 2) - 11*(5*sqrt(2) - 7)*sqrt(-17*sqrt(2) + 26) + 4"
"45*sqrt(2) - 630)*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3)*sqrt(sqrt(sq"
"rt(2) + 2) - 1) - 10*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2) + 2*((61*sqrt(2) "
"- 85)*sqrt(sqrt(2) + 2) - 3*((2*sqrt(2) - 3)*sqrt(sqrt(2) + 2) - 2*sqrt(2) "
"+ 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt(2) + 2) + 22*((5*"
"sqrt(2) - 7)*sqrt(sqrt(2) + 2) - 5*sqrt(2) + 7)*sqrt(-17*sqrt(2) + 26) + 89"
"0*sqrt(2) - 1260)*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt(2) + 2) - 2*sqrt(-17*sqrt"
"(2) + 26) + 24) + 4*((319*sqrt(2) - 452)*sqrt(sqrt(2) + 2) - 561*sqrt(2) + "
"794)*sqrt(-17*sqrt(2) + 26) + 17064*sqrt(2) - 24132) + ((5*(89*sqrt(2) - 12"
"6)*sqrt(sqrt(2) + 2) - ((61*sqrt(2) - 85)*sqrt(sqrt(2) + 2) - 3*((2*sqrt(2)"
" - 3)*sqrt(sqrt(2) + 2) - 4*sqrt(2) + 6)*sqrt(-17*sqrt(2) + 26) - 122*sqrt("
"2) + 170)*sqrt(-sqrt(2) + 2) - 11*((5*sqrt(2) - 7)*sqrt(sqrt(2) + 2) - 10*s"
"qrt(2) + 14)*sqrt(-17*sqrt(2) + 26) - 890*sqrt(2) + 1260)*sqrt(3*sqrt(2) + "
"sqrt(-sqrt(2) + 2) - 5)*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt(2) + 2) - 2*sqrt(-1"
"7*sqrt(2) + 26) + 24) + 2*(10*(63*sqrt(2) - 89)*sqrt(sqrt(2) + 2) - ((85*sq"
"rt(2) - 122)*sqrt(sqrt(2) + 2) - 3*((3*sqrt(2) - 4)*sqrt(sqrt(2) + 2) - 6*s"
"qrt(2) + 8)*sqrt(-17*sqrt(2) + 26) - 170*sqrt(2) + 244)*sqrt(-sqrt(2) + 2) "
"- 11*((7*sqrt(2) - 10)*sqrt(sqrt(2) + 2) - 14*sqrt(2) + 20)*sqrt(-17*sqrt(2"
") + 26) - 1260*sqrt(2) + 1780)*sqrt(3*sqrt(2) + sqrt(-sqrt(2) + 2) - 5))*(("
"sqrt(2) + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2) - 1) + 3*sqrt(2) - 5*sq"
"rt(sqrt(2) + 2) + 8)*sqrt(sqrt(sqrt(2) + 2) + 2)/((2*((3*(3*sqrt(2) - 4)*sq"
"rt(-17*sqrt(2) + 26) - 85*sqrt(2) + 122)*sqrt(-sqrt(2) + 2) - 11*(7*sqrt(2)"
" - 10)*sqrt(-17*sqrt(2) + 26) + 630*sqrt(2) - 890)*sqrt(3*sqrt(2) + sqrt(-1"
"7*sqrt(2) + 26) - 3)*sqrt(sqrt(sqrt(2) + 2) - 1) - 2*(4896*sqrt(2) - 6923)*"
"sqrt(sqrt(2) + 2) + (20*(79*sqrt(2) - 112)*sqrt(sqrt(2) + 2) - (7*(27*sqrt("
"2) - 38)*sqrt(sqrt(2) + 2) - 342*sqrt(2) + 484)*sqrt(-17*sqrt(2) + 26) - 28"
"20*sqrt(2) + 3992)*sqrt(-sqrt(2) + 2) + (((3*(2*sqrt(2) - 3)*sqrt(-17*sqrt("
"2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt(2) + 2) - 11*(5*sqrt(2) - 7)*sqrt(-1"
"7*sqrt(2) + 26) + 445*sqrt(2) - 630)*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26"
") - 3)*sqrt(sqrt(sqrt(2) + 2) - 1) - 10*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2"
") + 2*((61*sqrt(2) - 85)*sqrt(sqrt(2) + 2) - 3*((2*sqrt(2) - 3)*sqrt(sqrt(2"
") + 2) - 2*sqrt(2) + 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqr"
"t(2) + 2) + 22*((5*sqrt(2) - 7)*sqrt(sqrt(2) + 2) - 5*sqrt(2) + 7)*sqrt(-17"
"*sqrt(2) + 26) + 890*sqrt(2) - 1260)*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt(2) + 2"
") - 2*sqrt(-17*sqrt(2) + 26) + 24) + 4*((319*sqrt(2) - 452)*sqrt(sqrt(2) + "
"2) - 561*sqrt(2) + 794)*sqrt(-17*sqrt(2) + 26) + 17064*sqrt(2) - 24132)*((s"
"qrt(sqrt(2) + 2) + 1)*sqrt(sqrt(sqrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1)*"
"((sqrt(sqrt(2) + 2) + 1)*sqrt(sqrt(sqrt(2) + 2) - 1) - sqrt(sqrt(2) + 2) - "
"((sqrt(2) + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2) - 1) + 3*sqrt(2) - 5*"
"sqrt(sqrt(2) + 2) + 8)*((sqrt(2) + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2"
") - 1) - 3*sqrt(2) + 5*sqrt(sqrt(2) + 2) - 8)*(sqrt(sqrt(2) + 2) + 2)/(((sq"
"rt(sqrt(2) + 2) + 1)*sqrt(sqrt(sqrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1)*("
"sqrt(sqrt(2) + 2) - 2)^3) - 1)*(-1/4*sqrt(sqrt(2) + 2) + 1/2)^(3/2)))/(sqrt"
"(sqrt(2) + 2) - 2))/(2*((3*(3*sqrt(2) - 4)*sqrt(-17*sqrt(2) + 26) - 85*sqrt"
"(2) + 122)*sqrt(-sqrt(2) + 2) - 11*(7*sqrt(2) - 10)*sqrt(-17*sqrt(2) + 26) "
"+ 630*sqrt(2) - 890)*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3)*sqrt(sqrt"
"(sqrt(2) + 2) - 1) - 2*(4896*sqrt(2) - 6923)*sqrt(sqrt(2) + 2) + (20*(79*sq"
"rt(2) - 112)*sqrt(sqrt(2) + 2) - (7*(27*sqrt(2) - 38)*sqrt(sqrt(2) + 2) - 3"
"42*sqrt(2) + 484)*sqrt(-17*sqrt(2) + 26) - 2820*sqrt(2) + 3992)*sqrt(-sqrt("
"2) + 2) + (((3*(2*sqrt(2) - 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sq"
"rt(-sqrt(2) + 2) - 11*(5*sqrt(2) - 7)*sqrt(-17*sqrt(2) + 26) + 445*sqrt(2) "
"- 630)*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3)*sqrt(sqrt(sqrt(2) + 2) "
"- 1) - 10*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2) + 2*((61*sqrt(2) - 85)*sqrt("
"sqrt(2) + 2) - 3*((2*sqrt(2) - 3)*sqrt(sqrt(2) + 2) - 2*sqrt(2) + 3)*sqrt(-"
"17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt(2) + 2) + 22*((5*sqrt(2) - 7"
")*sqrt(sqrt(2) + 2) - 5*sqrt(2) + 7)*sqrt(-17*sqrt(2) + 26) + 890*sqrt(2) -"
" 1260)*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt(2) + 2) - 2*sqrt(-17*sqrt(2) + 26) +"
" 24) + 4*((319*sqrt(2) - 452)*sqrt(sqrt(2) + 2) - 561*sqrt(2) + 794)*sqrt(-"
"17*sqrt(2) + 26) + 17064*sqrt(2) - 24132) - 1/16*((5*(89*sqrt(2) - 126)*sqr"
"t(sqrt(2) + 2) - ((61*sqrt(2) - 85)*sqrt(sqrt(2) + 2) - 3*((2*sqrt(2) - 3)*"
"sqrt(sqrt(2) + 2) - 4*sqrt(2) + 6)*sqrt(-17*sqrt(2) + 26) - 122*sqrt(2) + 1"
"70)*sqrt(-sqrt(2) + 2) - 11*((5*sqrt(2) - 7)*sqrt(sqrt(2) + 2) - 10*sqrt(2)"
" + 14)*sqrt(-17*sqrt(2) + 26) - 890*sqrt(2) + 1260)*sqrt(3*sqrt(2) + sqrt(-"
"sqrt(2) + 2) - 5)*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt(2) + 2) - 2*sqrt(-17*sqrt"
"(2) + 26) + 24) + 2*(10*(63*sqrt(2) - 89)*sqrt(sqrt(2) + 2) - ((85*sqrt(2) "
"- 122)*sqrt(sqrt(2) + 2) - 3*((3*sqrt(2) - 4)*sqrt(sqrt(2) + 2) - 6*sqrt(2)"
" + 8)*sqrt(-17*sqrt(2) + 26) - 170*sqrt(2) + 244)*sqrt(-sqrt(2) + 2) - 11*("
"(7*sqrt(2) - 10)*sqrt(sqrt(2) + 2) - 14*sqrt(2) + 20)*sqrt(-17*sqrt(2) + 26"
") - 1260*sqrt(2) + 1780)*sqrt(3*sqrt(2) + sqrt(-sqrt(2) + 2) - 5))*(((sqrt("
"2)*sqrt(sqrt(2) + 2) - sqrt(2) - 1)*sqrt(sqrt(sqrt(2) + 2) - 1) - sqrt(sqrt"
"(2) + 2) + 1)*(8*(44*(7*sqrt(2) - 10)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + "
"26) - 2*(11*(7*sqrt(2) - 10)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - 10*"
"(63*sqrt(2) - 89)*sqrt(sqrt(2) + 2) - (3*(3*sqrt(2) - 4)*sqrt(sqrt(2) + 2)*"
"sqrt(-17*sqrt(2) + 26) - (85*sqrt(2) - 122)*sqrt(sqrt(2) + 2))*sqrt(-sqrt(2"
") + 2) - 2*((3*(3*sqrt(2) - 4)*sqrt(-17*sqrt(2) + 26) - 85*sqrt(2) + 122)*s"
"qrt(-sqrt(2) + 2) - 11*(7*sqrt(2) - 10)*sqrt(-17*sqrt(2) + 26) + 630*sqrt(2"
") - 890)*sqrt(sqrt(sqrt(2) + 2) - 1))*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 2"
"6) - 3) - 40*(63*sqrt(2) - 89)*sqrt(sqrt(2) + 2) - 4*(3*(3*sqrt(2) - 4)*sqr"
"t(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - (85*sqrt(2) - 122)*sqrt(sqrt(2) + 2"
"))*sqrt(-sqrt(2) + 2) + (22*(5*sqrt(2) - 7)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt"
"(2) + 26) - (11*(5*sqrt(2) - 7)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - "
"5*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2) - (3*(2*sqrt(2) - 3)*sqrt(sqrt(2) + "
"2)*sqrt(-17*sqrt(2) + 26) - (61*sqrt(2) - 85)*sqrt(sqrt(2) + 2))*sqrt(-sqrt"
"(2) + 2) - 2*((3*(2*sqrt(2) - 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)*"
"sqrt(-sqrt(2) + 2) - 11*(5*sqrt(2) - 7)*sqrt(-17*sqrt(2) + 26) + 445*sqrt(2"
") - 630)*sqrt(sqrt(sqrt(2) + 2) - 1))*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 2"
"6) - 3) - 10*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2) - 2*(3*(2*sqrt(2) - 3)*sq"
"rt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - (61*sqrt(2) - 85)*sqrt(sqrt(2) + 2"
"))*sqrt(-sqrt(2) + 2) - 4*((3*(2*sqrt(2) - 3)*sqrt(-17*sqrt(2) + 26) - 61*s"
"qrt(2) + 85)*sqrt(-sqrt(2) + 2) - 11*(5*sqrt(2) - 7)*sqrt(-17*sqrt(2) + 26)"
" + 445*sqrt(2) - 630)*sqrt(sqrt(sqrt(2) + 2) - 1))*sqrt(-12*sqrt(2) - 2*sqr"
"t(-sqrt(2) + 2) - 2*sqrt(-17*sqrt(2) + 26) + 24) - 8*((3*(3*sqrt(2) - 4)*sq"
"rt(-17*sqrt(2) + 26) - 85*sqrt(2) + 122)*sqrt(-sqrt(2) + 2) - 11*(7*sqrt(2)"
" - 10)*sqrt(-17*sqrt(2) + 26) + 630*sqrt(2) - 890)*sqrt(sqrt(sqrt(2) + 2) -"
" 1))*(1/((sqrt(sqrt(2) + 2) + 1)*sqrt(sqrt(sqrt(2) + 2) - 1) + sqrt(sqrt(2)"
" + 2) + 1) + ((sqrt(2) + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2) - 1) + 3"
"*sqrt(2) - 5*sqrt(sqrt(2) + 2) + 8)*((sqrt(2) + sqrt(sqrt(2) + 2))*sqrt(sqr"
"t(sqrt(2) + 2) - 1) - 3*sqrt(2) + 5*sqrt(sqrt(2) + 2) - 8)*(sqrt(sqrt(2) + "
"2) + 2)/(((sqrt(sqrt(2) + 2) + 1)*sqrt(sqrt(sqrt(2) + 2) - 1) + sqrt(sqrt(2"
") + 2) + 1)^2*((sqrt(sqrt(2) + 2) + 1)*sqrt(sqrt(sqrt(2) + 2) - 1) - sqrt(s"
"qrt(2) + 2) - ((sqrt(2) + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2) - 1) + "
"3*sqrt(2) - 5*sqrt(sqrt(2) + 2) + 8)*((sqrt(2) + sqrt(sqrt(2) + 2))*sqrt(sq"
"rt(sqrt(2) + 2) - 1) - 3*sqrt(2) + 5*sqrt(sqrt(2) + 2) - 8)*(sqrt(sqrt(2) +"
" 2) + 2)/(((sqrt(sqrt(2) + 2) + 1)*sqrt(sqrt(sqrt(2) + 2) - 1) + sqrt(sqrt("
"2) + 2) + 1)*(sqrt(sqrt(2) + 2) - 2)^3) - 1)*(sqrt(sqrt(2) + 2) - 2)^3))/(2"
"*((3*(3*sqrt(2) - 4)*sqrt(-17*sqrt(2) + 26) - 85*sqrt(2) + 122)*sqrt(-sqrt("
"2) + 2) - 11*(7*sqrt(2) - 10)*sqrt(-17*sqrt(2) + 26) + 630*sqrt(2) - 890)*s"
"qrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3)*sqrt(sqrt(sqrt(2) + 2) - 1) - 2"
"*(4896*sqrt(2) - 6923)*sqrt(sqrt(2) + 2) + (20*(79*sqrt(2) - 112)*sqrt(sqrt"
"(2) + 2) - (7*(27*sqrt(2) - 38)*sqrt(sqrt(2) + 2) - 342*sqrt(2) + 484)*sqrt"
"(-17*sqrt(2) + 26) - 2820*sqrt(2) + 3992)*sqrt(-sqrt(2) + 2) + (((3*(2*sqrt"
"(2) - 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt(2) + 2) - 11*"
"(5*sqrt(2) - 7)*sqrt(-17*sqrt(2) + 26) + 445*sqrt(2) - 630)*sqrt(3*sqrt(2) "
"+ sqrt(-17*sqrt(2) + 26) - 3)*sqrt(sqrt(sqrt(2) + 2) - 1) - 10*(89*sqrt(2) "
"- 126)*sqrt(sqrt(2) + 2) + 2*((61*sqrt(2) - 85)*sqrt(sqrt(2) + 2) - 3*((2*s"
"qrt(2) - 3)*sqrt(sqrt(2) + 2) - 2*sqrt(2) + 3)*sqrt(-17*sqrt(2) + 26) - 61*"
"sqrt(2) + 85)*sqrt(-sqrt(2) + 2) + 22*((5*sqrt(2) - 7)*sqrt(sqrt(2) + 2) - "
"5*sqrt(2) + 7)*sqrt(-17*sqrt(2) + 26) + 890*sqrt(2) - 1260)*sqrt(-12*sqrt(2"
") - 2*sqrt(-sqrt(2) + 2) - 2*sqrt(-17*sqrt(2) + 26) + 24) + 4*((319*sqrt(2)"
" - 452)*sqrt(sqrt(2) + 2) - 561*sqrt(2) + 794)*sqrt(-17*sqrt(2) + 26) + 170"
"64*sqrt(2) - 24132) + ((5*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2) - ((61*sqrt("
"2) - 85)*sqrt(sqrt(2) + 2) - 3*((2*sqrt(2) - 3)*sqrt(sqrt(2) + 2) - 4*sqrt("
"2) + 6)*sqrt(-17*sqrt(2) + 26) - 122*sqrt(2) + 170)*sqrt(-sqrt(2) + 2) - 11"
"*((5*sqrt(2) - 7)*sqrt(sqrt(2) + 2) - 10*sqrt(2) + 14)*sqrt(-17*sqrt(2) + 2"
"6) - 890*sqrt(2) + 1260)*sqrt(3*sqrt(2) + sqrt(-sqrt(2) + 2) - 5)*sqrt(-12*"
"sqrt(2) - 2*sqrt(-sqrt(2) + 2) - 2*sqrt(-17*sqrt(2) + 26) + 24) + 2*(10*(63"
"*sqrt(2) - 89)*sqrt(sqrt(2) + 2) - ((85*sqrt(2) - 122)*sqrt(sqrt(2) + 2) - "
"3*((3*sqrt(2) - 4)*sqrt(sqrt(2) + 2) - 6*sqrt(2) + 8)*sqrt(-17*sqrt(2) + 26"
") - 170*sqrt(2) + 244)*sqrt(-sqrt(2) + 2) - 11*((7*sqrt(2) - 10)*sqrt(sqrt("
"2) + 2) - 14*sqrt(2) + 20)*sqrt(-17*sqrt(2) + 26) - 1260*sqrt(2) + 1780)*sq"
"rt(3*sqrt(2) + sqrt(-sqrt(2) + 2) - 5))*((sqrt(2) + sqrt(sqrt(2) + 2))*sqrt"
"(sqrt(sqrt(2) + 2) - 1) + 3*sqrt(2) - 5*sqrt(sqrt(2) + 2) + 8)*sqrt(sqrt(sq"
"rt(2) + 2) + 2)/((2*((3*(3*sqrt(2) - 4)*sqrt(-17*sqrt(2) + 26) - 85*sqrt(2)"
" + 122)*sqrt(-sqrt(2) + 2) - 11*(7*sqrt(2) - 10)*sqrt(-17*sqrt(2) + 26) + 6"
"30*sqrt(2) - 890)*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3)*sqrt(sqrt(sq"
"rt(2) + 2) - 1) - 2*(4896*sqrt(2) - 6923)*sqrt(sqrt(2) + 2) + (20*(79*sqrt("
"2) - 112)*sqrt(sqrt(2) + 2) - (7*(27*sqrt(2) - 38)*sqrt(sqrt(2) + 2) - 342*"
"sqrt(2) + 484)*sqrt(-17*sqrt(2) + 26) - 2820*sqrt(2) + 3992)*sqrt(-sqrt(2) "
"+ 2) + (((3*(2*sqrt(2) - 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt("
"-sqrt(2) + 2) - 11*(5*sqrt(2) - 7)*sqrt(-17*sqrt(2) + 26) + 445*sqrt(2) - 6"
"30)*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3)*sqrt(sqrt(sqrt(2) + 2) - 1"
") - 10*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2) + 2*((61*sqrt(2) - 85)*sqrt(sqr"
"t(2) + 2) - 3*((2*sqrt(2) - 3)*sqrt(sqrt(2) + 2) - 2*sqrt(2) + 3)*sqrt(-17*"
"sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt(2) + 2) + 22*((5*sqrt(2) - 7)*s"
"qrt(sqrt(2) + 2) - 5*sqrt(2) + 7)*sqrt(-17*sqrt(2) + 26) + 890*sqrt(2) - 12"
"60)*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt(2) + 2) - 2*sqrt(-17*sqrt(2) + 26) + 24"
") + 4*((319*sqrt(2) - 452)*sqrt(sqrt(2) + 2) - 561*sqrt(2) + 794)*sqrt(-17*"
"sqrt(2) + 26) + 17064*sqrt(2) - 24132)*((sqrt(sqrt(2) + 2) + 1)*sqrt(sqrt(s"
"qrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1)*((sqrt(sqrt(2) + 2) + 1)*sqrt(sqr"
"t(sqrt(2) + 2) - 1) - sqrt(sqrt(2) + 2) - ((sqrt(2) + sqrt(sqrt(2) + 2))*sq"
"rt(sqrt(sqrt(2) + 2) - 1) + 3*sqrt(2) - 5*sqrt(sqrt(2) + 2) + 8)*((sqrt(2) "
"+ sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2) - 1) - 3*sqrt(2) + 5*sqrt(sqrt("
"2) + 2) - 8)*(sqrt(sqrt(2) + 2) + 2)/(((sqrt(sqrt(2) + 2) + 1)*sqrt(sqrt(sq"
"rt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1)*(sqrt(sqrt(2) + 2) - 2)^3) - 1)*(-"
"1/4*sqrt(sqrt(2) + 2) + 1/2)^(3/2)))*sqrt(sqrt(sqrt(2) + 2) + 2)/sqrt(-1/4*"
"sqrt(sqrt(2) + 2) + 1/2) - 2*(sqrt(sqrt(2) + 2)*(sqrt(2) - 1)*sqrt(sqrt(sqr"
"t(2) + 2) - 1) + sqrt(2) - sqrt(sqrt(2) + 2))*(8*((5*(89*sqrt(2) - 126)*sqr"
"t(sqrt(2) + 2) - ((61*sqrt(2) - 85)*sqrt(sqrt(2) + 2) - 3*((2*sqrt(2) - 3)*"
"sqrt(sqrt(2) + 2) - 4*sqrt(2) + 6)*sqrt(-17*sqrt(2) + 26) - 122*sqrt(2) + 1"
"70)*sqrt(-sqrt(2) + 2) - 11*((5*sqrt(2) - 7)*sqrt(sqrt(2) + 2) - 10*sqrt(2)"
" + 14)*sqrt(-17*sqrt(2) + 26) - 890*sqrt(2) + 1260)*sqrt(3*sqrt(2) + sqrt(-"
"sqrt(2) + 2) - 5)*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt(2) + 2) - 2*sqrt(-17*sqrt"
"(2) + 26) + 24) + 2*(10*(63*sqrt(2) - 89)*sqrt(sqrt(2) + 2) - ((85*sqrt(2) "
"- 122)*sqrt(sqrt(2) + 2) - 3*((3*sqrt(2) - 4)*sqrt(sqrt(2) + 2) - 6*sqrt(2)"
" + 8)*sqrt(-17*sqrt(2) + 26) - 170*sqrt(2) + 244)*sqrt(-sqrt(2) + 2) - 11*("
"(7*sqrt(2) - 10)*sqrt(sqrt(2) + 2) - 14*sqrt(2) + 20)*sqrt(-17*sqrt(2) + 26"
") - 1260*sqrt(2) + 1780)*sqrt(3*sqrt(2) + sqrt(-sqrt(2) + 2) - 5))*(1/((sqr"
"t(sqrt(2) + 2) + 1)*sqrt(sqrt(sqrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1) + "
"((sqrt(2) + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2) - 1) + 3*sqrt(2) - 5*"
"sqrt(sqrt(2) + 2) + 8)*((sqrt(2) + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2"
") - 1) - 3*sqrt(2) + 5*sqrt(sqrt(2) + 2) - 8)*(sqrt(sqrt(2) + 2) + 2)/(((sq"
"rt(sqrt(2) + 2) + 1)*sqrt(sqrt(sqrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1)^2"
"*((sqrt(sqrt(2) + 2) + 1)*sqrt(sqrt(sqrt(2) + 2) - 1) - sqrt(sqrt(2) + 2) -"
" ((sqrt(2) + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2) - 1) + 3*sqrt(2) - 5"
"*sqrt(sqrt(2) + 2) + 8)*((sqrt(2) + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + "
"2) - 1) - 3*sqrt(2) + 5*sqrt(sqrt(2) + 2) - 8)*(sqrt(sqrt(2) + 2) + 2)/(((s"
"qrt(sqrt(2) + 2) + 1)*sqrt(sqrt(sqrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1)*"
"(sqrt(sqrt(2) + 2) - 2)^3) - 1)*(sqrt(sqrt(2) + 2) - 2)^3))/(2*((3*(3*sqrt("
"2) - 4)*sqrt(-17*sqrt(2) + 26) - 85*sqrt(2) + 122)*sqrt(-sqrt(2) + 2) - 11*"
"(7*sqrt(2) - 10)*sqrt(-17*sqrt(2) + 26) + 630*sqrt(2) - 890)*sqrt(3*sqrt(2)"
" + sqrt(-17*sqrt(2) + 26) - 3)*sqrt(sqrt(sqrt(2) + 2) - 1) - 2*(4896*sqrt(2"
") - 6923)*sqrt(sqrt(2) + 2) + (20*(79*sqrt(2) - 112)*sqrt(sqrt(2) + 2) - (7"
"*(27*sqrt(2) - 38)*sqrt(sqrt(2) + 2) - 342*sqrt(2) + 484)*sqrt(-17*sqrt(2) "
"+ 26) - 2820*sqrt(2) + 3992)*sqrt(-sqrt(2) + 2) + (((3*(2*sqrt(2) - 3)*sqrt"
"(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt(2) + 2) - 11*(5*sqrt(2) - "
"7)*sqrt(-17*sqrt(2) + 26) + 445*sqrt(2) - 630)*sqrt(3*sqrt(2) + sqrt(-17*sq"
"rt(2) + 26) - 3)*sqrt(sqrt(sqrt(2) + 2) - 1) - 10*(89*sqrt(2) - 126)*sqrt(s"
"qrt(2) + 2) + 2*((61*sqrt(2) - 85)*sqrt(sqrt(2) + 2) - 3*((2*sqrt(2) - 3)*s"
"qrt(sqrt(2) + 2) - 2*sqrt(2) + 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)"
"*sqrt(-sqrt(2) + 2) + 22*((5*sqrt(2) - 7)*sqrt(sqrt(2) + 2) - 5*sqrt(2) + 7"
")*sqrt(-17*sqrt(2) + 26) + 890*sqrt(2) - 1260)*sqrt(-12*sqrt(2) - 2*sqrt(-s"
"qrt(2) + 2) - 2*sqrt(-17*sqrt(2) + 26) + 24) + 4*((319*sqrt(2) - 452)*sqrt("
"sqrt(2) + 2) - 561*sqrt(2) + 794)*sqrt(-17*sqrt(2) + 26) + 17064*sqrt(2) - "
"24132) + (44*(7*sqrt(2) - 10)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) + 2*"
"(11*(7*sqrt(2) - 10)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - 10*(63*sqrt"
"(2) - 89)*sqrt(sqrt(2) + 2) - (3*(3*sqrt(2) - 4)*sqrt(sqrt(2) + 2)*sqrt(-17"
"*sqrt(2) + 26) - (85*sqrt(2) - 122)*sqrt(sqrt(2) + 2))*sqrt(-sqrt(2) + 2) +"
" 2*((3*(3*sqrt(2) - 4)*sqrt(-17*sqrt(2) + 26) - 85*sqrt(2) + 122)*sqrt(-sqr"
"t(2) + 2) - 11*(7*sqrt(2) - 10)*sqrt(-17*sqrt(2) + 26) + 630*sqrt(2) - 890)"
"*sqrt(sqrt(sqrt(2) + 2) - 1))*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3) "
"- 40*(63*sqrt(2) - 89)*sqrt(sqrt(2) + 2) - 4*(3*(3*sqrt(2) - 4)*sqrt(sqrt(2"
") + 2)*sqrt(-17*sqrt(2) + 26) - (85*sqrt(2) - 122)*sqrt(sqrt(2) + 2))*sqrt("
"-sqrt(2) + 2) + (22*(5*sqrt(2) - 7)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26"
") + (11*(5*sqrt(2) - 7)*sqrt(sqrt(2) + 2)*sqrt(-17*sqrt(2) + 26) - 5*(89*sq"
"rt(2) - 126)*sqrt(sqrt(2) + 2) - (3*(2*sqrt(2) - 3)*sqrt(sqrt(2) + 2)*sqrt("
"-17*sqrt(2) + 26) - (61*sqrt(2) - 85)*sqrt(sqrt(2) + 2))*sqrt(-sqrt(2) + 2)"
" + 2*((3*(2*sqrt(2) - 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sq"
"rt(2) + 2) - 11*(5*sqrt(2) - 7)*sqrt(-17*sqrt(2) + 26) + 445*sqrt(2) - 630)"
"*sqrt(sqrt(sqrt(2) + 2) - 1))*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3) "
"- 10*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2) - 2*(3*(2*sqrt(2) - 3)*sqrt(sqrt("
"2) + 2)*sqrt(-17*sqrt(2) + 26) - (61*sqrt(2) - 85)*sqrt(sqrt(2) + 2))*sqrt("
"-sqrt(2) + 2) + 4*((3*(2*sqrt(2) - 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) +"
" 85)*sqrt(-sqrt(2) + 2) - 11*(5*sqrt(2) - 7)*sqrt(-17*sqrt(2) + 26) + 445*s"
"qrt(2) - 630)*sqrt(sqrt(sqrt(2) + 2) - 1))*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt("
"2) + 2) - 2*sqrt(-17*sqrt(2) + 26) + 24) + 8*((3*(3*sqrt(2) - 4)*sqrt(-17*s"
"qrt(2) + 26) - 85*sqrt(2) + 122)*sqrt(-sqrt(2) + 2) - 11*(7*sqrt(2) - 10)*s"
"qrt(-17*sqrt(2) + 26) + 630*sqrt(2) - 890)*sqrt(sqrt(sqrt(2) + 2) - 1))*((s"
"qrt(2) + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2) - 1) + 3*sqrt(2) - 5*sqr"
"t(sqrt(2) + 2) + 8)*sqrt(sqrt(sqrt(2) + 2) + 2)/((2*((3*(3*sqrt(2) - 4)*sqr"
"t(-17*sqrt(2) + 26) - 85*sqrt(2) + 122)*sqrt(-sqrt(2) + 2) - 11*(7*sqrt(2) "
"- 10)*sqrt(-17*sqrt(2) + 26) + 630*sqrt(2) - 890)*sqrt(3*sqrt(2) + sqrt(-17"
"*sqrt(2) + 26) - 3)*sqrt(sqrt(sqrt(2) + 2) - 1) - 2*(4896*sqrt(2) - 6923)*s"
"qrt(sqrt(2) + 2) + (20*(79*sqrt(2) - 112)*sqrt(sqrt(2) + 2) - (7*(27*sqrt(2"
") - 38)*sqrt(sqrt(2) + 2) - 342*sqrt(2) + 484)*sqrt(-17*sqrt(2) + 26) - 282"
"0*sqrt(2) + 3992)*sqrt(-sqrt(2) + 2) + (((3*(2*sqrt(2) - 3)*sqrt(-17*sqrt(2"
") + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt(2) + 2) - 11*(5*sqrt(2) - 7)*sqrt(-17"
"*sqrt(2) + 26) + 445*sqrt(2) - 630)*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26)"
" - 3)*sqrt(sqrt(sqrt(2) + 2) - 1) - 10*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2)"
" + 2*((61*sqrt(2) - 85)*sqrt(sqrt(2) + 2) - 3*((2*sqrt(2) - 3)*sqrt(sqrt(2)"
" + 2) - 2*sqrt(2) + 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt"
"(2) + 2) + 22*((5*sqrt(2) - 7)*sqrt(sqrt(2) + 2) - 5*sqrt(2) + 7)*sqrt(-17*"
"sqrt(2) + 26) + 890*sqrt(2) - 1260)*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt(2) + 2)"
" - 2*sqrt(-17*sqrt(2) + 26) + 24) + 4*((319*sqrt(2) - 452)*sqrt(sqrt(2) + 2"
") - 561*sqrt(2) + 794)*sqrt(-17*sqrt(2) + 26) + 17064*sqrt(2) - 24132)*((sq"
"rt(sqrt(2) + 2) + 1)*sqrt(sqrt(sqrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1)*("
"(sqrt(sqrt(2) + 2) + 1)*sqrt(sqrt(sqrt(2) + 2) - 1) - sqrt(sqrt(2) + 2) - ("
"(sqrt(2) + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2) - 1) + 3*sqrt(2) - 5*s"
"qrt(sqrt(2) + 2) + 8)*((sqrt(2) + sqrt(sqrt(2) + 2))*sqrt(sqrt(sqrt(2) + 2)"
" - 1) - 3*sqrt(2) + 5*sqrt(sqrt(2) + 2) - 8)*(sqrt(sqrt(2) + 2) + 2)/(((sqr"
"t(sqrt(2) + 2) + 1)*sqrt(sqrt(sqrt(2) + 2) - 1) + sqrt(sqrt(2) + 2) + 1)*(s"
"qrt(sqrt(2) + 2) - 2)^3) - 1)*(-1/4*sqrt(sqrt(2) + 2) + 1/2)^(3/2)))/(sqrt("
"sqrt(2) + 2) - 2))/(2*((3*(3*sqrt(2) - 4)*sqrt(-17*sqrt(2) + 26) - 85*sqrt("
"2) + 122)*sqrt(-sqrt(2) + 2) - 11*(7*sqrt(2) - 10)*sqrt(-17*sqrt(2) + 26) +"
" 630*sqrt(2) - 890)*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3)*sqrt(sqrt("
"sqrt(2) + 2) - 1) - 2*(4896*sqrt(2) - 6923)*sqrt(sqrt(2) + 2) + (20*(79*sqr"
"t(2) - 112)*sqrt(sqrt(2) + 2) - (7*(27*sqrt(2) - 38)*sqrt(sqrt(2) + 2) - 34"
"2*sqrt(2) + 484)*sqrt(-17*sqrt(2) + 26) - 2820*sqrt(2) + 3992)*sqrt(-sqrt(2"
") + 2) + (((3*(2*sqrt(2) - 3)*sqrt(-17*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqr"
"t(-sqrt(2) + 2) - 11*(5*sqrt(2) - 7)*sqrt(-17*sqrt(2) + 26) + 445*sqrt(2) -"
" 630)*sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3)*sqrt(sqrt(sqrt(2) + 2) -"
" 1) - 10*(89*sqrt(2) - 126)*sqrt(sqrt(2) + 2) + 2*((61*sqrt(2) - 85)*sqrt(s"
"qrt(2) + 2) - 3*((2*sqrt(2) - 3)*sqrt(sqrt(2) + 2) - 2*sqrt(2) + 3)*sqrt(-1"
"7*sqrt(2) + 26) - 61*sqrt(2) + 85)*sqrt(-sqrt(2) + 2) + 22*((5*sqrt(2) - 7)"
"*sqrt(sqrt(2) + 2) - 5*sqrt(2) + 7)*sqrt(-17*sqrt(2) + 26) + 890*sqrt(2) - "
"1260)*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt(2) + 2) - 2*sqrt(-17*sqrt(2) + 26) + "
"24) + 4*((319*sqrt(2) - 452)*sqrt(sqrt(2) + 2) - 561*sqrt(2) + 794)*sqrt(-1"
"7*sqrt(2) + 26) + 17064*sqrt(2) - 24132)";

#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

const char * EXPR_M =
"-(4*(6*sqrt(2) + sqrt(-sqrt(2) + 2) + sqrt(-17*sqrt(2) + 26) - 8)*sqrt(3*sq"
"rt(2) + sqrt(-sqrt(2) + 2) - 5) - sqrt(3*sqrt(2) + sqrt(-17*sqrt(2) + 26) -"
" 3)*(-24*I*sqrt(2) - 4*I*sqrt(-sqrt(2) + 2) - 4*I*sqrt(-17*sqrt(2) + 26) + "
"32*I) - ((sqrt(2)*sqrt(-sqrt(2) + 2) + sqrt(2)*sqrt(-17*sqrt(2) + 26) - 8*s"
"qrt(2) + 12)*sqrt(3*sqrt(2) + sqrt(-sqrt(2) + 2) - 5) + (I*sqrt(2)*sqrt(-sq"
"rt(2) + 2) + I*sqrt(2)*sqrt(-17*sqrt(2) + 26) - 8*I*sqrt(2) + 12*I)*sqrt(3*"
"sqrt(2) + sqrt(-17*sqrt(2) + 26) - 3))*sqrt(-12*sqrt(2) - 2*sqrt(-sqrt(2) +"
" 2) - 2*sqrt(-17*sqrt(2) + 26) + 24) - ((24*I*sqrt(2) + 4*I*sqrt(-17*sqrt(2"
") + 26) - 32*I)*sqrt(-sqrt(2) + 2) + 8*I*(3*sqrt(2) - 4)*sqrt(-17*sqrt(2) +"
" 26) - 228*I*sqrt(2) + 328*I)*sqrt(sqrt(sqrt(2) + 2) - 1))/(4*(6*sqrt(2) + "
"sqrt(-sqrt(2) + 2) + sqrt(-17*sqrt(2) + 26) - 8)*sqrt(3*sqrt(2) + sqrt(-17*"
"sqrt(2) + 26) - 3)*sqrt(sqrt(sqrt(2) + 2) - 1) + sqrt(3*sqrt(2) + sqrt(-sqr"
"t(2) + 2) - 5)*(-24*I*sqrt(2) - 4*I*sqrt(-sqrt(2) + 2) - 4*I*sqrt(-17*sqrt("
"2) + 26) + 32*I)*sqrt(sqrt(sqrt(2) + 2) - 1) - 4*(6*sqrt(2) + sqrt(-17*sqrt"
"(2) + 26) - 8)*sqrt(-sqrt(2) + 2) + ((I*sqrt(2)*sqrt(-sqrt(2) + 2) + I*sqrt"
"(2)*sqrt(-17*sqrt(2) + 26) - 8*I*sqrt(2) + 12*I)*sqrt(3*sqrt(2) + sqrt(-sqr"
"t(2) + 2) - 5)*sqrt(sqrt(sqrt(2) + 2) - 1) - (sqrt(2)*sqrt(-sqrt(2) + 2) + "
"sqrt(2)*sqrt(-17*sqrt(2) + 26) - 8*sqrt(2) + 12)*sqrt(3*sqrt(2) + sqrt(-17*"
"sqrt(2) + 26) - 3)*sqrt(sqrt(sqrt(2) + 2) - 1))*sqrt(-12*sqrt(2) - 2*sqrt(-"
"sqrt(2) + 2) - 2*sqrt(-17*sqrt(2) + 26) + 24) - 8*(3*sqrt(2) - 4)*sqrt(-17*"
"sqrt(2) + 26) + 228*sqrt(2) - 328)";


void doit(gr_ctx_t ctx)
{
    gr_ptr N, M, E;
    truth_t equal;
    int status = GR_SUCCESS;

    gr_ctx_println(ctx);

    GR_TMP_INIT3(N, M, E, ctx);

    flint_printf("Evaluating N...\n");
    TIMEIT_ONCE_START
    GR_MUST_SUCCEED(gr_set_str(N, EXPR_N, ctx));
    TIMEIT_ONCE_STOP

    flint_printf("Evaluating M...\n");
    TIMEIT_ONCE_START
    GR_MUST_SUCCEED(gr_set_str(M, EXPR_M, ctx));
    TIMEIT_ONCE_STOP

    flint_printf("Evaluating E = -(1-|M|^2)^2...\n");
    TIMEIT_ONCE_START
    GR_MUST_SUCCEED(gr_abs(E, M, ctx));
    GR_MUST_SUCCEED(gr_pow_ui(E, E, 2, ctx));
    GR_MUST_SUCCEED(gr_sub_si(E, E, 1, ctx));
    GR_MUST_SUCCEED(gr_neg(E, E, ctx));
    GR_MUST_SUCCEED(gr_pow_ui(E, E, 2, ctx));
    GR_MUST_SUCCEED(gr_neg(E, E, ctx));
    TIMEIT_ONCE_STOP

    if (status != GR_SUCCESS)
    {
        flint_printf("evaluation failed\n");
        flint_abort();
    }

    flint_printf("N ~ "); gr_println(N, ctx);
    flint_printf("E ~ "); gr_println(E, ctx);

    flint_printf("Testing E = N...\n");
    TIMEIT_ONCE_START
    equal = gr_equal(E, N, ctx);
    TIMEIT_ONCE_STOP

    flint_printf("\nEqual = ");
    truth_print(equal);
    flint_printf("\n");

    GR_TMP_CLEAR3(N, M, E, ctx);
}

int main(int argc, char *argv[])
{
    TIMEIT_ONCE_START

    if (argc >= 2 && strcmp(argv[1], "-ca") == 0)
    {
        gr_ctx_t ctx;
        gr_ctx_init_complex_algebraic_ca(ctx);
        gr_ctx_ca_set_option(ctx, CA_OPT_QQBAR_DEG_LIMIT, 10000);
        doit(ctx);
        gr_ctx_clear(ctx);
    }
    else
    {
        gr_ctx_t ctx;
        gr_ctx_init_complex_qqbar(ctx);
        doit(ctx);
        gr_ctx_clear(ctx);
    }

    flint_printf("\n");
    flint_printf("Total: ");
    TIMEIT_ONCE_STOP
    SHOW_MEMORY_USAGE

    flint_cleanup_master();
    return 0;
}
