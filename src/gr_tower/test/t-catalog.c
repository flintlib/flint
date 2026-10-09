/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdio.h>
#include "test_helpers.h"
#include "fmpz_poly.h"
#include "acb.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

/*
    A catalog of zero-testing problems collected from the literature and
    from the test suites of other systems (Calcium, SymPy, mpmath, Sage,
    Wester's problem set, denesting papers, near-miss collections), parsed
    with gr_set_str in the lazy tower field (src/python/gr_tower_profile.py
    has the sources and the cases not included here: RootOf, relations,
    most special functions (those with several arguments, which the parser
    does not read: hypergeometric and theta functions, Eisenstein series),
    programs, the big expressions). The prefix of each name says where the case comes from
    (ca: Calcium/FLINT examples, sp: SymPy, mp: SymPy minpoly, dn:
    denesting, w: Wester, tr: trigonometric values, el: elementary
    functions, nm: near misses, sage: Sage QQbar).

        'Z'  the expression is zero
        'N'  the expression is not zero (the zero test must say so, never
             that it is zero), and its value is within 5% of re + im i
        'P'  P(x) = 0 for the expression x and the polynomial P with the
             given coefficients (in ascending order)

    The cases run once in a fresh field each, then in sequences which
    went wrong, then together in a shared field in a random order (so
    that each one meets generators left over by the others), followed by
    scalable families: Gauss sums, products of sines, Swinnerton-Dyer
    polynomials and the exact DFT benchmark.
*/

typedef struct
{
    const char * name;
    char kind;
    const char * re;
    const char * im;
    const char * poly;
    const char * expr;
}
catalog_case_struct;

static const catalog_case_struct catalog_cases[] =
{
    {"ca.euler", 'Z', NULL, NULL, NULL,
     "(exp((pi*i))+1)"},
    {"ca.log_m1", 'Z', NULL, NULL, NULL,
     "((log((-1))/(pi*i))-1)"},
    {"ca.log_mi", 'Z', NULL, NULL, NULL,
     "((log((-i))/(pi*i))+(1/2))"},
    {"ca.log_pow10", 'Z', NULL, NULL, NULL,
     "((log((1/(10^123)))/log(100))+(123/2))"},
    {"ca.log_unit", 'Z', NULL, NULL, NULL,
     "((log((1+sqrt(2)))/log((3+(2*sqrt(2)))))-(1/2))"},
    {"ca.sqrt6", 'Z', NULL, NULL, NULL,
     "((sqrt(2)*sqrt(3))-sqrt(6))"},
    {"ca.exp_sum", 'Z', NULL, NULL, NULL,
     "(((exp((1+sqrt(2)))*exp((1-sqrt(2))))/(exp(1)^2))-1)"},
    {"ca.i_to_i", 'Z', NULL, NULL, NULL,
     "((i^i)-exp(((-pi)/2)))"},
    {"ca.exp_sqrt12", 'Z', NULL, NULL, NULL,
     "((exp(sqrt(3))^2)-exp(sqrt(12)))"},
    {"ca.log_pi_i", 'Z', NULL, NULL, NULL,
     "(((2*log((pi*i)))-(4*log(sqrt(pi))))-(pi*i))"},
    {"ca.bbk_ex1", 'Z', NULL, NULL, NULL,
     "((((((((((-i)*pi)/8)*(log(((2/3)-((2*i)/3)))^2))+(((i*pi)/8)*(log(((2/3)+((2*i)/3)))^2)))+(((pi^2)/1"
     "2)*log(((-1)-i))))+(((pi^2)/12)*log(((-1)+i))))+(((pi^2)/12)*log(((1/3)-(i/3)))))+(((pi^2)/12)*log(("
     "(1/3)+(i/3)))))+(((pi^2)/48)*log(18)))"},
    {"ca.denest5_2_6", 'Z', NULL, NULL, NULL,
     "((sqrt((5+(2*sqrt(6))))-sqrt(2))-sqrt(3))"},
    {"ca.sqrt_i", 'Z', NULL, NULL, NULL,
     "(sqrt(i)-((1+i)/sqrt(2)))"},
    {"ca.erf_logs", 'Z', NULL, NULL, NULL,
     "(erf(((2*log(sqrt(((1/2)-(sqrt(2)/4)))))+log(4)))-erf(log((2-sqrt(2)))))"},
    {"ca.iter_pow", 'Z', NULL, NULL, NULL,
     "((i^i)-exp((pi/((sqrt((-2))^sqrt(2))^sqrt(2)))))"},
    {"ca.trig_tr1", 'Z', NULL, NULL, NULL,
     "(((sin((sqrt(2)/2))^2)+(cos((1/sqrt(2)))^2))-1)"},
    {"ca.trig_tr2", 'Z', NULL, NULL, NULL,
     "(sin((3+pi))+sin(3))"},
    {"ca.trig_tr3", 'Z', NULL, NULL, NULL,
     "(tan((1+pi))-tan(1))"},
    {"ca.gd", 'Z', NULL, NULL, NULL,
     "(sin((2*atan(tanh((1/2)))))-tanh(1))"},
    {"ca.atan_alg", 'Z', NULL, NULL, NULL,
     "(atan((1-sqrt(2)))+(pi/8))"},
    {"ca.atan_tan", 'Z', NULL, NULL, NULL,
     "(atan(tan(((23*pi)/27)))+((4*pi)/27))"},
    {"ca.asin_sin", 'Z', NULL, NULL, NULL,
     "(asin(sin((sqrt(2)-1)))-(sqrt(2)-1))"},
    {"ca.gamma_fe", 'Z', NULL, NULL, NULL,
     "((gamma((pi+1))/gamma(pi))-pi)"},
    {"ca.erf_erfc", 'Z', NULL, NULL, NULL,
     "((erf(exp(((pi*i)/3)))-erfc(exp(((((-2)*pi)*i)/3))))+1)"},
    {"ca.mixed1", 'Z', NULL, NULL, NULL,
     "((((pi+sqrt(2))+sqrt(3))/(pi+sqrt((5+(2*sqrt(6))))))-1)"},
    {"ca.mixed2", 'Z', NULL, NULL, NULL,
     "((log((1/exp((sqrt(2)+1))))+sqrt(2))+1)"},
    {"ca.mixed3", 'Z', NULL, NULL, NULL,
     "(arg(sqrt(((-pi)*i)))+(pi/4))"},
    {"ca.mixed4", 'Z', NULL, NULL, NULL,
     "(sin(((1+sqrt(2))/2))-sqrt(((1-cos((1+sqrt(2))))/2)))"},
    {"ca.gosper", 'Z', NULL, NULL, NULL,
     "((((sqrt(((36+((3*(((-54)+((35*i)*sqrt(3)))^(1/3)))*(3^(1/3))))+(117/(((-162)+((105*i)*sqrt(3)))^(1/"
     "3)))))/3)+((sqrt(5)*(((((1296*i)+(840*sqrt(3)))-((35*(3^(5/6)))*(((-54)+((35*i)*sqrt(3)))^(1/3))))-("
     "(54*i)*(((-162)+((105*i)*sqrt(3)))^(1/3))))+((13*i)*(((-162)+((105*i)*sqrt(3)))^(2/3)))))/(5*((162*i"
     ")+(105*sqrt(3))))))-sqrt(5))-sqrt(7))"},
    {"ca.ramanujan_163", 'N', "-7.4993e-13", "0", NULL,
     "(exp((pi*sqrt(163)))-((640320^3)+744))"},
    {"ca.exp_tiny", 'N', "1.0e-10000", "0", NULL,
     "(exp((1/(10^10000)))-1)"},
    {"ca.machin", 'Z', NULL, NULL, NULL,
     "(((4*atan((1/5)))-atan((1/239)))-(pi/4))"},
    {"ca.machin2", 'Z', NULL, NULL, NULL,
     "((atan((1/2))+atan((1/3)))-(pi/4))"},
    {"ca.machin3", 'Z', NULL, NULL, NULL,
     "(((2*atan((1/2)))-atan((1/7)))-(pi/4))"},
    {"ca.machin4", 'Z', NULL, NULL, NULL,
     "(((2*atan((1/3)))+atan((1/7)))-(pi/4))"},
    {"ca.machin5", 'Z', NULL, NULL, NULL,
     "(((atan((1/2))+atan((1/5)))+atan((1/8)))-(pi/4))"},
    {"ca.machin6", 'Z', NULL, NULL, NULL,
     "((((atan((1/3))+atan((1/4)))+atan((1/7)))+atan((1/13)))-(pi/4))"},
    {"ca.machin7", 'Z', NULL, NULL, NULL,
     "(((((12*atan((1/49)))+(32*atan((1/57))))-(5*atan((1/239))))+(12*atan((1/110443))))-(pi/4))"},
    {"ca.hmachin2", 'Z', NULL, NULL, NULL,
     "((((14*atanh((1/31)))+(10*atanh((1/49))))+(6*atanh((1/161))))-log(2))"},
    {"ca.hmachin3", 'Z', NULL, NULL, NULL,
     "((((22*atanh((1/31)))+(16*atanh((1/49))))+(10*atanh((1/161))))-log(3))"},
    {"ca.hmachin5", 'Z', NULL, NULL, NULL,
     "((((32*atanh((1/31)))+(24*atanh((1/49))))+(14*atanh((1/161))))-log(5))"},
    {"ca.hmachin7", 'Z', NULL, NULL, NULL,
     "(((((404*atanh((1/251)))+(152*atanh((1/449))))-(106*atanh((1/4801))))+(174*atanh((1/8749))))-log(7))"},
    {"ca.hmachin2b", 'Z', NULL, NULL, NULL,
     "(((((144*atanh((1/251)))+(54*atanh((1/449))))-(38*atanh((1/4801))))+(62*atanh((1/8749))))-log(2))"},
    {"sp.equals1", 'Z', NULL, NULL, NULL,
     "(((-3)-sqrt(5))+((((-sqrt(10))/2)-(sqrt(2)/2))^2))"},
    {"sp.equals2", 'Z', NULL, NULL, NULL,
     "(((-((-1)^(3/4)))*(6^(1/4)))+(((-6)^(1/4))*i))"},
    {"sp.equals3", 'Z', NULL, NULL, NULL,
     "((sqrt((1+sqrt(3)))+sqrt((3+(3*sqrt(3)))))-sqrt((10+(6*sqrt(3)))))"},
    {"sp.equals4", 'Z', NULL, NULL, NULL,
     "(((((3^(1/3))+3)^3)^(1/3))-((3^(1/3))+3))"},
    {"sp.equals_branch_zero", 'Z', NULL, NULL, NULL,
     "((((((2*sqrt(2))*((-1)^(5/2)))*((1+(1/(2*(-1))))^(5/2)))/5)+((((2*sqrt(2))*((-1)^(3/2)))*((1+(1/(2*("
     "-1))))^(5/2)))/((-6)-(3/(-1)))))-((sqrt(((2*(-1))+1))*(((6*((-1)^2))+(-1))-1))/15))"},
    {"sp.equals_branch_nonzero", 'Z', NULL, NULL, NULL,
     "(((((((2*sqrt(2))*(((-1)/4)^(5/2)))*((1+(1/(2*((-1)/4))))^(5/2)))/5)+((((2*sqrt(2))*(((-1)/4)^(3/2))"
     ")*((1+(1/(2*((-1)/4))))^(5/2)))/((-6)-(3/((-1)/4)))))-((sqrt(((2*((-1)/4))+1))*(((6*(((-1)/4)^2))+(("
     "-1)/4))-1))/15))-((7*sqrt(2))/120))"},
    {"sp.cardano93a", 'Z', NULL, NULL, NULL,
     "(((((((-(2^(1/3)))*(((3*sqrt(93))+29)^2))-(4*(((3*sqrt(93))+29)^(4/3))))+((12*sqrt(93))*(((3*sqrt(93"
     "))+29)^(1/3))))+(116*(((3*sqrt(93))+29)^(1/3))))+((174*(2^(1/3)))*sqrt(93)))+(1678*(2^(1/3))))"},
    {"sp.cardano93b", 'Z', NULL, NULL, NULL,
     "((((9*(((3*sqrt(93))+29)^(2/3)))*((((((3*sqrt(93))+29)^(1/3))*(((-(2^(2/3)))*(((3*sqrt(93))+29)^(1/3"
     ")))-2))-(2*(2^(1/3))))^3))+((72*(((3*sqrt(93))+29)^(2/3)))*((81*sqrt(93))+783)))+(((162*sqrt(93))+15"
     "66)*((((((3*sqrt(93))+29)^(1/3))*(((-(2^(2/3)))*(((3*sqrt(93))+29)^(1/3)))-2))-(2*(2^(1/3))))^2)))"},
    {"sp.trig90", 'Z', NULL, NULL, NULL,
     "((((2*(((-3)*tan(((19*pi)/90)))+sqrt(3)))*cos(((11*pi)/90)))*cos(((19*pi)/90)))-(sqrt(3)*((-3)+(4*(c"
     "os(((19*pi)/90))^2)))))"},
    {"sp.issue4956_num", 'Z', NULL, NULL, NULL,
     "(((((-27)*(12^(1/3)))*sqrt(31))*i)+((((27*(2^(2/3)))*(3^(1/3)))*sqrt(31))*i))"},
    {"sp.issue4956_den", 'N', "-6.9122e4", "3.9971e4", NULL,
     "((((-2511)*(2^(2/3)))*(3^(1/3)))+((((29*(18^(1/3)))+((((9*(2^(1/3)))*(3^(2/3)))*sqrt(31))*i))+(((87*"
     "(2^(1/3)))*(3^(1/6)))*i))^2))"},
    {"sp.hyperbolic_nullspace", 'Z', NULL, NULL, NULL,
     "((((-exp(1))-(2*cosh((1/3))))*(((-2)*cosh((1/3)))-exp((-1))))-(((4*(cosh((1/3))^2))-1)^2))"},
    {"sp.binet_exact", 'Z', NULL, NULL, NULL,
     "(434665576869374564356885276750406258025646605173717804024817290895365554179490518904038798400792551"
     "6929592259308032263477520968962323987332247116164299644090653318793829896964992851600370447613779516"
     "6849228875-(((((1+sqrt(5))/2)^1000)-((((1+sqrt(5))/2)-1)^1000))/sqrt(5)))"},
    {"sp.binet_near", 'N', "-4.6012e-210", "0", NULL,
     "(434665576869374564356885276750406258025646605173717804024817290895365554179490518904038798400792551"
     "6929592259308032263477520968962323987332247116164299644090653318793829896964992851600370447613779516"
     "6849228875-((((1+sqrt(5))/2)^1000)/sqrt(5)))"},
    {"sp.binet5000", 'N', "5.156e-1046", "0", NULL,
     "((((1+sqrt(5))^5000)/((2^5000)*sqrt(5)))-38789684543883256337019163083259053120821277146462451061605"
     "9721489555013904403709701082291646221066947929345285888297381348310200895498294036143015691147893836"
     "4216563944106910214505634133706558656238254656700712525929903854933813928836378347518908762970712033"
     "3370529231076930085180938498018038478139967488817655546537882916442689129803846137789690215022930824"
     "7566634622492307188332480328037503913035290330450584270114763524227021093463769910400671417488329842"
     "2891491273104054328753298044273676822977244987749874555691907703880637046832794811358973739993110106"
     "2193081490185708153978543791953056175107610530756887837660336673554452588448862416192105534574936758"
     "9784902798823435102359984466393485325641195222185956306047536464547076033090242080638258492915645287"
     "6291575759142343809142302917491088984155209854432486594079793571316841692868039545309545388698114665"
     "0820668628974206393234384884652409887423958738019769938203171742089322654688793640026307977800587591"
     "29671389634214252579116872755600360311370547754724604639987588046985178408674382863125)"},
    {"sp.j_series", 'N', "1.6137e-59", "0", NULL,
     "(((1-(262537412640768744*exp(((-pi)*sqrt(163)))))-(196884*exp((((-2)*pi)*sqrt(163)))))+(103378831900"
     "730205293632*exp((((-3)*pi)*sqrt(163)))))"},
    {"sp.sin_163", 'N', "-2.356e-12", "0", NULL,
     "sin((pi*exp((pi*sqrt(163)))))"},
    {"sp.e_rational", 'N', "6.0376e-11", "0", NULL,
     "((45-((613*exp(1))/37))+(35/991))"},
    {"sp.sin2017", 'N', "2.1432e-17", "0", NULL,
     "(1+sin((2017*(2^(1/5)))))"},
    {"sp.atan_tan_pole", 'Z', NULL, NULL, NULL,
     "((atan(tan(260515))-260515)+(82924*pi))"},
    {"sp.floor_50e", 'Z', NULL, NULL, NULL,
     "(floor((30414093201713378043612608166064768844377641568960512000000000000/exp(1)))-11188719610782480"
     "504630258070757734324011354208865721592720336800)"},
    {"sp.floor_binet", 'Z', NULL, NULL, NULL,
     "(floor((((((1+sqrt(5))/2)^1000)/sqrt(5))+(1/2)))-434665576869374564356885276750406258025646605173717"
     "8040248172908953655541794905189040387984007925516929592259308032263477520968962323987332247116164299"
     "6440906533187938298969649928516003704476137795166849228875)"},
    {"sp.ceiling_pyth", 'Z', NULL, NULL, NULL,
     "(ceil((10*((sin(1)^2)+(cos(1)^2))))-10)"},
    {"sp.pi_109", 'N', "4.854e-11", "0", NULL,
     "((((((2*(2^(22/109)))*(3^(42/109)))*(5^(90/109)))*(7^(71/109)))/15)-pi)"},
    {"sp.nsimplify_cosatan", 'Z', NULL, NULL, NULL,
     "(cos(atan((1/3)))-((3*sqrt(10))/10))"},
    {"sp.nsimplify_exp_atan", 'Z', NULL, NULL, NULL,
     "(((2+exp(((2*atan((1/4)))*i)))-(49/17))-((8*i)/17))"},
    {"sp.nsimplify_root10", 'Z', NULL, NULL, NULL,
     "((1/(exp((((3*pi)*i)/5))+1))-((1/2)-(i*sqrt(((sqrt(5)/10)+(1/4))))))"},
    {"sp.gamma_reflect", 'Z', NULL, NULL, NULL,
     "((gamma((1/4))*gamma((3/4)))-(sqrt(2)*pi))"},
    {"sp.golden", 'Z', NULL, NULL, NULL,
     "((4/(1+sqrt(5)))-((-2)+(2*((1+sqrt(5))/2))))"},
    {"mp.issue26903", 'Z', NULL, NULL, NULL,
     "(sqrt(((10000000000000061^2)*10000000000000069))-(10000000000000061*sqrt(10000000000000069)))"},
    {"mp.issue14831", 'Z', NULL, NULL, NULL,
     "(((((-3)*sqrt(((12*sqrt(2))+17)))+(12*sqrt(2)))+17)-((2*sqrt(2))*sqrt(((12*sqrt(2))+17))))"},
    {"mp.issue19760", 'P', NULL, NULL, "-2 0 4 -4 1",
     "((1/(sqrt((1+sqrt(2)))-(sqrt(2)*sqrt((1+sqrt(2))))))+1)"},
    {"mp.issue6868", 'P', NULL, NULL, "71999 -48000 8000",
     "(((-1)/(800*sqrt(((((-1)/240)+(1/(18000*((((-1)/17280000)+((sqrt(15)*i)/28800000))^(1/3)))))+(2*(((("
     "-1)/17280000)+((sqrt(15)*i)/28800000))^(1/3)))))))+3)"},
    {"mp.issue5934_den", 'Z', NULL, NULL, NULL,
     "(((-36000)-(7200*sqrt(5)))+((((12*sqrt(10))*sqrt((sqrt(5)+5)))+((24*sqrt(10))*sqrt(((-sqrt(5))+5))))"
     "^2))"},
    {"mp.sin7_sqrt2", 'P', NULL, NULL, "28561 0 -268432 0 770912 0 -826496 0 351488 0 -63488 0 4096",
     "(sin((pi/7))+sqrt(2))"},
    {"mp.exp7_sqrt2", 'P', NULL, NULL, "127 142 -37 -212 211 126 -97 -70 43 16 -9 -2 1",
     "(exp(((i*pi)/7))+sqrt(2))"},
    {"mp.cos7_ratio", 'Z', NULL, NULL, NULL,
     "((((5*cos(((2*pi)/7)))-7)/((9*cos((pi/7)))-(5*cos(((3*pi)/7)))))+(1/(2*cos((pi/7)))))"},
    {"mp.cube_sq", 'P', NULL, NULL, "-3008 7424 1984 -5056 480 448 -56 -8 1",
     "(((((1+sqrt(2))-(2*sqrt(3)))+sqrt(7))^3)^(1/3))"},
    {"mp.cube_sq2", 'Z', NULL, NULL, NULL,
     "(((((1+(5*sqrt(2)))+(2*sqrt(3)))^3)^(1/3))-((1+(5*sqrt(2)))+(2*sqrt(3))))"},
    {"mp.mixed48", 'N', "4.3971", "0", NULL,
     "((sqrt((1+(2^(1/3))))+sqrt((1+(2^(1/4)))))+sqrt(2))"},
    {"mp.hi_prec", 'N', "1.588", "0", NULL,
     "(1/sqrt((((1-(9*sqrt(2)))+(7*sqrt(3)))+(1/(10^30)))))"},
    {"mp.cos15", 'P', NULL, NULL, "1 -8 -16 8 16",
     "cos((pi/15))"},
    {"mp.sin11", 'P', NULL, NULL, "-11 0 220 0 -1232 0 2816 0 -2816 0 1024",
     "sin((pi/11))"},
    {"mp.sin21", 'P', NULL, NULL, "1 0 -64 0 960 0 -4992 0 11264 0 -11264 0 4096",
     "sin((pi/21))"},
    {"mp.tan5", 'P', NULL, NULL, "5 0 -10 0 1",
     "tan((pi/5))"},
    {"mp.tan10", 'P', NULL, NULL, "1 0 -10 0 5",
     "tan((pi/10))"},
    {"mp.root_rot", 'P', NULL, NULL, "-2 0 0 1",
     "((2^(1/3))*exp((((2*i)*pi)/3)))"},
    {"mp.cbrt_m1", 'Z', NULL, NULL, NULL,
     "(((-((-1)^(1/3)))+((-1)^(2/3)))+1)"},
    {"mp.exp3ipi", 'Z', NULL, NULL, NULL,
     "(exp(((3*i)*pi))+1)"},
    {"mp.not_alg1", 'N', "-0.2663", "0", NULL,
     "cos((pi*sqrt(2)))"},
    {"mp.not_alg2", 'N', "-0.26626", "-0.96390", NULL,
     "exp(((i*pi)*sqrt(2)))"},
    {"dn.shanks29a", 'Z', NULL, NULL, NULL,
     "((sqrt(((16-(2*sqrt(29)))+(2*sqrt((55-(10*sqrt(29)))))))-sqrt(5))-sqrt((11-(2*sqrt(29)))))"},
    {"dn.shanks29b", 'Z', NULL, NULL, NULL,
     "(sqrt(((-sqrt(5))+sqrt(((((-2)*sqrt(29))+(2*sqrt((((-10)*sqrt(29))+55))))+16))))-((11-(2*sqrt(29)))^"
     "(1/4)))"},
    {"dn.jr43", 'Z', NULL, NULL, NULL,
     "((sqrt(((5*sqrt(3))+(6*sqrt(2))))-(sqrt(2)*(3^(1/4))))-(3^(3/4)))"},
    {"dn.mq1", 'Z', NULL, NULL, NULL,
     "(sqrt((((((-4)*sqrt(14))-(2*sqrt(6)))+(4*sqrt(21)))+33))-(((-sqrt(2))+sqrt(3))+(2*sqrt(7))))"},
    {"dn.mq1_ctrl", 'N', "0.088440", "0", NULL,
     "(sqrt((((((-4)*sqrt(14))-(2*sqrt(6)))+(4*sqrt(21)))+34))-(((-sqrt(2))+sqrt(3))+(2*sqrt(7))))"},
    {"dn.mq2", 'Z', NULL, NULL, NULL,
     "(sqrt((((((-28)*sqrt(7))-(14*sqrt(5)))+(4*sqrt(35)))+82))-(((-7)+sqrt(5))+(2*sqrt(7))))"},
    {"dn.mq3", 'Z', NULL, NULL, NULL,
     "(sqrt(((((468*sqrt(3))+(3024*sqrt(2)))+(2912*sqrt(6)))+19735))-(((9*sqrt(3))+26)+(56*sqrt(6))))"},
    {"dn.mq4", 'Z', NULL, NULL, NULL,
     "(sqrt((((((-490)*sqrt(3))-(98*sqrt(115)))-(98*sqrt(345)))-2107))-(i*(((7*sqrt(5))+(7*sqrt(15)))+(7*s"
     "qrt(23)))))"},
    {"dn.mq5", 'Z', NULL, NULL, NULL,
     "(sqrt(((((4*sqrt(15))+(8*sqrt(5)))+(12*sqrt(3)))+24))-(((1+sqrt(3))+sqrt(5))+sqrt(15)))"},
    {"dn.mq6", 'Z', NULL, NULL, NULL,
     "(sqrt((sqrt(((2*sqrt(6))+5))+sqrt(((2*sqrt(7))+8))))-sqrt((((1+sqrt(2))+sqrt(3))+sqrt(7))))"},
    {"dn.nonga1", 'Z', NULL, NULL, NULL,
     "(sqrt(((13-(2*sqrt(10)))+((2*sqrt(2))*sqrt((((-2)*sqrt(10))+11)))))-(((-1)+sqrt(2))+sqrt(10)))"},
    {"dn.nonga2", 'Z', NULL, NULL, NULL,
     "(sqrt(((112+(70*sqrt(2)))+((46+(34*sqrt(2)))*sqrt(5))))-(((sqrt(10)+5)+(4*sqrt(2)))+(3*sqrt(5))))"},
    {"dn.nonga3", 'Z', NULL, NULL, NULL,
     "(sqrt((((((2*sqrt(2))*sqrt((sqrt(2)+2)))+(5*sqrt(2)))+(4*sqrt((sqrt(2)+2))))+8))-((sqrt(2)+sqrt((sqr"
     "t(2)+2)))+2))"},
    {"dn.c55", 'Z', NULL, NULL, NULL,
     "(sqrt(((8-(sqrt(2)*sqrt((5-sqrt(5)))))-(sqrt(3)*(1+sqrt(5)))))-((((((((((-sqrt(15))*sqrt((5-sqrt(5))"
     "))-(sqrt(3)*sqrt((5-sqrt(5)))))+sqrt((5-sqrt(5))))+(sqrt(5)*sqrt((5-sqrt(5)))))-sqrt(6))-sqrt(2))+sq"
     "rt(10))+sqrt(30))/4))"},
    {"dn.complex1", 'Z', NULL, NULL, NULL,
     "((((3-(sqrt(2)*sqrt((4+(3*i)))))+(3*i))/2)-i)"},
    {"dn.complex2", 'Z', NULL, NULL, NULL,
     "((-sqrt(((-2)+((2*sqrt(3))*i))))-((-1)-(sqrt(3)*i)))"},
    {"dn.complex3", 'Z', NULL, NULL, NULL,
     "(sqrt(((-8)-sqrt(63)))-((i*(sqrt(14)+(3*sqrt(2))))/2))"},
    {"dn.recip", 'Z', NULL, NULL, NULL,
     "(sqrt(((1/((4*sqrt(3))+7))+1))-((sqrt(2)+sqrt(6))/(sqrt(3)+2)))"},
    {"dn.ctrl1", 'N', "2.4511", "0", NULL,
     "sqrt(((15-(2*sqrt(31)))+(2*sqrt((55-(10*sqrt(29)))))))"},
    {"dn.ctrl2", 'P', NULL, NULL, "2 0 -16 0 20 0 -8 0 1",
     "sqrt((2+sqrt((2+sqrt(2)))))"},
    {"dn.ram1", 'Z', NULL, NULL, NULL,
     "((((2^(1/3))-1)^(1/3))-((((1/9)^(1/3))-((2/9)^(1/3)))+((4/9)^(1/3))))"},
    {"dn.ram2", 'Z', NULL, NULL, NULL,
     "(sqrt(((5^(1/3))-(4^(1/3))))-((((2^(1/3))+(20^(1/3)))-(25^(1/3)))/3))"},
    {"dn.ram3", 'Z', NULL, NULL, NULL,
     "(sqrt(((28^(1/3))-(27^(1/3))))-((((98^(1/3))-(28^(1/3)))-1)/3))"},
    {"dn.ram4", 'Z', NULL, NULL, NULL,
     "((((3+(2*(5^(1/4))))/(3-(2*(5^(1/4)))))^(1/4))-(((5^(1/4))+1)/((5^(1/4))-1)))"},
    {"dn.ram6", 'Z', NULL, NULL, NULL,
     "(((((32/5)^(1/5))-((27/5)^(1/5)))^(1/3))-((((1/25)^(1/5))+((3/25)^(1/5)))-((9/25)^(1/5))))"},
    {"dn.ram7", 'Z', NULL, NULL, NULL,
     "((((49+(20*sqrt(6)))^(1/4))+((49-(20*sqrt(6)))^(1/4)))-(2*sqrt(3)))"},
    {"dn.ram_cos9", 'Z', NULL, NULL, NULL,
     "((((sgn(cos(((2*pi)/9)))*abs(cos(((2*pi)/9)))^(1/3))+(sgn(cos(((4*pi)/9)))*abs(cos(((4*pi)/9)))^(1/3"
     ")))+(sgn(cos(((8*pi)/9)))*abs(cos(((8*pi)/9)))^(1/3)))-(sgn((((3*(9^(1/3)))-6)/2))*abs((((3*(9^(1/3)"
     "))-6)/2))^(1/3)))"},
    {"dn.ram_cos7", 'Z', NULL, NULL, NULL,
     "((((sgn(cos(((2*pi)/7)))*abs(cos(((2*pi)/7)))^(1/3))+(sgn(cos(((4*pi)/7)))*abs(cos(((4*pi)/7)))^(1/3"
     ")))+(sgn(cos(((8*pi)/7)))*abs(cos(((8*pi)/7)))^(1/3)))-(sgn(((5-(3*(7^(1/3))))/2))*abs(((5-(3*(7^(1/"
     "3))))/2))^(1/3)))"},
    {"dn.cavallo", 'Z', NULL, NULL, NULL,
     "((sgn((7-(5*sqrt(2))))*abs((7-(5*sqrt(2))))^(1/3))-(1-sqrt(2)))"},
    {"dn.cavallo_principal", 'N', "0.6213", "0.3587", NULL,
     "(((7-(5*sqrt(2)))^(1/3))-(1-sqrt(2)))"},
    {"dn.shanks", 'Z', NULL, NULL, NULL,
     "(((sqrt(5)+sqrt((22+(2*sqrt(5)))))-sqrt((11+(2*sqrt(29)))))-sqrt(((16-(2*sqrt(29)))+(2*sqrt((55-(10*"
     "sqrt(29))))))))"},
    {"dn.jr44", 'Z', NULL, NULL, NULL,
     "(((sqrt((((12+(2*sqrt(6)))+(2*sqrt(14)))+(2*sqrt(21))))-sqrt(2))-sqrt(3))-sqrt(7))"},
    {"dn.bombelli", 'Z', NULL, NULL, NULL,
     "((((2+(11*i))^(1/3))+((2-(11*i))^(1/3)))-4)"},
    {"w.C14", 'Z', NULL, NULL, NULL,
     "((sqrt(((2*sqrt(3))+4))-1)-sqrt(3))"},
    {"w.C15", 'Z', NULL, NULL, NULL,
     "((sqrt((14+(3*sqrt((3+(2*sqrt((5-(12*sqrt((3-(2*sqrt(2)))))))))))))-3)-sqrt(2))"},
    {"w.C16", 'Z', NULL, NULL, NULL,
     "(((sqrt((((10+(2*sqrt(6)))+(2*sqrt(10)))+(2*sqrt(15))))-sqrt(2))-sqrt(3))-sqrt(5))"},
    {"w.C18", 'Z', NULL, NULL, NULL,
     "((sqrt(((-2)+sqrt((-5))))*sqrt(((-2)-sqrt((-5)))))-3)"},
    {"w.C19", 'Z', NULL, NULL, NULL,
     "((((90+(34*sqrt(7)))^(1/3))-3)-sqrt(7))"},
    {"w.C20", 'Z', NULL, NULL, NULL,
     "((((((135+(78*sqrt(3)))^(2/3))+3)*sqrt(3))/((135+(78*sqrt(3)))^(1/3)))-12)"},
    {"w.C21", 'Z', NULL, NULL, NULL,
     "((((41+(29*sqrt(2)))^(1/5))-1)-sqrt(2))"},
    {"w.C22", 'Z', NULL, NULL, NULL,
     "(((((((6-(4*sqrt(2)))*log((3-(2*sqrt(2)))))+((3-(2*sqrt(2)))*log((17-(12*sqrt(2))))))+32)-(24*sqrt(2"
     ")))/((48*sqrt(2))-72))-((sqrt(2)/3)-(log((sqrt(2)-1))/3)))"},
    {"w.C13", 'Z', NULL, NULL, NULL,
     "(((10*((1+(29/1000))^(1/3)))/7)-(3^(1/3)))"},
    {"w.K2", 'Z', NULL, NULL, NULL,
     "(abs(((3-sqrt(7))+(i*sqrt(((6*sqrt(7))-15)))))-1)"},
    {"w.K4", 'Z', NULL, NULL, NULL,
     "((log((3+(4*i)))-log(5))-(i*atan((4/3))))"},
    {"w.L1", 'Z', NULL, NULL, NULL,
     "(sqrt(997)-((997^3)^(1/6)))"},
    {"w.L2", 'Z', NULL, NULL, NULL,
     "(sqrt(999983)-((999983^3)^(1/6)))"},
    {"w.L3", 'Z', NULL, NULL, NULL,
     "(((((2^(1/3))+(4^(1/3)))^3)-(6*((2^(1/3))+(4^(1/3)))))-6)"},
    {"w.I1", 'Z', NULL, NULL, NULL,
     "(tan(((7*pi)/10))+sqrt((1+(2/sqrt(5)))))"},
    {"w.I2", 'Z', NULL, NULL, NULL,
     "(sqrt(((1+cos(6))/2))+cos(3))"},
    {"w.K8", 'N', "0", "2", NULL,
     "(sqrt((1/(-1)))-(1/sqrt((-1))))"},
    {"tr.gauss11", 'Z', NULL, NULL, NULL,
     "((tan(((3*pi)/11))+(4*sin(((2*pi)/11))))-sqrt(11))"},
    {"tr.prod_tan11", 'Z', NULL, NULL, NULL,
     "(((((tan((pi/11))*tan(((2*pi)/11)))*tan(((3*pi)/11)))*tan(((4*pi)/11)))*tan(((5*pi)/11)))-sqrt(11))"},
    {"tr.prod_sin7", 'Z', NULL, NULL, NULL,
     "(((sin((pi/7))*sin(((2*pi)/7)))*sin(((3*pi)/7)))-(sqrt(7)/8))"},
    {"tr.sum_sin7", 'Z', NULL, NULL, NULL,
     "(((sin((pi/7))-sin(((2*pi)/7)))-sin(((4*pi)/7)))+(sqrt(7)/2))"},
    {"tr.alt_cos7", 'Z', NULL, NULL, NULL,
     "(((cos((pi/7))-cos(((2*pi)/7)))+cos(((3*pi)/7)))-(1/2))"},
    {"tr.morrie", 'Z', NULL, NULL, NULL,
     "(((cos((pi/9))*cos(((2*pi)/9)))*cos(((4*pi)/9)))-(1/8))"},
    {"tr.gauss17", 'Z', NULL, NULL, NULL,
     "(cos(((2*pi)/17))-(((((-1)+sqrt(17))+sqrt((34-(2*sqrt(17)))))+(2*sqrt((((17+(3*sqrt(17)))-sqrt((34-("
     "2*sqrt(17)))))-(2*sqrt((34+(2*sqrt(17)))))))))/16))"},
    {"tr.cos153", 'Z', NULL, NULL, NULL,
     "((cos((pi/153))-(sin((pi/9))*sin(((2*pi)/17))))-(cos((pi/9))*cos(((2*pi)/17))))"},
    {"tr.tan2pi5", 'Z', NULL, NULL, NULL,
     "(tan(((2*pi)/5))-(sqrt(((sqrt(5)/8)+(5/8)))/(((-1)/4)+(sqrt(5)/4))))"},
    {"tr.zeta10", 'Z', NULL, NULL, NULL,
     "((exp(((i*pi)/5))-((sqrt(5)+1)/4))-(i*sqrt(((5-sqrt(5))/8))))"},
    /* cos(pi/85) as nested square roots (SymPy's rewrite(sqrt)): a
       difference split into its two parts, decided by a minimal
       polynomial (GR_TOWER_OPT_SPLIT_DEGREE_LIMIT) */
    {"tr.cos85", 'Z', NULL, NULL, NULL,
     "((((((-1)/4)+(sqrt(5)/4))*((((((((((((((((-3)*sqrt((17-sqrt(17))))*sqrt((sqrt(17)+17)))*sqrt(((((sqr"
     "t(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17))"
     ")-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))"
     "/16)-((((3*sqrt(34))*sqrt((sqrt(17)+17)))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+"
     "((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqr"
     "t((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))/32))-((((3*sqrt((17-sqrt(17))))*sqrt(((((sqrt"
     "(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))"
     "-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))*"
     "sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt"
     "(17)))))+(6*sqrt(17)))+34)))/32))-(((sqrt((sqrt(17)+17))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sq"
     "rt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17))))"
     ")+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))*sqrt((((((((-8)*sqrt(2))*sqrt("
     "(sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/1"
     "6))-((9*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))"
     "*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+3"
     "4)))/32))+(15/32))))/8))-(((sqrt(34)*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqr"
     "t(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17"
     "-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqr"
     "t(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))-(((sqrt(34)*sqrt"
     "((17-sqrt(17))))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)"
     "*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqr"
     "t(17)))+34)))/32))+(15/32))))/32))+((((7*sqrt(2))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))"
     "))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt"
     "(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(1"
     "7)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+((("
     "(sqrt(17)*sqrt((17-sqrt(17))))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*s"
     "qrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt("
     "17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*s"
     "qrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+((((11*sqrt(2))*sqrt(("
     "sqrt(17)+17)))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*s"
     "qrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt("
     "17)))+34)))/32))+(15/32))))/32))+(((5*sqrt(17))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))"
     "/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(3"
     "4)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))/8))+((((19*sqrt(2))*sqrt((17-sqrt(17)))"
     ")*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt("
     "(sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/3"
     "2))+(15/32))))/32)))+(sqrt(((sqrt(5)/8)+(5/8)))*sqrt((((((((-9)*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt"
     "((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt"
     "(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))/16)+((sqrt(17)*sqrt(((("
     "(sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+"
     "17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32"
     "))))/16))+(((sqrt(2)*sqrt((17-sqrt(17))))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+"
     "((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqr"
     "t((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))/16))+(((sqrt(2)*sqrt(((((sqrt(17)/32)+((sqrt("
     "2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt(("
     "17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))*sqrt((((((((-8)*"
     "sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt"
     "(17)))+34)))/16))+(1/2))))))-cos(pi/85)"},
    {"tr.cos85_3", 'N', "5.457923e-03", "0", NULL,
     "((((((-1)/4)+(sqrt(5)/4))*((((((((((((((((-3)*sqrt((17-sqrt(17))))*sqrt((sqrt(17)+17)))*sqrt(((((sqr"
     "t(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17))"
     ")-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))"
     "/16)-((((3*sqrt(34))*sqrt((sqrt(17)+17)))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+"
     "((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqr"
     "t((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))/32))-((((3*sqrt((17-sqrt(17))))*sqrt(((((sqrt"
     "(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))"
     "-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))*"
     "sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt"
     "(17)))))+(6*sqrt(17)))+34)))/32))-(((sqrt((sqrt(17)+17))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sq"
     "rt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17))))"
     ")+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))*sqrt((((((((-8)*sqrt(2))*sqrt("
     "(sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/1"
     "6))-((9*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))"
     "*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+3"
     "4)))/32))+(15/32))))/8))-(((sqrt(34)*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqr"
     "t(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17"
     "-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqr"
     "t(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))-(((sqrt(34)*sqrt"
     "((17-sqrt(17))))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)"
     "*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqr"
     "t(17)))+34)))/32))+(15/32))))/32))+((((7*sqrt(2))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))"
     "))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt"
     "(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(1"
     "7)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+((("
     "(sqrt(17)*sqrt((17-sqrt(17))))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*s"
     "qrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt("
     "17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*s"
     "qrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+((((11*sqrt(2))*sqrt(("
     "sqrt(17)+17)))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*s"
     "qrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt("
     "17)))+34)))/32))+(15/32))))/32))+(((5*sqrt(17))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))"
     "/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(3"
     "4)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))/8))+((((19*sqrt(2))*sqrt((17-sqrt(17)))"
     ")*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt("
     "(sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/3"
     "2))+(15/32))))/32)))+(sqrt(((sqrt(5)/8)+(5/8)))*sqrt((((((((-9)*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt"
     "((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt"
     "(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))/16)+((sqrt(17)*sqrt(((("
     "(sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+"
     "17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32"
     "))))/16))+(((sqrt(2)*sqrt((17-sqrt(17))))*sqrt(((((sqrt(17)/32)+((sqrt(2)*sqrt((17-sqrt(17))))/32))+"
     "((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqr"
     "t((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))/16))+(((sqrt(2)*sqrt(((((sqrt(17)/32)+((sqrt("
     "2)*sqrt((17-sqrt(17))))/32))+((sqrt(2)*sqrt((((((((-8)*sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt(("
     "17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt(17)))+34)))/32))+(15/32))))*sqrt((((((((-8)*"
     "sqrt(2))*sqrt((sqrt(17)+17)))-(sqrt(2)*sqrt((17-sqrt(17)))))+(sqrt(34)*sqrt((17-sqrt(17)))))+(6*sqrt"
     "(17)))+34)))/16))+(1/2))))))-cos(3*pi/85)"},
    {"tr.cheb5", 'Z', NULL, NULL, NULL,
     "(cos((5*acos((1/3))))-(241/243))"},
    {"tr.cos_big", 'Z', NULL, NULL, NULL,
     "(cos(((10^30)*pi/3))+(1/2))"},
    {"tr.exp_big", 'Z', NULL, NULL, NULL,
     "(exp(((10^40)*pi*i/7))-(exp((2*pi*i/7))^2))"},
    {"tr.sin_big", 'Z', NULL, NULL, NULL,
     "sin(((10^30)*pi))"},
    {"el.asinh", 'Z', NULL, NULL, NULL,
     "(asinh(1)-log((1+sqrt(2))))"},
    {"el.acosh3", 'Z', NULL, NULL, NULL,
     "(acosh(3)-(2*log((1+sqrt(2)))))"},
    {"el.gs_tower", 'Z', NULL, NULL, NULL,
     "(((sqrt(2)^sqrt(2))^sqrt(2))-2)"},
    {"el.pow_exp", 'Z', NULL, NULL, NULL,
     "(exp((sqrt(2)*log(2)))-(2^sqrt(2)))"},
    {"el.log_fac", 'Z', NULL, NULL, NULL,
     "((log((sqrt((sqrt(2)/3))/2))-log(sqrt((sqrt(2)/3))))-log((1/2)))"},
    {"el.lambertw1", 'Z', NULL, NULL, NULL,
     "(lambertw(((-log(2))/2))+log(2))"},
    {"el.lambertw2", 'Z', NULL, NULL, NULL,
     "(lambertw((2*log(2)))-log(2))"},
    {"el.omega", 'Z', NULL, NULL, NULL,
     "((lambertw(1)*exp(lambertw(1)))-1)"},
    {"nm.h163_3", 'N', "6.05e-10", "0", NULL,
     "(exp(((pi*sqrt(163))/3))-640320)"},
    {"nm.h67", 'N', "-1.33e-6", "0", NULL,
     "(exp((pi*sqrt(67)))-((5280^3)+744))"},
    {"nm.h43", 'N', "-2.22e-4", "0", NULL,
     "(exp((pi*sqrt(43)))-((960^3)+744))"},
    {"nm.r58", 'N', "-1.78e-7", "0", NULL,
     "(exp((pi*sqrt(58)))-((396^4)-104))"},
    {"nm.epi", 'N', "-9.0e-4", "0", NULL,
     "((exp(pi)-pi)-20)"},
    {"nm.pi4", 'N', "2.75e-6", "0", NULL,
     "((22*(pi^4))-2143)"},
    {"nm.pi45", 'N', "-1.767e-5", "0", NULL,
     "(((pi^4)+(pi^5))-exp(6))"},
    {"nm.log163", 'N', "-1.26e-6", "0", NULL,
     "((163/log(163))-32)"},
    {"nm.phi17", 'N', "2.8e-4", "0", NULL,
     "((((1+sqrt(5))/2)^17)-3571)"},
    {"nm.rootred", 'N', "1e-20000", "0", NULL,
     "(((sqrt(2)+sqrt(3))-sqrt((5+(2*sqrt(6)))))+(1/(10^20000)))"},
    {"nm.sq1", 'N', "1.94e-4", "0", NULL,
     "(((sqrt(10)+sqrt(11))-sqrt(5))-sqrt(18))"},
    {"nm.sq2", 'N', "-4.8e-6", "0", NULL,
     "((((sqrt(5)+sqrt(6))+sqrt(18))-sqrt(4))-(2*sqrt(12)))"},
    {"nm.sq3", 'N', "2.84e-20", "0", NULL,
     "(((((sqrt(29)+sqrt(1097))+sqrt(3153))-sqrt(226))-sqrt(2324))-sqrt(987))"},
    {"nm.sq4", 'N', "-1.26e-15", "0", NULL,
     "(((sqrt(11075)+sqrt(27187))+sqrt(68057))-531)"},
    {"nm.sq5", 'N', "1.83e-5", "0", NULL,
     "(((sqrt(3)+sqrt(20))+sqrt(23))-11)"},
    {"sage.34gon", 'Z', NULL, NULL, NULL,
     "(((sqrt(2)*sqrt(((15+sqrt(17))+(sqrt(2)*(sqrt((((34+(6*sqrt(17)))+((sqrt(2)*(sqrt(17)-1))*sqrt((17-s"
     "qrt(17)))))-((8*sqrt(2))*sqrt((17+sqrt(17))))))+sqrt((17-sqrt(17))))))))/8)-cos((pi/17)))"},
    {"sage.cardano", 'Z', NULL, NULL, NULL,
     "((((((2/(3*sqrt(3)))+(10/27))^(1/3))-(2/(9*(((2/(3*sqrt(3)))+(10/27))^(1/3)))))+(1/3))-1)"},
    /* special functions (hg, mf: elliptic integrals, modular functions;
       rich.torsion: exponentials related modulo roots of unity) */
    {"hg.K_half", 'Z', NULL, NULL, NULL,
     "(elliptic_k((1/2))-((gamma((1/4))^2)/(4*sqrt(pi))))"},
    {"hg.K_imag", 'Z', NULL, NULL, NULL,
     "(elliptic_k((-3))-(elliptic_k((3/4))/2))"},
    {"hg.legendre", 'Z', NULL, NULL, NULL,
     "((((elliptic_e((1/5))*elliptic_k((4/5)))+(elliptic_e((4/5))*elliptic_k((1/5))))-(elliptic_k((1/5))*e"
     "lliptic_k((4/5))))-(pi/2))"},
    {"hg.landen", 'Z', NULL, NULL, NULL,
     "(elliptic_k(((4*sqrt((1/3)))/((1+sqrt((1/3)))^2)))-((1+sqrt((1/3)))*elliptic_k((1/3))))"},
    {"hg.K_singular3", 'Z', NULL, NULL, NULL,
     "(elliptic_k(((2-sqrt(3))/4))-(((3^(1/4))*(gamma((1/3))^3))/((2^(7/3))*pi)))"},
    {"mf.j_i", 'Z', NULL, NULL, NULL,
     "(modular_j(i)-1728)"},
    {"mf.j_rho_tau", 'Z', NULL, NULL, NULL,
     "modular_j((((2*((1+sqrt((-3)))/2))+1)/((3*((1+sqrt((-3)))/2))+2)))"},
    {"mf.j_163", 'Z', NULL, NULL, NULL,
     "(modular_j(((1+sqrt((-163)))/2))+(640320^3))"},
    {"mf.j_163_near", 'N', "-7.4993e-13", "0", NULL,
     "((modular_j(((1+sqrt((-163)))/2))+exp((pi*sqrt(163))))-744)"},
    {"mf.j_sqrtm5", 'Z', NULL, NULL, NULL,
     "((modular_j(sqrt((-5)))-632000)-(282880*sqrt(5)))"},
    {"mf.j_sqrtm14", 'P', NULL, NULL, "10064086044321563803648 2257767342088912896 2059647197077504 -16220384512 1",
     "modular_j(sqrt((-14)))"},
    {"mf.lambda_sqrtm2", 'Z', NULL, NULL, NULL,
     "(modular_lambda(sqrt((-2)))-((sqrt(2)-1)^2))"},
    {"mf.eta_i", 'Z', NULL, NULL, NULL,
     "(dedekind_eta(i)-(gamma((1/4))/(2*(pi^(3/4)))))"},
    {"mf.eta_rho", 'Z', NULL, NULL, NULL,
     "((dedekind_eta(((1+sqrt((-3)))/2))^24)+((27*(gamma((1/3))^36))/((2^24)*(pi^24))))"},
    {"mf.eta_S", 'Z', NULL, NULL, NULL,
     "(dedekind_eta(((-1)/((1/3)+((pi*i)/4))))-(sqrt(((-i)*((1/3)+((pi*i)/4))))*dedekind_eta(((1/3)+((pi*i"
     ")/4)))))"},
    {"mf.eta_gamma", 'Z', NULL, NULL, NULL,
     "((dedekind_eta((((2*((1/3)+((pi*i)/4)))+1)/((7*((1/3)+((pi*i)/4)))+4)))^24)-((((7*((1/3)+((pi*i)/4))"
     ")+4)^12)*(dedekind_eta(((1/3)+((pi*i)/4)))^24)))"},
    {"mf.phi2", 'Z', NULL, NULL, NULL,
     "((((((((modular_j(((1/5)+((pi*i)/3)))^3)+(modular_j((2*((1/5)+((pi*i)/3))))^3))-((modular_j(((1/5)+("
     "(pi*i)/3)))^2)*(modular_j((2*((1/5)+((pi*i)/3))))^2)))+(1488*(((modular_j(((1/5)+((pi*i)/3)))^2)*mod"
     "ular_j((2*((1/5)+((pi*i)/3)))))+(modular_j(((1/5)+((pi*i)/3)))*(modular_j((2*((1/5)+((pi*i)/3))))^2)"
     "))))-(162000*((modular_j(((1/5)+((pi*i)/3)))^2)+(modular_j((2*((1/5)+((pi*i)/3))))^2))))+((40773375*"
     "modular_j(((1/5)+((pi*i)/3))))*modular_j((2*((1/5)+((pi*i)/3))))))+(8748000000*(modular_j(((1/5)+((p"
     "i*i)/3)))+modular_j((2*((1/5)+((pi*i)/3)))))))-157464000000000)"},
    {"rich.torsion2", 'Z', NULL, NULL, NULL,
     "((exp(((16+((30*pi)*i))/225))^75)-((exp(((5+((42*pi)*i))/60))^64)*exp(((((-174)*pi)*i)/5))))"},
    {"rich.torsion2_pi", 'Z', NULL, NULL, NULL,
     "((exp(((pi*(16+(30*i)))/225))^75)-((exp(((pi*(5+(42*i)))/60))^64)*exp(((((-174)*pi)*i)/5))))"},
    {"rich.torsion2_sqrt2", 'Z', NULL, NULL, NULL,
     "((exp((((16*sqrt(2))+((30*pi)*i))/225))^75)-((exp((((5*sqrt(2))+((42*pi)*i))/60))^64)*exp(((((-174)*"
     "pi)*i)/5))))"},
    {"rich.torsion3", 'Z', NULL, NULL, NULL,
     "((((exp(((16+((30*pi)*i))/225))^2)*exp(((1/12)+(((7*pi)*i)/10))))*exp(((3/40)+(((5*pi)*i)/8))))-exp("
     "((((32/225)+(1/12))+(3/40))+((pi*i)*(((4/15)+(7/10))+(5/8))))))"},
    {"rich.torsion_near", 'N', "1e-40", "0", NULL,
     "(((exp(((16+((30*pi)*i))/225))^75)-((exp(((5+((42*pi)*i))/60))^64)*exp(((((-174)*pi)*i)/5))))+(1/(10"
     "^40)))"},
};

#define NUM_CATALOG_CASES (sizeof(catalog_cases) / sizeof(catalog_case_struct))

/* in the shared field: the cases run before the current one */
static const slong * _cat_order = NULL;
static slong _cat_order_len = 0;

static void
_cat_fail(const catalog_case_struct * c, const char * what, gr_srcptr x, gr_ctx_t K)
{
    slong i;

    flint_printf("FAIL: catalog case %s: %s\n", c->name, what);
    if (_cat_order != NULL)
    {
        flint_printf("in a shared field, after:");
        for (i = 0; i < _cat_order_len; i++)
            flint_printf(" %s", catalog_cases[_cat_order[i]].name);
        flint_printf("\n");
    }
    flint_printf("expr = %s\n", c->expr);
    if (x != NULL)
    {
        flint_printf("x = ");
        gr_println(x, K);
    }
    flint_abort();
}

/* whether the value of x is within 5% of re + im i */
static int
_cat_close(gr_srcptr x, const char * re, const char * im, gr_ctx_t K)
{
    acb_t z, a;
    arb_t d, e;
    slong prec;
    int ok;

    acb_init(z); acb_init(a); arb_init(d); arb_init(e);

    /* the value may be tiny (10^-20000, say): increase the precision
       until it is known to a few bits */
    for (prec = 64; ; prec *= 2)
    {
        ok = (gr_tower_lazy_get_acb(z, x, prec, K) == GR_SUCCESS);
        if (!ok || acb_rel_accuracy_bits(z) >= 16 || prec > 400000)
            break;
    }

    ok = ok && arb_set_str(acb_realref(a), re, 64) == 0 &&
         arb_set_str(acb_imagref(a), im, 64) == 0;

    if (ok)
    {
        acb_sub(z, z, a, 64);
        acb_abs(d, z, 64);
        acb_abs(e, a, 64);
        arb_mul_2exp_si(e, e, -4);   /* 1/16 < 1/20 + rounding of the data */
        ok = arb_le(d, e);
    }

    acb_clear(z); acb_clear(a); arb_clear(d); arb_clear(e);
    return ok;
}

static void
_cat_check(const catalog_case_struct * c, gr_ctx_t K)
{
    gr_ptr x, y;
    truth_t t;

    GR_TMP_INIT2(x, y, K);

    if (gr_set_str(x, c->expr, K) != GR_SUCCESS)
        _cat_fail(c, "evaluation", NULL, K);

    if (c->kind == 'Z')
    {
        t = gr_is_zero(x, K);
        if (t != T_TRUE)
            _cat_fail(c, (t == T_FALSE) ? "zero decided nonzero" : "zero not decided", x, K);
    }
    else if (c->kind == 'N')
    {
        t = gr_is_zero(x, K);
        if (t != T_FALSE)
            _cat_fail(c, (t == T_TRUE) ? "nonzero decided zero" : "nonzero not decided", x, K);
        if (!_cat_close(x, c->re, c->im, K))
            _cat_fail(c, "value", x, K);
    }
    else
    {
        /* Horner's rule; the coefficients are separated by spaces */
        const char * s = c->poly;
        const char * end;
        slong n = 0, i;
        fmpz * coeffs;
        char buf[1024];

        for (end = s; *end; end++)
            n += (*end == ' ');
        n++;
        coeffs = _fmpz_vec_init(n);
        for (i = 0; i < n; i++)
        {
            end = strchr(s, ' ');
            if (end == NULL)
                end = s + strlen(s);
            if (end - s >= (slong) sizeof(buf))
                _cat_fail(c, "coefficient too long", NULL, K);
            memcpy(buf, s, end - s);
            buf[end - s] = '\0';
            if (fmpz_set_str(coeffs + i, buf, 10) != 0)
                _cat_fail(c, "coefficient", NULL, K);
            s = end + (*end == ' ');
        }

        GR_MUST_SUCCEED(gr_set_fmpz(y, coeffs + n - 1, K));
        for (i = n - 2; i >= 0; i--)
        {
            GR_MUST_SUCCEED(gr_mul(y, y, x, K));
            GR_MUST_SUCCEED(gr_add_fmpz(y, y, coeffs + i, K));
        }

        t = gr_is_zero(y, K);
        if (t != T_TRUE)
            _cat_fail(c, (t == T_FALSE) ? "P(x) decided nonzero" : "P(x) = 0 not decided", x, K);

        _fmpz_vec_clear(coeffs, n);
    }

    GR_TMP_CLEAR2(x, y, K);
}

/*
    Sequences in a shared field which went wrong (found by the shared
    runs): each entry is either the name of a case, which is checked, or
    an expression, which is evaluated and tested for zero (the generators
    it leaves behind matter, not the answer).
*/
static const char * catalog_sequences[][6] =
{
    /* exp(i) around: no relation exp(260515 i) = exp(i)^260515 */
    { "sin(3+pi)+sin(3)", "sp.atan_tan_pole", NULL },
    /* root_3(70 zeta_3 - 19) and root_3(3) generate the same field: no
       modular proof for the steps above a dynamic one */
    { "sin(pi/21)", "ca.gosper", NULL },
    { "tr.cos153", "ca.gosper", NULL },
    /* sqrt(2) as 2^(1/16)^8: rationalization by powers of the generator */
    { "2^(1/16)", "tr.cheb5", NULL },
    /* (2 sqrt 2 + i)/3 in a tower with a logarithm: log of a root of unity */
    { "cos(5*acos(1/3))", "ca.log_mi", NULL },
    /* exp((pi i + 2 log 5 - 4 log(2 + i))/4) with sqrt 5, sqrt 2 around */
    { "nm.sq2", "sp.nsimplify_cosatan", NULL },
    /* a conjectural logarithm in the tower, not involved */
    { "sp.trig90", "ca.mixed3", "ca.atan_alg", "nm.rootred", NULL },
    /* sqrt of a negative real number with a complex representation */
    { "sp.trig90", "nm.sq5", "dn.mq4", NULL },
    /* a primitive exponential created by the relation search, unnamed
       (a clash of names in the parser) */
    { "el.omega", "sp.ceiling_pyth", "mp.not_alg1", "ca.ramanujan_163", "nm.h163_3", "w.C16" },
};

static void
_cat_run_sequence(const char * const * seq, slong len, gr_ctx_t K)
{
    gr_ptr x;
    slong i, j;

    GR_TMP_INIT(x, K);
    for (i = 0; i < len && seq[i] != NULL; i++)
    {
        for (j = 0; j < (slong) NUM_CATALOG_CASES; j++)
            if (strcmp(catalog_cases[j].name, seq[i]) == 0)
                break;

        if (j < (slong) NUM_CATALOG_CASES)
            _cat_check(catalog_cases + j, K);
        else
        {
            if (gr_set_str(x, seq[i], K) != GR_SUCCESS)
            {
                flint_printf("FAIL: evaluation of %s\n", seq[i]);
                flint_abort();
            }
            (void) gr_is_zero(x, K);
        }
    }
    GR_TMP_CLEAR(x, K);
}

static void
_cat_expect_zero(gr_srcptr x, const char * what, gr_ctx_t K)
{
    truth_t t = gr_is_zero(x, K);

    if (t != T_TRUE)
    {
        flint_printf("FAIL: %s (%d)\n", what, t);
        flint_printf("x = "); gr_println(x, K);
        flint_abort();
    }
}

static void
_cat_set_str(gr_ptr x, const char * s, gr_ctx_t K)
{
    if (gr_set_str(x, s, K) != GR_SUCCESS)
    {
        flint_printf("FAIL: evaluation of %s\n", s);
        flint_abort();
    }
}

/* sum_{k=1}^{p-1} (k|p) cos(2 pi k/p) = sqrt(p) for a prime p = 1 mod 4 */
static void
_cat_gauss_sum(ulong p, gr_ctx_t K)
{
    gr_ptr s, t;
    ulong k;
    char buf[64];

    GR_TMP_INIT2(s, t, K);
    for (k = 1; k < p; k++)
    {
        flint_sprintf(buf, "cos(2*pi*%wu/%wu)", k, p);
        _cat_set_str(t, buf, K);
        if (n_jacobi(k, p) == 1)
            GR_MUST_SUCCEED(gr_add(s, s, t, K));
        else
            GR_MUST_SUCCEED(gr_sub(s, s, t, K));
    }
    GR_MUST_SUCCEED(gr_set_ui(t, p, K));
    GR_MUST_SUCCEED(gr_sqrt(t, t, K));
    GR_MUST_SUCCEED(gr_sub(s, s, t, K));
    flint_sprintf(buf, "Gauss sum, p = %wu", p);
    _cat_expect_zero(s, buf, K);
    GR_TMP_CLEAR2(s, t, K);
}

/* prod_{k=1}^{n-1} 2 sin(k pi/n) = n */
static void
_cat_sine_product(ulong n, gr_ctx_t K)
{
    gr_ptr s, t;
    ulong k;
    char buf[64];

    GR_TMP_INIT2(s, t, K);
    GR_MUST_SUCCEED(gr_one(s, K));
    for (k = 1; k < n; k++)
    {
        flint_sprintf(buf, "2*sin(%wu*pi/%wu)", k, n);
        _cat_set_str(t, buf, K);
        GR_MUST_SUCCEED(gr_mul(s, s, t, K));
    }
    GR_MUST_SUCCEED(gr_sub_ui(s, s, n, K));
    flint_sprintf(buf, "sine product, n = %wu", n);
    _cat_expect_zero(s, buf, K);
    GR_TMP_CLEAR2(s, t, K);
}

/* S_n(sqrt(2) + sqrt(3) + ... + sqrt(p_n)) = 0 for the Swinnerton-Dyer
   polynomial S_n (degree 2^n) */
static void
_cat_swinnerton_dyer(ulong n, gr_ctx_t K)
{
    fmpz_poly_t S;
    gr_ptr a, y, t;
    slong i;
    ulong p;
    char buf[64];

    fmpz_poly_init(S);
    fmpz_poly_swinnerton_dyer(S, n);
    GR_TMP_INIT3(a, y, t, K);

    for (i = 0, p = 2; i < (slong) n; i++, p = n_nextprime(p, 1))
    {
        GR_MUST_SUCCEED(gr_set_ui(t, p, K));
        GR_MUST_SUCCEED(gr_sqrt(t, t, K));
        GR_MUST_SUCCEED(gr_add(a, a, t, K));
    }

    GR_MUST_SUCCEED(gr_set_fmpz(y, S->coeffs + S->length - 1, K));
    for (i = S->length - 2; i >= 0; i--)
    {
        GR_MUST_SUCCEED(gr_mul(y, y, a, K));
        GR_MUST_SUCCEED(gr_add_fmpz(y, y, S->coeffs + i, K));
    }

    flint_sprintf(buf, "Swinnerton-Dyer, n = %wu", n);
    _cat_expect_zero(y, buf, K);
    GR_TMP_CLEAR3(a, y, t, K);
    fmpz_poly_clear(S);
}

/* x = IDFT(DFT(x)) for the inputs of the exact DFT benchmark
   (https://fredrikj.net/blog/2020/09/benchmarking-exact-dft-computation/) */
static void
_cat_dft(slong N, int kind, gr_ctx_t K)
{
    static const char * inputs[] = { "%wd", "sqrt(%wd)", "log(%wd)",
        "exp(2*pi*i/%wd)", "1/(1+%wd*pi)", "1/(1+sqrt(%wd)*pi)" };
    gr_ptr x, X, y, w, t;
    slong j, k, sz = K->sizeof_elem;
    char buf[64];

    GR_TMP_INIT_VEC(x, N, K);
    GR_TMP_INIT_VEC(X, N, K);
    GR_TMP_INIT_VEC(y, N, K);
    GR_TMP_INIT2(w, t, K);

    for (j = 0; j < N; j++)
    {
        flint_sprintf(buf, inputs[kind], j + 2);
        _cat_set_str(GR_ENTRY(x, j, sz), buf, K);
    }

    /* X_k = sum_j x_j w^(jk), w = exp(-2 pi i/N) */
    flint_sprintf(buf, "exp(-2*pi*i/%wd)", N);
    _cat_set_str(w, buf, K);
    for (k = 0; k < N; k++)
        for (j = 0; j < N; j++)
        {
            GR_MUST_SUCCEED(gr_pow_ui(t, w, (j * k) % N, K));
            GR_MUST_SUCCEED(gr_mul(t, t, GR_ENTRY(x, j, sz), K));
            GR_MUST_SUCCEED(gr_add(GR_ENTRY(X, k, sz), GR_ENTRY(X, k, sz), t, K));
        }

    /* y_j = (1/N) sum_k X_k w^(-jk) */
    for (j = 0; j < N; j++)
    {
        for (k = 0; k < N; k++)
        {
            GR_MUST_SUCCEED(gr_pow_ui(t, w, (N - (j * k) % N) % N, K));
            GR_MUST_SUCCEED(gr_mul(t, t, GR_ENTRY(X, k, sz), K));
            GR_MUST_SUCCEED(gr_add(GR_ENTRY(y, j, sz), GR_ENTRY(y, j, sz), t, K));
        }
        GR_MUST_SUCCEED(gr_div_si(GR_ENTRY(y, j, sz), GR_ENTRY(y, j, sz), N, K));
        GR_MUST_SUCCEED(gr_sub(t, GR_ENTRY(y, j, sz), GR_ENTRY(x, j, sz), K));
        flint_sprintf(buf, "DFT, N = %wd, input %d, entry %wd", N, kind, j);
        _cat_expect_zero(t, buf, K);
    }

    GR_TMP_CLEAR_VEC(x, N, K);
    GR_TMP_CLEAR_VEC(X, N, K);
    GR_TMP_CLEAR_VEC(y, N, K);
    GR_TMP_CLEAR2(w, t, K);
}

static int
_cat_is_special(const catalog_case_struct * c)
{
    return strncmp(c->name, "hg.", 3) == 0 || strncmp(c->name, "mf.", 3) == 0 ||
           strncmp(c->name, "rich.torsion", 12) == 0;
}

TEST_FUNCTION_START(gr_tower_catalog, state)
{
    gr_ctx_t QQ, K;
    slong i, j, iter;
    slong perm[NUM_CATALOG_CASES];

    gr_ctx_init_fmpq(QQ);

    /* each case in a fresh field */
    for (i = 0; i < (slong) NUM_CATALOG_CASES; i++)
    {
        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        _cat_check(catalog_cases + i, K);
        gr_ctx_clear(K);
    }

    /* sequences which went wrong */
    for (i = 0; i < (slong) (sizeof(catalog_sequences) / sizeof(catalog_sequences[0])); i++)
    {
        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        _cat_run_sequence(catalog_sequences[i], 6, K);
        gr_ctx_clear(K);
    }

    /* all cases in a shared field, in a random order: the algebraic and
       elementary cases together, and the special function cases (hg, mf,
       rich.torsion) together. (Mixed, the cyclotomic fields of the
       latter absorb the radicals of the former: the near miss nm.rootred,
       which needs 66000 bits when sqrt(5 + 2 sqrt(6)) is not denested,
       then takes minutes in a field of degree 128.) Each pass runs the
       whole catalog: their number grows slowly with the multiplier. */
    for (iter = 0; iter < 2 + flint_test_multiplier() / 2; iter++)
    {
        slong n = 0, k;
        int special = iter % 2;

        for (i = 0; i < (slong) NUM_CATALOG_CASES; i++)
            if (_cat_is_special(catalog_cases + i) == special)
                perm[n++] = i;
        for (k = n - 1; k > 0; k--)
        {
            slong r = n_randint(state, k + 1);
            FLINT_SWAP(slong, perm[k], perm[r]);
        }

        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        _cat_order = perm;
        for (i = 0; i < n; i++)
        {
            _cat_order_len = i;
            _cat_check(catalog_cases + perm[i], K);
        }
        _cat_order = NULL;
        gr_ctx_clear(K);
    }

    /* scalable families */
    {
        static const ulong gauss_p[] = { 5, 13, 17, 29, 37, 41, 53, 61 };
        static const ulong sine_n[] = { 7, 12, 30, 60 };

        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        for (i = 0; i < 8; i++)
            _cat_gauss_sum(gauss_p[i], K);
        for (i = 0; i < 4; i++)
            _cat_sine_product(sine_n[i], K);
        gr_ctx_clear(K);

        for (i = 2; i <= 5; i++)
        {
            gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
            _cat_swinnerton_dyer(i, K);
            gr_ctx_clear(K);
        }

        for (j = 0; j < 6; j++)
        {
            gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
            for (i = 2; i <= 8; i++)
                _cat_dft(i, j, K);
            gr_ctx_clear(K);
        }
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
