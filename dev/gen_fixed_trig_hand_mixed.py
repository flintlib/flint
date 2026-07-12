#!/usr/bin/env python3
"""Emit mixed-precision hand-written series for the odd families used by
the bitwise log, atan and sin/cos, in the shape of exp_rs_opt_hand.inc.

All four are series in z = x^2, so an error at Horner level k is
attenuated by TWICE as many powers of x as in the exponential:

    atanh(x) = x + x^3 V_1,   V_1 = 1/3 + z(1/5 + z(1/7 + ...))
    atan(x)  = x - x^3 V_1,   V_1 = 1/3 - z(1/5 - z(1/7 - ...))
    sin(x)   = x - x^3 V_1,   V_1 = 1/3! - z(1/5! - z(1/7! - ...))
    cos(x)   = 1 - z  V_0,    V_0 = 1/2! - z(1/4! - z(1/6! - ...))

so that err(atanh) = x^3 z^(k-1) err(V_k), giving

    bits(V_k) = 64 n - 6 - r (2k + 1)          (odd families)
    bits(V_k) = 64 n - 6 - 2 r (k + 1)         (cos)

Every V_k is positive and decreasing (z < 2^-2r swamps the alternating
term), so the recurrences need no sign handling.  Levels of one and two
limbs stay in registers; wider ones use flint_mpn_mulhigh_n.
"""
import sys
from fractions import Fraction
from math import factorial, log2

FUNCS = ("atanh", "atan", "sin", "cos", "tan")


# tangent numbers T_m: tan y = sum_{m>=1} T_m y^(2m-1).  They have no
# tidy hypergeometric recurrence, which is why the tan series is a poor
# choice at high precision -- but at these sizes only a handful are
# needed and they are simply tabulated.
_TAN = [
        Fraction(1, 1),
        Fraction(1, 3),
        Fraction(2, 15),
        Fraction(17, 315),
        Fraction(62, 2835),
        Fraction(1382, 155925),
        Fraction(21844, 6081075),
        Fraction(929569, 638512875),
        Fraction(6404582, 10854718875),
        Fraction(443861162, 1856156927625),
        Fraction(18888466084, 194896477400625),
        Fraction(113927491862, 2900518163668125),
        Fraction(58870668456604, 3698160658676859375),
        Fraction(8374643517010684, 1298054391195577640625),
        Fraction(689005380505609448, 263505041412702261046875),
        Fraction(129848163681107301953, 122529844256906551386796875),
        Fraction(1736640792209901647222, 4043484860477916195764296875),
        Fraction(418781231495293038913922, 2405873491984360136479756640625),
        Fraction(56518638202982204522669764, 801155872830791925447758961328125),
        Fraction(32207686319158956594455462, 1126482925555250126673224649609375),
        Fraction(1410211493828985228276049834684, 121699582862361447435141825020548828125),
        Fraction(516098083439704913515348955653804, 109894723324712387033933067993555591796875),
        Fraction(103537504005512749467288942622106408, 54397888045732631581796868656810017939453125),
        Fraction(45361105584983995647044252937847808918, 58804116977436974739922415018011629392548828125),
        Fraction(87176517890549500795745183943750553204, 278845328893007589895761129278958371635634765625),
        Fraction(1396470103398938597044980843456514101088564, 11021361624496124990629958634750829638898464111328125),
        Fraction(389951962465960362323362101491789115193414088, 7593718159277830118544041499343321621201041772705078125),
        Fraction(321055735622680218266276441690024211623548948, 15426363155304483893348702635464887287939188826904296875),
        Fraction(37951675284166717133668639194471627545621910540728, 4499391915144503512689122748983408210385935265954349365234375),
        Fraction(641885182338872430017276041951405742741480681339128, 187767306507615744151489976183185645072447200976777847900390625),
        Fraction(9759387159544076997817707959584835600439155088202173136, 7044090503633204641843146456512209474892856744643820963983154296875),
        Fraction(7724760729208487305545342963324697288405380586579904269441, 13757108753595648665519665029568345104465749222289382342659100341796875),
        Fraction(203497294113685566585581532155905318177648366860283186794902, 894212068983717163258778226921942431790273699448809852272841522216796875),
        Fraction(4928399482219367719407723316603414444925848984020131054275466, 53435213095216179674734017830389586937521490526522123875006827178955078125),
        Fraction(65033291600604926267730204296537787351363679128255646162810228, 1739791210461723491420203381737988344092550795455053787171171272613525390625),
        Fraction(5135746785881293900665825337251063912099812018333491820410850621042, 339003946094735963173345929589646769396958527805814290898437688022862701416015625),
        Fraction(23247600823869669181617874661621842533234313612312895049759227683259644, 3786335073932105972683100687586764767394629797063139815044650537527353512115478515625),
        Fraction(26145766198741584025528698683516199629583197662307446174936102767991445644, 10507079830161594074195604408053272229520097686850212986748905241638405996120452880859375),
        Fraction(15502650114137077692879282322945671281504464245814336939030473090718667221928, 15371857791526412130548169248981937271787902915861861599613648368516987972324222564697265625),
        Fraction(624447880395344915327575701011165339822237764093445803722686391628220033123558, 1527764317925576637878029337293978991431565447863561148013214536238736772346159023284912109375),
        Fraction(3177409273870478888667047675148588648707554958779304377142474306979107045033591884, 19181081011555614688558658329725906237423304197927010213305908502477340176806026537342071533203125),
        Fraction(4382231878630427838203781834402677719332895725070795917855918399166899403692359715604, 65273218682323756785165114296057258925951504185545615755880006633930388621670908306575069427490234375),
        Fraction(3170252255497465850441721151634611336636763918620746439504134202802452839465833038132808, 116512695347947905861519729018462207182823434971198924124245811841565743689682571327236498928070068359375),
        Fraction(282743645351878196427175381372737723603024007898834815645760821235876739738607970302830604, 25639646664510183283996782721062771592408380601603245596988446005841026302535441137364220146465301513671875),
        Fraction(26125334033648605299760215770950877021992145130248999739712305337867298482114198536057875576, 5845488211471821649254225408584215172773324360992915294118886395550852065110922559577434464350986480712890625),
        Fraction(3165183288800001305552844563646295861445487620103190581092424979876576485235244926077451042788728, 1747421018496329004719811872515122362672993717853417133447429304653993962083933635347280371600762143611907958984375),
        Fraction(2743910203329295441771249659819135452881820248213784360868316092063767619714125477901042867863808976, 3737733558563647741095677595309846733757533562488459248444051282654893084897534046007832714854030225185871124267578125),
        Fraction(4965369860827668851290623237994135971062634507133697075181739578164128821938437552932557133865324691926, 16688980338986687163992200463058465666227387356510970544302688977054097624067489515424973071823244955454914569854736328125),
        Fraction(13618722892337243626196029509843989171050070079708515121451737330953176267769526553238319205989401876548, 112941704154537813133063496156977058345864412110341684381211220751691683921014870906713189858152657721799538135528564453125),
        Fraction(905838570048586218745173742117616558174626778700773083971608582082083300800057692087180696588351163326044, 18535679696858777383843519947971924100345314960922504303800151196111426769580058982725432267570131653944398944377899169921875),
        Fraction(19969444131507966133609609358428263100759985484071470685128036140914404082310201166506827762184000170499728152, 1008238872188753437093803601373349096900605265613806271179292037794923714338034640538748502381459902774786625800311565399169921875),
        Fraction(3987814874802810721470124591769744923850042632192202120449412189599538513629892760456398486843262214286476453252, 496789266996571997612354958287500740151214448578521679185666841649709573948559717071943890457199878040192837161229193210601806640625),
        Fraction(1387645500322043443555953489772401338297804785107018959194246870273556387596028040241465510555069360983643762020966008, 426535812804251768570003781861273822982730457332790035923625765237815894548624144682214945188194957287138967094073966852724552154541015625),
        Fraction(43689348130876350871102057334620023477738061122122001995347183306853388813969411635784979517513701885490729169745889816, 33135405402916599720006732149798408905959786623756880735929886502241834766921198965655355536469227435278973731376622822216451168060302734375),
        Fraction(1378787166833423254782605792405958387519190426245239460019055283860744326482790137539609993625992339841019464231605232944, 2580202015714823982329295515886386506323847262943568140785226852575256810253926591004868042080287188462449959215375958862689435482025146484375),
        Fraction(218890237082106444667594780636101674529884218866910738369780428071950740720724222345527066614706066991430244099517042667958292, 1010700832350830282256745191499022068944987314502471255234963293269478843223023536076180018977849370113233164320610045614318319325149059295654296875),
        Fraction(3017390902943153954924646545974094232016457803239196621061774645271945050351139675598054703144762373636903084176538255589895569688, 34376967410748790390398674198456237631025853528172554804306806493974783894544699532559110985493590625661399618036909481479808995206294953823089599609375),
        Fraction(8016125698658910896026338407135505037983056097657327298061786976140893111868803952754263980667101276996408494081931256051026542566888, 225341021377458321009063309370880637671374469877171096742231116568004708428740505435924972509910486551210474496231941651100147963577263422310352325439453125),
        Fraction(11023223784956415765092452292756708395157132122068077128614030905176550757470755157331736772521746819562138390511717012267104100078048656, 764582085533716083183751808695398003618973576293241531246390178515239975698716534944093431726126280868257139965714978022182802040417654791899025440216064453125),
        Fraction(1845095806549047480007724487024229333095613891915687099660442773288389958102660070705854697790626890590455116433817177131329560008520762488, 315772401325424742354889496991199375494636087009108752404759143726794109963569928931910587302890153998590198805840285923161497242692491429054297506809234619140625),
        Fraction(23073034661769080715240040747769917624661098665012841744677519691544270422033599652214622890947230491842854950787523478548053096878670074836336, 9743157442895980425360115429663456730886996464666050555448843379690232262925950157194101171230675701626500584154202022159147997423276823043470349572598934173583984375),
        Fraction(1047189088219697125924090306823548244526397427545278624872533528017788773666576764953340621837052680416448057356215085373869618547517384956119248, 1091088213344008076589208150280073370923061708572975780858696595191280786100498567603393150563339698347815431088193698093434140666669343332763552729003131389617919921875),
        Fraction(14654526103938947853592798037203550117196617305224693478424136546593326352491371399756954287697926193752853204535120636684297175083289795800197728, 37674395184125295502656359358487039357082044964524383176677630473006553277209993431843456139780166136926663292495317045818305462240195273675852823150344192981719970703125),
        Fraction(357302767470032900576643605538835088084055212588960920085261795996340330997333306469144562500392344758421560010463942134842407723273904635849262137252097, 2266473492892894213172669801297728855456125068602158955680618267315816729004826866852861604666771560537596201951474583228930218003442810130412667755516771334223449230194091796875),
        Fraction(602593777047525947875926259058537881757401870925833129781936737837405766159506515572089830001168920042332342928487661101418141142568864762750484783869802, 9431454212360753338686271108626032979156133349989629202670959886572269614245892445936101516193984880946771291991620039888128971691745887316878520660053661358542740345001220703125),
]


def coeff(func, k):
    """the coefficient of z^k in V (as a Fraction), all positive"""
    if func == "tan":
        return _TAN[k + 1]                     # V_0 starts at T_2 = 1/3
    if func in ("atanh", "atan"):
        return Fraction(1, 2 * k + 3)          # V_1 starts at 1/3
    if func == "sin":
        return Fraction(1, factorial(2 * k + 3))
    return Fraction(1, factorial(2 * k + 2))   # cos: V_0 starts at 1/2!


def minN(func, n, r):
    """number of V levels needed"""
    k = 0
    while True:
        c = coeff(func, k)
        # the term contributes x^3 z^k (odd) or z^(k+1) (cos)
        if func == "cos":
            mag = -2 * r * (k + 1) + log2(float(c))
        else:
            mag = -r * (2 * k + 3) + log2(float(c))
        if mag < -64 * n:
            return max(1, k)
        k += 1


def width(func, n, r, k):
    if func == "cos":
        bits = 64 * n - 6 - 2 * r * (k + 1)
    else:
        bits = 64 * n - 6 - r * (2 * k + 3)
    if bits <= 0:
        return 1
    return max(1, -(-bits // 64))


def cbits(c, L):
    """Fraction c as an L-limb fixed value, low word first"""
    v = (c.numerator << (64 * L)) // c.denominator
    return [(v >> (64 * i)) & ((1 << 64) - 1) for i in range(L)]


def emit_sqrhigh(o, dst, src, L):
    """dst = high L limbs of src^2"""
    if L == 1:
        o("    %s[0] = n_mulhi(%s[0], %s[0]);" % (dst, src, src))
    elif L == 2:
        o("    _fixed_sqrhi_2x2(&%s[1], &%s[0], %s[1], %s[0]);"
          % (dst, dst, src, src))
    else:
        o("    flint_mpn_sqrhigh(%s, %s, %d);" % (dst, src, L))


def emit_mulhigh(o, dst, a, b, L):
    """dst = high L limbs of a*b (a, b are L-limb pointers)"""
    if L == 1:
        o("    %s[0] = n_mulhi(%s[0], %s[0]);" % (dst, a, b))
    elif L == 2:
        o("    _fixed_mulhi_2x2_sloppy(&%s[1], &%s[0], %s[1], %s[0], %s[1], %s[0]);"
          % (dst, dst, a, a, b, b))
    else:
        o("    flint_mpn_mulhigh_n(%s, %s, %s, %d);" % (dst, a, b, L))


def emit_addsub(o, sub, dst, u, v, L):
    """dst = u -+ v over L limbs.  The longlong chains stop at width 8;
    wider sums fall back to one mpn call (their functions are dozens of
    mulhigh rows long -- one dispatched add is noise there)."""
    if L == 1:
        o("    %s[0] = %s[0] %s %s[0];" % (dst, u, "-" if sub else "+", v))
    elif L <= 8:
        o("    NN_%s_%d(%s, %s, %s);"
          % ("SUB" if sub else "ADD", L, dst, u, v))
    elif dst == u:
        o("    mpn_%s_n(%s, %s, %s, %d);"
          % ("sub" if sub else "add", dst, u, v, L))
    else:
        o("    mpn_%s_n(%s, %s, %s, %d);"
          % ("sub" if sub else "add", dst, u, v, L))


def emit(o, func, n, r, name):
    N = minN(func, n, r)
    L = [min(n, width(func, n, r, k)) for k in range(N)]
    for k in range(N - 1, 0, -1):              # inner levels never wider
        L[k - 1] = max(L[k - 1], L[k])

    sub = func in ("atan", "sin", "cos")       # alternating?
    o("/* %s, n = %d, r = %d: %d levels, widths %s */"
      % (func, n, r, N, L))
    o("static void")
    o("%s(nn_ptr res, nn_srcptr x)" % name)
    o("{")
    o("    ulong z[%d], V[%d], T[%d], P[%d], Z[2];" % (n, n, n, n))
    o("    ulong v, w1, w0, h, l;")
    o("    slong i;")
    o("")
    o("    (void) v; (void) w1; (void) w0; (void) h; (void) l;")
    o("    (void) i; (void) Z; (void) V; (void) T; (void) P;")
    o("")
    emit_sqrhigh(o, "z", "x", n)
    o("")

    # innermost level
    k = N - 1
    cur = L[k]
    c = cbits(coeff(func, k), cur)
    if cur == 1:
        o("    v = UWORD(0x%016x);" % c[0])
    elif cur == 2:
        o("    w1 = UWORD(0x%016x); w0 = UWORD(0x%016x);" % (c[1], c[0]))
    else:
        for i in range(cur):
            o("    V[%d] = UWORD(0x%016x);" % (i, c[i]))

    for k in range(N - 2, -1, -1):
        prev, cur = cur, L[k]
        c = cbits(coeff(func, k), cur)
        op = "-" if sub else "+"
        o("")
        o("    /* V_%d = c %s z V_%d   (%d x %d limbs) */" % (k, op, k + 1, cur, prev))
        if cur == 1:
            o("    v = n_mulhi(z[%d], v);" % (n - 1))
            o("    v = UWORD(0x%016x) %s v;" % (c[0], op))
        elif cur == 2 and prev == 1:
            o("    _fixed_mulhi_2x1(&h, &l, z[%d], z[%d], v);" % (n - 1, n - 2))
            if sub:
                o("    sub_ddmmss(w1, w0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
                  % (c[1], c[0]))
            else:
                o("    add_ssaaaa(w1, w0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
                  % (c[1], c[0]))
        elif cur == 2:
            o("    _fixed_mulhi_2x2_sloppy(&h, &l, z[%d], z[%d], w1, w0);"
              % (n - 1, n - 2))
            if sub:
                o("    sub_ddmmss(w1, w0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
                  % (c[1], c[0]))
            else:
                o("    add_ssaaaa(w1, w0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
                  % (c[1], c[0]))
        else:
            # promote the previous level into a cur-limb array (top aligned)
            if prev == 1:
                for i in range(cur - 1):
                    o("    P[%d] = UWORD(0);" % i)
                o("    P[%d] = v;" % (cur - 1))
            elif prev == 2:
                for i in range(cur - 2):
                    o("    P[%d] = UWORD(0);" % i)
                o("    P[%d] = w0; P[%d] = w1;" % (cur - 2, cur - 1))
            else:
                for i in range(cur - prev):
                    o("    P[%d] = UWORD(0);" % i)
                for i in range(prev):
                    o("    P[%d] = V[%d];" % (cur - prev + i, i))
            if cur >= 3:
                emit_mulhigh(o, "T", "z + %d" % (n - cur), "P", cur)
            else:
                for i in range(cur):
                    o("    Z[%d] = z[%d];" % (i, n - cur + i))
                emit_mulhigh(o, "T", "Z", "P", cur)
            for i in range(cur):
                o("    P[%d] = UWORD(0x%016x);" % (i, c[i]))
            emit_addsub(o, sub, "V", "P", "T", cur)

    # tail
    o("")
    if func == "cos":
        o("    /* cos = 1 - z V_0 */")
    else:
        o("    /* %s = x %s x^3 V_0 */" % (func, "-" if sub else "+"))
    if cur == 1:
        for i in range(n - 1):
            o("    P[%d] = UWORD(0);" % i)
        o("    P[%d] = v;" % (n - 1))
    elif cur == 2:
        for i in range(n - 2):
            o("    P[%d] = UWORD(0);" % i)
        o("    P[%d] = w0; P[%d] = w1;" % (n - 2, n - 1))
    else:
        for i in range(n - cur):
            o("    P[%d] = UWORD(0);" % i)
        for i in range(cur):
            o("    P[%d] = V[%d];" % (n - cur + i, i))

    if func == "cos":
        emit_mulhigh(o, "T", "z", "P", n)
        o("    /* 1 - T: below one the unit limb is zero */")
        o("    if (mpn_zero_p(T, %d))" % n)
        o("    {")
        o("        flint_mpn_zero(res, %d);" % n)
        o("        res[%d] = 1;" % n)
        o("    }")
        o("    else")
        o("    {")
        o("        mpn_neg(res, T, %d);" % n)
        o("        res[%d] = 0;" % n)
        o("    }")
    else:
        emit_mulhigh(o, "T", "z", "x", n)      # x^3
        emit_mulhigh(o, "V", "T", "P", n)
        emit_addsub(o, sub, "res", "x", "V", n)
    o("}")
    o("")
    return N


def emit_chain(o, func, n, r, zname):
    """emit the Horner chain for one family on an already-computed z,
    leaving the result in the register/array named by the return value"""
    N = minN(func, n, r)
    L = [min(n, width(func, n, r, k)) for k in range(N)]
    for k in range(N - 1, 0, -1):
        L[k - 1] = max(L[k - 1], L[k])
    sub = True                              # sin and cos both alternate
    pre = "s" if func == "sin" else "c"
    k = N - 1
    cur = L[k]
    c = cbits(coeff(func, k), cur)
    if cur == 1:
        o("    %sv = UWORD(0x%016x);" % (pre, c[0]))
    elif cur == 2:
        o("    %sw1 = UWORD(0x%016x); %sw0 = UWORD(0x%016x);"
          % (pre, c[1], pre, c[0]))
    else:
        for i in range(cur):
            o("    %sV[%d] = UWORD(0x%016x);" % (pre, i, c[i]))
    for k in range(N - 2, -1, -1):
        prev, cur = cur, L[k]
        c = cbits(coeff(func, k), cur)
        o("")
        o("    /* %s V_%d  (%d x %d) */" % (func, k, cur, prev))
        if cur == 1:
            o("    %sv = n_mulhi(%s[%d], %sv);" % (pre, zname, n - 1, pre))
            o("    %sv = UWORD(0x%016x) - %sv;" % (pre, c[0], pre))
        elif cur == 2 and prev == 1:
            o("    _fixed_mulhi_2x1(&h, &l, %s[%d], %s[%d], %sv);"
              % (zname, n - 1, zname, n - 2, pre))
            o("    sub_ddmmss(%sw1, %sw0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
              % (pre, pre, c[1], c[0]))
        elif cur == 2:
            o("    _fixed_mulhi_2x2_sloppy(&h, &l, %s[%d], %s[%d], %sw1, %sw0);"
              % (zname, n - 1, zname, n - 2, pre, pre))
            o("    sub_ddmmss(%sw1, %sw0, UWORD(0x%016x), UWORD(0x%016x), h, l);"
              % (pre, pre, c[1], c[0]))
        else:
            if prev == 1:
                for i in range(cur - 1):
                    o("    P[%d] = UWORD(0);" % i)
                o("    P[%d] = %sv;" % (cur - 1, pre))
            elif prev == 2:
                for i in range(cur - 2):
                    o("    P[%d] = UWORD(0);" % i)
                o("    P[%d] = %sw0; P[%d] = %sw1;" % (cur - 2, pre, cur - 1, pre))
            else:
                for i in range(cur - prev):
                    o("    P[%d] = UWORD(0);" % i)
                for i in range(prev):
                    o("    P[%d] = %sV[%d];" % (cur - prev + i, pre, i))
            if cur >= 3:
                emit_mulhigh(o, "T", "%s + %d" % (zname, n - cur), "P", cur)
            else:
                for i in range(cur):
                    o("    Z[%d] = %s[%d];" % (i, zname, n - cur + i))
                emit_mulhigh(o, "T", "Z", "P", cur)
            for i in range(cur):
                o("    P[%d] = UWORD(0x%016x);" % (i, c[i]))
            emit_addsub(o, True, "%sV" % pre, "P", "T", cur)
    return cur


def emit_sincos(o, n, r, name):
    """sin and cos together: they are both series in z = x^2, so ONE
    squaring serves both"""
    o("/* sin_cos, n = %d, r = %d */" % (n, r))
    o("static void")
    o("%s(nn_ptr ysin, nn_ptr ycos, nn_srcptr x)" % name)
    o("{")
    o("    ulong z[%d], sV[%d], cV[%d], T[%d], P[%d], Z[2];" % (n, n, n, n, n))
    o("    ulong sv, sw1, sw0, cv, cw1, cw0, h, l;")
    o("")
    o("    (void) sv; (void) sw1; (void) sw0; (void) cv; (void) cw1;")
    o("    (void) cw0; (void) h; (void) l; (void) Z; (void) T; (void) P;")
    o("    (void) sV; (void) cV;")
    o("")
    o("    /* the one squaring shared by both series */")
    emit_sqrhigh(o, "z", "x", n)
    o("")
    cs = emit_chain(o, "sin", n, r, "z")
    o("")
    cc = emit_chain(o, "cos", n, r, "z")
    o("")
    o("    /* sin = x - x^3 V_0 */")
    if cs == 1:
        for i in range(n - 1):
            o("    P[%d] = UWORD(0);" % i)
        o("    P[%d] = sv;" % (n - 1))
    elif cs == 2:
        for i in range(n - 2):
            o("    P[%d] = UWORD(0);" % i)
        o("    P[%d] = sw0; P[%d] = sw1;" % (n - 2, n - 1))
    else:
        for i in range(n - cs):
            o("    P[%d] = UWORD(0);" % i)
        for i in range(cs):
            o("    P[%d] = sV[%d];" % (n - cs + i, i))
    emit_mulhigh(o, "T", "z", "x", n)
    emit_mulhigh(o, "sV", "T", "P", n)
    emit_addsub(o, True, "ysin", "x", "sV", n)
    o("    ysin[%d] = 0;" % n)
    o("")
    o("    /* cos = 1 - z V_0 */")
    if cc == 1:
        for i in range(n - 1):
            o("    P[%d] = UWORD(0);" % i)
        o("    P[%d] = cv;" % (n - 1))
    elif cc == 2:
        for i in range(n - 2):
            o("    P[%d] = UWORD(0);" % i)
        o("    P[%d] = cw0; P[%d] = cw1;" % (n - 2, n - 1))
    else:
        for i in range(n - cc):
            o("    P[%d] = UWORD(0);" % i)
        for i in range(cc):
            o("    P[%d] = cV[%d];" % (n - cc + i, i))
    emit_mulhigh(o, "T", "z", "P", n)
    o("    if (mpn_zero_p(T, %d))" % n)
    o("    {")
    o("        flint_mpn_zero(ycos, %d);" % n)
    o("        ycos[%d] = 1;" % n)
    o("    }")
    o("    else")
    o("    {")
    o("        mpn_neg(ycos, T, %d);" % n)
    o("        ycos[%d] = 0;" % n)
    o("    }")
    o("}")
    o("")

def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    out = []
    o = out.append
    ent = []
    for func in FUNCS:
        for n in range(1, 5):
            for r in range(4, 49):
                if r > 64 * n - 16:
                    continue
                nm = "hs_%s_%d_%d" % (func, n, r)
                N = emit(o, func, n, r, nm)
                ent.append((nm, func, n, r, N))
    o("typedef void (*hs_fn)(nn_ptr, nn_srcptr);")
    o("typedef struct { hs_fn f; const char * name; int n; int r; int N; } hs_entry;")
    o("static hs_entry hs_all[] = {")
    for nm, func, n, r, N in ent:
        o('    { %s, "%s", %d, %d, %d },' % (nm, func, n, r, N))
    o("};")
    o("static int hs_count = %d;" % len(ent))
    print("\n".join(out))


if __name__ == "__main__":
    main()
