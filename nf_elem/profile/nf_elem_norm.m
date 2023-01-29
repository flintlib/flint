/*============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

===============================================================================*/

/*
Profiling MAGMA polynomial multiplication in Z[x].

Usage: run magma with the -b flag to prevent the start up banner, i.e.

   magma -b magma-profile.m > output.prof

(C) 2007 David Harvey + Bill Hart, GPL

*/

target_name := "NFElemNorm";
target_description := "MAGMA number field norm over various degree";

BITS := 10;        // bits per coefficient
ratio := 1.1;      // ratio between consecutive lengths
monic := false;    // monic integral case or generic rational case

R<x>:=PolynomialRing(Rationals());
S<y>:=PolynomialRing(Integers());

// Timing runs need to last at least this many microseconds to be counted:
DURATION_THRESHOLD := 200000;
// Microseconds per timing run that the prof2d_sample function aims for:
DURATION_TARGET := 300000;


forward prof2d_sample;

/*
This function should run count iterations at position (x, y),
and return the total time in seconds, using the Cputime() function.
*/
function sampler(length, bits, count)

    // first time random element generation + multiplication

    countmod := 100;
    if length ge 50 then
       countmod := 10;
    end if;
    if length ge 500 then
       countmod := 4;
    end if;
    time1 := Cputime();
    for i := 1 to count do
      if (i-1) mod countmod eq 0 then
         p:=Polynomial([(-1)^Random(1)*RandomBits(bits): x in [1..length-1]]);
         if monic then
            p1:=S!p + x^(length-1);
         else
            den:=RandomBits(bits);
            if den eq 0 then
               den:=1;
            end if;
            p1:=((R!p)+(-1)^Random(1)*RandomBits(bits)*x^(length-1))/den;
         end if;
         if Degree(p1) eq 0 then
            if monic then
               p1:=x+1;
            else
               p1:=y+1;
            end if;
         end if;
         K:=NumberField(p1:Check:=false);
         a:=Random(K,2^bits);
         if monic then
            a:=Numerator(a);
         end if;
      end if;
      c:=Norm(a);
    end for;

    time2 := Cputime();

    // now time just the random element generation

    for i := 1 to count do
      if (i-1) mod countmod eq 0 then
         p:=Polynomial([(-1)^Random(1)*RandomBits(bits): x in [1..length-1]]);
         if monic then
            p1:=S!p + x^(length-1);
         else
            den:=RandomBits(bits);
            if den eq 0 then
               den:=1;
            end if;
            p1:=((R!p)+(-1)^Random(1)*RandomBits(bits)*x^(length-1))/den;
         end if;
         if Degree(p1) eq 0 then
            if monic then
               p1:=x+1;
            else
               p1:=y+1;
            end if;
         end if;
         K:=NumberField(p1:Check:=false);
         a:=Random(K,2^bits);
         if monic then
            a:=Numerator(a);
         end if;
      end if;
      c:=a;
    end for;

    time3 := Cputime();
    return (time2 - time1) - (time3 - time2);
end function;


/*
This function should loop over appropriate combinations of (x, y),
and call prof2d_sample(x, y) for each one.
*/
procedure driver()
   length:=2;
   while length lt 10000 do
      monic:=0;
      prof2d_sample(length, BITS);
      monic:=1;
      prof2d_sample(length, BITS);
      length := Ceiling(length*ratio);
   end while;
end procedure;



/************************************************************************

 This last section is the generic profiling code. Just leave this
 stuff alone.

************************************************************************/

/*
Formats in scientific notation with 3 decimal places
*/
function format_sci(x)
    L := Floor(Log(10, x));
    x := x / 10^L;
    s := Sprintf("%.3oe", x);
    if L lt 0 then
      s := s cat "-";
    else
      s := s cat "+";
    end if;

    s := s cat Sprintf("%o", Floor(Abs(L / 10)));
    s := s cat Sprintf("%o", (Abs(L) mod 10));

    return s;
end function;


procedure prof2d_sample(x, y)
    // number of timings that were at least DURATION_THRESHOLD microseconds:
    good_count := 0;

    // first try just a few loops
    num_trials := 4;
    last_time := sampler(x, y, num_trials) * 1000000.0;

    max_time := 0;
    min_time := 0;

    // loop until we have enough good times
    while true do
      per_trial := last_time / num_trials;

      // if the last recorded time was long enough, record it
      if last_time gt DURATION_THRESHOLD then
          if good_count gt 0 then
            max_time := Max(max_time, per_trial);
            min_time := Min(min_time, per_trial);
          else
            max_time := per_trial;
            min_time := per_trial;
          end if;

          good_count := good_count + 1;
          if good_count eq 5 then
            // we've got enough data
            // print it out and return
            print Sprintf("%o\t%o\t%o\t%o", x, y, format_sci(min_time), format_sci(max_time));
            return;
          end if;
      end if;

      // adjust num_trials so that the elapsed time gravitates towards
      // DURATION_TARGET; num_trials can be changed by a factor of
      // at most 25%, and must be at least 1
      if last_time lt 0.0001 then
          last_time := 0.0001;
      end if;
      adjust_ratio := 1.0 * DURATION_TARGET / last_time;
      if adjust_ratio gt 1.25 then
          adjust_ratio := 1.25;
      end if;
      if adjust_ratio lt 0.75 then
          adjust_ratio := 0.75;
      end if;
      num_trials := Ceiling(adjust_ratio * num_trials);
      // just to be safe:
      if num_trials eq 0 then
          num_trials := 1;
      end if;

      // run another trial
      last_time := sampler(x, y, num_trials) * 1000000.0;
    end while;
end procedure;


procedure print_header()
    print "FLINT profile output";
    print "";
    print "TIMESTAMP: (todo: write code to generate timestamp)";
    print "MACHINE: (todo: write code to get machine from environment var)";

    print "";
    print "MODULE: magma";
    print "TARGET:", target_name;
    print "";
    print "DESCRIPTION:";
    print target_description;
    print "";
    print "============================================== begin data";
    
end procedure;


print_header();
driver();

quit;

// ------------- end of file ------------------------------------
