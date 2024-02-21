dnl  X64-64 mpn_mullo_basecase optimised for Intel Broadwell.

dnl  Contributed to the GNU project by Torbjorn Granlund.

dnl  Copyright 2017 Free Software Foundation, Inc.

dnl  This file is part of the GNU MP Library.
dnl
dnl  The GNU MP Library is free software; you can redistribute it and/or modify
dnl  it under the terms of either:
dnl
dnl    * the GNU Lesser General Public License as published by the Free
dnl      Software Foundation; either version 3 of the License, or (at your
dnl      option) any later version.
dnl
dnl  or
dnl
dnl    * the GNU General Public License as published by the Free Software
dnl      Foundation; either version 2 of the License, or (at your option) any
dnl      later version.
dnl
dnl  or both in parallel, as here.
dnl
dnl  The GNU MP Library is distributed in the hope that it will be useful, but
dnl  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
dnl  for more details.
dnl
dnl  You should have received copies of the GNU General Public License and the
dnl  GNU Lesser General Public License along with the GNU MP Library.  If not,
dnl  see https://www.gnu.org/licenses/.
dnl
dnl   Copyright (C) 2024 Albin Ahlbäck
dnl
dnl   This file is part of FLINT.
dnl
dnl   FLINT is free software: you can redistribute it and/or modify it under
dnl   the terms of the GNU Lesser General Public License (LGPL) as published
dnl   by the Free Software Foundation; either version 3 of the License, or
dnl   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
dnl

include(`config.m4')

define(`rp',	   `%rdi')
define(`ap',	   `%rsi')
define(`bp_param', `%rdx')
define(`n',	   `%rcx')

define(`bp',	   `%r8')
define(`jmpreg',   `%r9')
define(`nn',	   `%r10')
define(`m',	   `%r13')
define(`mm',	   `%r14')

define(`rx',	   `%rax')

define(`r0',	   `%r11')
define(`r1',	   `%rbx')
define(`r2',	   `%rbp')
define(`r3',	   `%r12')

dnl Idea: Do similar to mpn_mullo_basecase for Skylake.

	TEXT
	ALIGN(32)
PROLOGUE(_flint_mpn_mulhigh_basecase)
	mov	bp_param, bp
	lea	-1*8(ap,n,8), ap	C ap += n - 1

	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	C Initial triangle
	C       h
	C     h x
	C   h x x
	C h x x x
	C x x x x
define(`s0', `jmpreg')
define(`s1', `m')
define(`s2', `mm')
define(`s3', `nn')
	mov	0*8(bp), %rdx
	xor	R32(s3), R32(s3)
	mulx	-1*8(ap), rx, rx
	mulx	0*8(ap), s0, r0
	add	s0, rx
	adc	s3, r0

	mov	1*8(bp), %rdx
	mulx	-2*8(ap), s1, s1
	mulx	-1*8(ap), r3, r2
	mulx	0*8(ap), s0, r1
	add	r3, rx
	adc	s0, r0
	adc	s3, r1
	add	s1, rx
	adc	r2, r0
	adc	s3, r1

	mov	2*8(bp), %rdx
	mulx	-3*8(ap), s0, s0
	mulx	-2*8(ap), r3, s1
	add	s0, rx
	adc	s1, r0
	mulx	-1*8(ap), s0, s1
	mulx	0*8(ap), %rdx, r2
	adc	s1, r1
	adc	s3, r2
	add	r3, rx
	adc	s0, r0
	adc	%rdx, r1
	adc	s3, r2

	mov	3*8(bp), %rdx
	mulx	-4*8(ap), s1, s1
	mulx	-3*8(ap), s0, s2
	add	s1, rx
	adc	s2, r0
	mulx	-2*8(ap), s1, r3
	mulx	-1*8(ap), s2, s3
	adc	r3, r1
	adc	s3, r2
	mulx	0*8(ap), %rdx, r3
	adc	$0, r3
	add	s0, rx
	adc	s1, r0
	adc	s2, r1
	mov	r0, 0*8(rp)
	mov	r1, 1*8(rp)
	adc	%rdx, r2
	adc	$0, r3
	mov	r2, 2*8(rp)
	mov	r3, 3*8(rp)
undefine(`s0')
undefine(`s1')
undefine(`s2')
undefine(`s3')

	C Addmul chains
	C - m = -8 * n_cur	(n_cur is the 4 at the start)
	C - mm = -8 * (n - 1)	(where n is the original n)
	C - n keeps track of how many loops to do in the addmul-loop.
	C - nn keeps track of initial n between loops.
	lea	-1*8(,n,8), R32(mm)
	lea	4*8(bp), bp
	lea	-3*8(ap), ap
	mov	$-4*8, m		C m <- -8 * 4
	neg	mm			C mm <- -8 * (n - 1)
	mov	0*8(bp), %rdx
	xor	R32(nn), R32(nn)	C nn <- 0
	xor	R32(n), R32(n)		C n <- 0
	mulx	-2*8(ap), r1, r1
	adcx	r1, rx

L(f4):	mulx	-1*8(ap), r2, r3
	mulx	0*8(ap), r0, r1
	adox	r2, rx
	adcx	r3, r0
	lea	3*8(ap), ap
	lea	-5*8(rp), rp
	lea	L(f5)(%rip), jmpreg
	jmp	L(b4)

L(f0):	mulx	-1*8(ap), r2, r3
	mulx	0*8(ap), r0, r1
	adox	r2, rx
	adcx	r3, r0
	lea	-1*8(ap), ap
	lea	-1*8(rp), rp
	lea	L(f1)(%rip), jmpreg
	jmp	L(b0)

L(f1):	mulx	-1*8(ap), r0, r1
	mulx	0*8(ap), r2, r3
	adox	r0, rx
	adcx	r1, r2
	lea	1(nn), R32(nn)
	lea	1(n), R32(n)
	lea	L(f2)(%rip), jmpreg
	jmp	L(b1)

L(f7):	mulx	-1*8(ap), r0, r1
	mulx	0*8(ap), r2, r3
	adox	r0, rx
	adcx	r1, r2
	lea	-2*8(ap), ap
	lea	-2*8(rp), rp
	lea	L(f0)(%rip), jmpreg
	jmp	L(b7)

L(f2):	mulx	-1*8(ap), r2, r3
	mulx	0*8(ap), r0, r1
	adox	r2, rx
	adcx	r3, r0
	lea	1*8(ap), ap
	lea	1*8(rp), rp
	mulx	0*8(ap), r2, r3
	lea	L(f3)(%rip), jmpreg
	jmp	L(b2)

L(end):	adox	0*8(rp), r2
	mov	r2, 0*8(rp)
	adox	n, r3		C n = 0
	adc	n, r3		C n = 0
	add	m, ap		C Reset ap
	mov	r3, 1*8(rp)
	lea	-1*8(m), m
	lea	1*8(bp), bp	C Increase bp
	lea	2*8(rp,m), rp	C Reset rp
	mov	0*8(bp), %rdx	C Load bp
	cmp	R32(m), R32(mm)
	jge	L(jmp)
	C If |m| < |mm|: goto jmpreg, but first do high part
	or	R32(nn), R32(n)	C Reset n, CF and OF
	mulx	-2*8(ap), r1, r1
	adcx	r1, rx
	jmp	*jmpreg
	C If |m| > |mm|: goto fin
L(jmp):	jg	L(fin)
	C If |m| = |mm|: goto jmpreg
	or	R32(nn), R32(n)	C Reset n, clear CF and OF
	jmp	*jmpreg

	ALIGN(32)
L(b2):	adox	-1*8(rp), r0
	adcx	r1, r2
	mov	r0, -1*8(rp)
	jrcxz	L(end)	C Jump if n = 0
L(b1):	mulx	1*8(ap), r0, r1
	adox	0*8(rp), r2
	lea	-1(n), R32(n)
	mov	r2, 0*8(rp)
	adcx	r3, r0
L(b0):	mulx	2*8(ap), r2, r3
	adcx	r1, r2
	adox	1*8(rp), r0
	mov	r0, 1*8(rp)
L(b7):	mulx	3*8(ap), r0, r1
	lea	8*8(ap), ap
	adcx	r3, r0
	adox	2*8(rp), r2
	mov	r2, 2*8(rp)
L(b6):	mulx	-4*8(ap), r2, r3
	adox	3*8(rp), r0
	adcx	r1, r2
	mov	r0, 3*8(rp)
L(b5):	mulx	-3*8(ap), r0, r1
	adcx	r3, r0
	adox	4*8(rp), r2
	mov	r2, 4*8(rp)
L(b4):	mulx	-2*8(ap), r2, r3
	adox	5*8(rp), r0
	adcx	r1, r2
	mov	r0, 5*8(rp)
L(b3):	adox	6*8(rp), r2
	mulx	-1*8(ap), r0, r1
	mov	r2, 6*8(rp)
	lea	8*8(rp), rp
	adcx	r3, r0
	mulx	0*8(ap), r2, r3
	jmp	L(b2)

L(f6):	mulx	-1*8(ap), r2, r3
	mulx	0*8(ap), r0, r1
	adox	r2, rx
	adcx	r3, r0
	lea	5*8(ap), ap
	lea	-3*8(rp), rp
	lea	L(f7)(%rip), jmpreg
	jmp	L(b6)

L(f5):	mulx	-1*8(ap), r0, r1
	mulx	0*8(ap), r2, r3
	adox	r0, rx
	adcx	r1, r2
	lea	4*8(ap), ap
	lea	-4*8(rp), rp
	lea	L(f6)(%rip), jmpreg
	jmp	L(b5)

L(f3):	mulx	-1*8(ap), r0, r1
	mulx	0*8(ap), r2, r3
	adox	r0, rx
	adcx	r1, r2
	lea	2*8(ap), ap
	lea	-6*8(rp), rp
	lea	L(f4)(%rip), jmpreg
	jmp	L(b3)

L(fin):	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()
