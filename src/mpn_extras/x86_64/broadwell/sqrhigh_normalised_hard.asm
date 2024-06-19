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

dnl TODO:
dnl  * Use SF flag from adc instead of using test

include(`config.m4')

define(`rp', `%rdi')
define(`ap', `%rsi')

define(`s0', `%rax')
define(`s1', `%rcx')
define(`s2', `%r8')
define(`s3', `%r9')
define(`s4', `%r10')
define(`s5', `%r11')
define(`s6', `%rbx')
define(`s7', `%rbp')
define(`s8', `%r12')
define(`s9', `%r13')
define(`s10', `%r14')
define(`s11', `%r15')

	TEXT
	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_normalised_1)
	mov	0*8(ap), %rdx
	mulx	%rdx, s0, s1
	xor	%edx, %edx
	test	s1, s1
	js	L(1)
	add	s0, s0
	inc	%edx
	adc	s1, s1
L(1):	mov	s1, 0*8(rp)
	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_normalised_2)
	C   0 1
	C 0 e
	C 1 x d
	mov	0*8(ap), %rdx
	mulx	1*8(ap), s2, s3		C a1 a0
	xor	R32(s1), R32(s1)
	mulx	%rdx, s0, s0		C a0^2
	C 2, 3, 1
	C 0

	mov	1*8(ap), %rdx
	add	s2, s2
	adc	s3, s3
	adc	R32(s1), R32(s1)
	mulx	%rdx, s4, s5		C a1^2
	C 2, 3, 1
	C 0, 4, 5

	xor	%edx, %edx
	add	s2, s0
	adc	s4, s3
	adc	s5, s1
	C 0, 3, 1

	js	L(2)
	add	s0, s0
	adc	s3, s3
	inc	%edx
	adc	s1, s1
L(2):	mov	s3, 0*8(rp)
	mov	s1, 1*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_normalised_3)
	C   0 1 2
	C 0
	C 1 h d
	C 2 x x d
	mov	1*8(ap), %rdx
	mulx	0*8(ap), s1, s1		C a1 a0
	mulx	%rdx, s0, s2		C a1^2
	C 1
	C 0, 2

	mov	2*8(ap), %rdx
	mulx	0*8(ap), s3, s4		C a2 a0
	mulx	1*8(ap), s5, ap		C a2 a1
	C (1, 3), (4, 5), ap
	C 0, 2

	add	s1, s3
	adc	s5, s4
	mulx	%rdx, s1, s5		C a2^2
	adc	$0, ap
	xor	%edx, %edx
	C 3, 4, ap, rdx
	C 0, 2, 1, 5

	add	s3, s3
	adc	s4, s4
	adc	ap, ap
	adc	%edx, %edx
	C 3, 4, ap, rdx
	C 0, 2, 1, 5

	add	s3, s0
	adc	s4, s2
	adc	ap, s1
	adc	%rdx, s5
	mov	$0, %edx
	C 0, 2, 1, 5

	js	L(3)
	add	s0, s0
	adc	s2, s2
	inc	%edx
	adc	s1, s1
	adc	s5, s5
L(3):	mov	s2, 0*8(rp)
	mov	s1, 1*8(rp)
	mov	s5, 2*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_normalised_4)
	C   0 1 2 3
	C 0
	C 1   e
	C 2 h x d
	C 3 x x x d
	push	s6
	push	s7

	mov	2*8(ap), %rdx
	mulx	0*8(ap), s0, s0		C a2 a0
	mulx	1*8(ap), s1, s2		C a2 a1
	C (0, 1), 2

	mov	3*8(ap), %rdx
	add	s1, s0
	mulx	1*8(ap), s5, s6		C a3 a1
	mulx	0*8(ap), s3, s4		C a3 a0
	adc	s5, s2
	mulx	2*8(ap), s5, s1		C a3 a2
	adc	$0, s6
	C (0, 3), (2, 4), (6, 5), 1

	add	s3, s0
	adc	s4, s2
	adc	s6, s5
	adc	$0, s1
	C 0, 2, 6, 1

	mov	1*8(ap), %rdx
	mulx	%rdx, s3, s3		C a1^2
	xor	R32(s6), R32(s6)
	C 0, 2, 5, 1, 6
	C 3

	add	s0, s0
	adc	s2, s2
	mov	2*8(ap), %rdx
	mulx	%rdx, s7, s4		C a2^2
	adc	s5, s5
	adc	s1, s1
	mov	3*8(ap), %rdx
	mulx	%rdx, %rdx, ap		C a3^2
	adc	R32(s6), R32(s6)
	C 0, 2, 5,   1,  6
	C 3, 7, 4, rdx, ap

	add	s3, s0
	adc	s7, s2
	adc	s4, s5
	adc	%rdx, s1
	adc	s6, ap
	C 0, 2, 5, 1, ap

	mov	$0, %edx
	js	L(4)
	add	s0, s0
	adc	s2, s2
	adc	s5, s5
	inc	%edx
	adc	s1, s1
	adc	ap, ap
L(4):	mov	s2, 0*8(rp)
	mov	s5, 1*8(rp)
	pop	s7
	mov	s1, 2*8(rp)
	mov	ap, 3*8(rp)
	pop	s6

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_normalised_5)
	C   0 1 2 3 4
	C 0
	C 1
	C 2   h d
	C 3 h x x d
	C 4 x x x x d

	mov	3*8(ap), %rdx
	mulx	0*8(ap), s0, s0		C a3 a0
	push	s6
	push	s7
	mulx	1*8(ap), s1, s2		C a3 a1
	push	s8
	mulx	2*8(ap), s3, s4		C a3 a2
	add	s1, s0
	adc	s3, s2
	C 0, 2, 4 + c

	mov	2*8(ap), %rdx
	mulx	1*8(ap), s5, s5		C a2 a1
	adc	$0, s4
	xor	R32(s1), R32(s1)
	C (0, 5), 2, 4

	mov	4*8(ap), %rdx
	mulx	0*8(ap), s6, s7		C a4 a0
	adcx	s5, s0
	mulx	1*8(ap), s8, s3		C a4 a1
	adox	s6, s0
	adcx	s7, s2
	mulx	2*8(ap), s5, s6		C a4 a2
	adox	s8, s2
	adcx	s3, s4
	mulx	3*8(ap), s7, s8		C a4 a3
	C 0, 2, (4 + o, 5), (6 + c, 7), 8

	adox	s5, s4
	mov	2*8(ap), %rdx
	adcx	s1, s6
	adox	s7, s6
	adcx	s0, s0
	mulx	%rdx, s3, s5		C a2^2
	adox	s1, s8
	adcx	s2, s2
	adcx	s4, s4
	adcx	s6, s6
	mov	3*8(ap), %rdx
	adcx	s8, s8
	adc	R32(s1), R32(s1)
	C 0, 2, 4, 6, 8, 1
	C 3, 5

	add	s3, s0
	mulx	%rdx, s7, s3		C a3^2
	adc	s5, s2
	adc	s7, s4
	mov	4*8(ap), %rdx
	mulx	%rdx, s5, s7		C a4^2
	adc	s6, s3
	adc	s8, s5
	mov	$0, %edx
	pop	s8
	adc	s7, s1
	pop	s7
	C 0, 2, 4, 3, 5, 1

	js	L(5)
	add	s0, s0
	adc	s2, s2
	adc	s4, s4
	inc	%edx
	adc	s3, s3
	adc	s5, s5
	adc	s1, s1
L(5):	mov	s2, 0*8(rp)
	mov	s4, 1*8(rp)
	pop	s6
	mov	s3, 2*8(rp)
	mov	s5, 3*8(rp)
	mov	s1, 4*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_normalised_6)
	C   0 1 2 3 4 5
	C 0
	C 1
	C 2     e
	C 3   h x d
	C 4 h x x x d
	C 5 x x x x x d

	mov	3*8(ap), %rdx
	mulx	1*8(ap), s0, s0		C a3 a1
	push	s6
	xor	R32(s1), R32(s1)
	mulx	2*8(ap), s2, s3		C a3 a2
	push	s7
	adcx	s2, s0
	C 0, 3 + c

	mov	4*8(ap), %rdx
	mulx	1*8(ap), s4, s5		C a4 a1
	push	s8
	adox	s4, s0
	adcx	s5, s3
	mulx	2*8(ap), s6, s7		C a4 a2
	push	s9
	adox	s6, s3
	adcx	s1, s7
	mulx	3*8(ap), s8, s2		C a4 a3
	adox	s8, s7
	mulx	0*8(ap), s9, s9		C a4 a0
	adox	s1, s2
	C (0, 9), 3, 7, 2

	mov	5*8(ap), %rdx
	mulx	0*8(ap), s4, s5		C a5 a0
	adcx	s9, s0
	mulx	1*8(ap), s6, s8		C a5 a1
	adox	s4, s0
	adcx	s5, s3
	mulx	2*8(ap), s9, s4		C a5 a2
	adox	s6, s3
	adcx	s8, s7
	mulx	3*8(ap), s5, s6		C a5 a3
	adox	s9, s7
	adcx	s4, s2
	mulx	4*8(ap), s8, s9		C a5 a4
	adox	s5, s2
	adcx	s1, s6
	C 0, 3, 7, 2, (6 + o, 8), 9

	adox	s8, s6
	adcx	s0, s0
	mov	2*8(ap), %rdx
	adox	s1, s9
	adc	s3, s3
	mulx	%rdx, s4, s4		C a2^2
	adc	s7, s7
	mov	3*8(ap), %rdx
	adc	s2, s2
	mulx	%rdx, s5, s8		C a3^2
	adc	s6, s6
	adc	s9, s9
	mov	4*8(ap), %rdx
	adc	R32(s1), R32(s1)
	C 0, 3, 7, 2, 6, 9, 1
	C 4, 5, 8

	add	s4, s0
	adc	s5, s3
	mulx	%rdx, s4, s5		C a4^2
	adc	s8, s7
	mov	5*8(ap), %rdx
	mulx	%rdx, s8, ap		C a5^2
	adc	s4, s2
	adc	s6, s5
	mov	$0, %edx
	adc	s9, s8
	adc	ap, s1
	pop	s9
	C 0, 3, 7, 2, 5, 8, 1

	js	L(6)
	add	s0, s0
	adc	s3, s3
	adc	s7, s7
	inc	%edx
	adc	s2, s2
	adc	s5, s5
	adc	s8, s8
	adc	s1, s1
L(6):	mov	s3, 0*8(rp)
	mov	s7, 1*8(rp)
	mov	s2, 2*8(rp)
	mov	s5, 3*8(rp)
	mov	s8, 4*8(rp)
	pop	s8
	pop	s7
	mov	s1, 5*8(rp)
	pop	s6

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_normalised_7)
	C   0 1 2 3 4 5 6
	C 0
	C 1
	C 2
	C 3     h d
	C 4   h x x d
	C 5 h x x x x d
	C 6 x x x x x x d
	push	s6

	mov	3*8(ap), %rdx
	mulx	2*8(ap), s0, s0		C a3 a2
	push	s7
	C 0

	xor	R32(s1), R32(s1)

	mov	4*8(ap), %rdx
	mulx	2*8(ap), s2, s3		C a4 a2
	push	s8
	adcx	s2, s0
	mulx	3*8(ap), s4, s5		C a4 a3
	push	s9
	adcx	s4, s3
	mulx	1*8(ap), s6, s6		C a4 a1
	push	s10
	adcx	s1, s5
	C (0, 6), 3, 5

	mov	5*8(ap), %rdx
	mulx	0*8(ap), s7, s7		C a5 a0
	adcx	s6, s0
	mulx	1*8(ap), s8, s9		C a5 a1
	adox	s7, s0
	adcx	s9, s3
	mulx	2*8(ap), s10, s2	C a5 a2
	adox	s10, s3
	adcx	s2, s5
	mulx	3*8(ap), s4, s6		C a5 a3
	adox	s4, s5
	adcx	s1, s6
	mulx	4*8(ap), s7, s9		C a5 a4
	adox	s7, s6
	C (0, 8), 3, 5, 6, 9 + c

	mov	6*8(ap), %rdx
	mulx	0*8(ap), s10, s2	C a6 a0
	adox	s1, s9
	adcx	s8, s0
	mulx	1*8(ap), s4, s7		C a6 a1
	adox	s10, s0
	adcx	s2, s3
	mulx	2*8(ap), s8, s10	C a6 a2
	adox	s4, s3
	adcx	s7, s5
	mulx	3*8(ap), s2, s4		C a6 a3
	adox	s8, s5
	adcx	s10, s6
	mulx	4*8(ap), s7, s8		C a6 a4
	adox	s2, s6
	adcx	s4, s9
	mulx	5*8(ap), s10, s2	C a6 a5
	adox	s7, s9
	adcx	s1, s8
	C 0, 3, 5, 6, 9, (8 + o, 10), 2

	mov	3*8(ap), %rdx
	adox	s10, s8
	adcx	s0, s0
	adox	s1, s2
	mulx	%rdx, s4, s7		C a3^2
	adc	s3, s3
	adc	s5, s5
	mov	4*8(ap), %rdx
	adc	s6, s6
	adc	s9, s9
	adc	s8, s8
	adc	s2, s2
	adc	R32(s1), R32(s1)
	C 0, 3, 5, 6, 9, 8, 2, 1
	C 4, 7

	add	s4, s0
	mulx	%rdx, s10, s4		C a4^2
	mov	5*8(ap), %rdx
	adc	s7, s3
	adc	s10, s5
	mulx	%rdx, s7, s10		C a5^2
	mov	6*8(ap), %rdx
	adc	s4, s6
	adc	s9, s7
	mulx	%rdx, s4, ap		C a6^2
	adc	s10, s8
	adc	s4, s2
	pop	s10
	mov	$0, %edx
	adc	ap, s1
	pop	s9
	C 0, 3, 5, 6, 7, 8, 2, 1

	js	L(7)
	add	s0, s0
	adc	s3, s3
	adc	s5, s5
	adc	s6, s6
	inc	%edx
	adc	s7, s7
	adc	s8, s8
	adc	s2, s2
	adc	s1, s1
L(7):	mov	s3, 0*8(rp)
	mov	s5, 1*8(rp)
	mov	s6, 2*8(rp)
	mov	s7, 3*8(rp)
	mov	s8, 4*8(rp)
	pop	s8
	mov	s2, 5*8(rp)
	pop	s7
	mov	s1, 6*8(rp)
	pop	s6

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_normalised_8)
	C   0 1 2 3 4 5 6 7
	C 0
	C 1
	C 2
	C 3       e
	C 4     h x d
	C 5   h x x x d
	C 6 h x x x x x d
	C 7 x x x x x x x d

	xor	R32(s1), R32(s1)
	mov	4*8(ap), %rdx
	mulx	2*8(ap), s0, s0		C a4 a2
	push	s6
	mulx	3*8(ap), s2, s3		C a4 a3
	push	s7
	adcx	s2, s0
	C 0, 3 + c

	mov	5*8(ap), %rdx
	mulx	2*8(ap), s4, s5		C a5 a2
	push	s8
	adcx	s5, s3
	mulx	3*8(ap), s6, s7		C a5 a3
	push	s9
	adcx	s1, s7
	adox	s4, s0
	mulx	4*8(ap), s8, s9		C a5 a4
	push	s10
	adox	s6, s3
	mulx	1*8(ap), s10, s10	C a5 a1
	push	s11
	adox	s8, s7
	C (0, 10), 3, 7, 9 + o

	mov	6*8(ap), %rdx
	mulx	1*8(ap), s11, s2	C a6 a1
	adox	s1, s9
	adcx	s10, s0
	mulx	2*8(ap), s4, s5		C a6 a2
	adox	s11, s0
	adcx	s2, s3
	mulx	3*8(ap), s6, s8		C a6 a3
	adox	s4, s3
	adcx	s5, s7
	mulx	4*8(ap), s10, s11	C a6 a4
	adox	s6, s7
	adcx	s8, s9
	mulx	5*8(ap), s2, s4		C a6 a5
	adox	s10, s9
	adcx	s1, s11
	mulx	0*8(ap), s5, s5		C a6 a0
	adox	s2, s11
	C (0, 5), 3, 7, 9, 11, 4 + o

	mov	7*8(ap), %rdx
	mulx	0*8(ap), s2, s6		C a7 a0
	adox	s1, s4
	adcx	s5, s0
	mulx	1*8(ap), s8, s10	C a7 a1
	adox	s2, s0
	adcx	s6, s3
	mulx	2*8(ap), s5, s2		C a7 a2
	adox	s8, s3
	adcx	s10, s7
	mulx	3*8(ap), s6, s8		C a7 a3
	adox	s5, s7
	adcx	s2, s9
	mulx	4*8(ap), s10, s5	C a7 a4
	adox	s6, s9
	adcx	s8, s11
	mulx	5*8(ap), s2, s6		C a7 a5
	adox	s10, s11
	adcx	s5, s4
	mulx	6*8(ap), s8, s10	C a7 a6
	adox	s2, s4
	adcx	s1, s6
	C 0, 3, 7, 9, 11, 4, (6 + o, 8), 10

	adox	s8, s6
	adcx	s0, s0
	mov	3*8(ap), %rdx
	adox	s1, s10
	adc	s3, s3
	mulx	%rdx, s2, s2		C a3^2
	adc	s7, s7
	mov	4*8(ap), %rdx
	adc	s9, s9
	mulx	%rdx, s5, s8		C a4^2
	adc	s11, s11
	adc	s4, s4
	adc	s6, s6
	mov	5*8(ap), %rdx
	adc	s10, s10
	adc	R32(s1), R32(s1)
	C 0, 3, 7, 9, 11, 4, 6, 10, 1
	C 2, 5, 8

	add	s2, s0
	adc	s5, s3
	mulx	%rdx, s2, s5		C a5^2
	adc	s8, s7
	mov	6*8(ap), %rdx
	adc	s2, s9
	mulx	%rdx, s8, s2		C a6^2
	adc	s5, s11
	mov	7*8(ap), %rdx
	mulx	%rdx, ap, s5		C a7^2
	adc	s8, s4
	adc	s6, s2
	mov	$0, %edx
	adc	s10, ap
	adc	s5, s1
	C 0, 3, 7, 9, 11, 4, 2, ap, 1

	js	L(8)
	add	s0, s0
	adc	s3, s3
	adc	s7, s7
	adc	s9, s9
	inc	%edx
	adc	s11, s11
	adc	s4, s4
	adc	s2, s2
	adc	ap, ap
	adc	s1, s1
L(8):	mov	s3, 0*8(rp)
	mov	s7, 1*8(rp)
	mov	s9, 2*8(rp)
	mov	s11, 3*8(rp)
	pop	s11
	pop	s10
	mov	s4, 4*8(rp)
	pop	s9
	mov	s2, 5*8(rp)
	pop	s8
	pop	s7
	mov	ap, 6*8(rp)
	pop	s6
	mov	s1, 7*8(rp)

	ret
EPILOGUE()
