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

	TEXT
	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_1)
	mov	0*8(ap), %rdx
	mulx	%rdx, s0, s1
	mov	s1, 0*8(rp)
	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_2)
	mov	0*8(ap), %rdx
	mulx	1*8(ap), s3, s4
	mulx	%rdx, s0, s0
	mov	1*8(ap), %rdx
	add	s3, s0
	mulx	%rdx, s1, s2
	adc	s4, s1
	adc	$0, s2
	add	s3, s0
	adc	s4, s1
	adc	$0, s2

	mov	s1, 0*8(rp)
	mov	s2, 1*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_3)
	C   0 1 2
	C 0
	C 1 h d
	C 2 x x d
	mov	1*8(ap), %rdx
	mulx	0*8(ap), s1, s1		C a1 a0
	mulx	%rdx, s0, s2		C a1^2

	mov	2*8(ap), %rdx
	mulx	0*8(ap), s3, s4		C a2 a0
	mulx	1*8(ap), s5, ap		C a2 a1

	add	s1, s3
	adc	s5, s4
	adc	$0, ap
	xor	R32(s1), R32(s1)

	add	s3, s3
	adc	s4, s4
	adc	ap, ap
	adc	s1, s1

	add	s3, s0
	adc	s4, s2

	mulx	%rdx, s3, s4		C a2^2

	mov	s2, 0*8(rp)
	adc	ap, s3
	adc	s1, s4
	mov	s3, 1*8(rp)
	mov	s4, 2*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_4)
	C   0 1 2 3
	C 0
	C 1   e
	C 2 h x d
	C 3 x x x d
	mov	2*8(ap), %rdx
	mulx	0*8(ap), s0, s0		C a2 a0
	mulx	1*8(ap), s3, s2		C a2 a1
	add	s3, s0
	adc	$0, s2
	C 0, 2

	mov	3*8(ap), %rdx
	mulx	0*8(ap), s1, s4		C a3 a0
	mulx	1*8(ap), s3, s5		C a3 a1
	add	s1, s0
	adc	s4, s2
	adc	$0, s5
	xor	R32(s1), R32(s1)
	add	s3, s2
	mulx	2*8(ap), %rdx, s4	C a3 a2
	adc	%rdx, s5
	adc	s1, s4
	C 0, 2, 5, 4

	add	s0, s0
	adc	s2, s2
	adc	s5, s5
	adc	s4, s4
	adc	R32(s1), R32(s1)
	C 0, 2, 5, 4, 1

	mov	1*8(ap), %rdx
	mulx	%rdx, s3, s3		C a1^2
	add	s3, s0

	mov	2*8(ap), %rdx
	mulx	%rdx, %rdx, s3		C a2^2
	adc	%rdx, s2
	adc	s3, s5
	mov	3*8(ap), %rdx
	mov	s2, 0*8(rp)
	mov	s5, 1*8(rp)

	mulx	%rdx, %rdx, s3		C a3^2
	adc	%rdx, s4
	adc	s3, s1
	mov	s4, 2*8(rp)
	mov	s1, 3*8(rp)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_5)
	C   0 1 2 3 4
	C 0
	C 1
	C 2   h d
	C 3 h x x d
	C 4 x x x x d
	push	%rbx
	push	%rbp

	mov	2*8(ap), %rdx
	mulx	1*8(ap), s0, s0		C a2 a1
	C 0

	xor	R32(s1), R32(s1)

	mov	3*8(ap), %rdx
	mulx	0*8(ap), s5, s5		C a3 a0
	mulx	1*8(ap), s2, s3		C a3 a1
	mulx	2*8(ap), s6, s4		C a3 a2
	add	s5, s0
	adc	s6, s3
	adc	s1, s4
	add	s2, s0
	adc	s1, s3
	adc	s1, s4
	C 0, 3, 4

	mov	4*8(ap), %rdx
	mulx	0*8(ap), s7, s2		C a4 a0
	mulx	1*8(ap), s6, s5		C a4 a1
	add	s7, s0
	adc	s2, s3
	mulx	2*8(ap), s7, s2		C a4 a2
	adc	s5, s4
	adc	s1, s2
	add	s6, s3
	adc	s7, s4
	mulx	3*8(ap), %rdx, s5	C a4 a3
	adc	%rdx, s2
	adc	s1, s5
	C 0, 3, 4, 2, 5

	add	s0, s0
	adc	s3, s3
	adc	s4, s4
	adc	s2, s2
	adc	s5, s5
	adc	R32(s1), R32(s1)
	C 0, 3, 4, 2, 5, 1

	mov	2*8(ap), %rdx
	mulx	%rdx, s6, s7		C a2^2
	add	s6, s0
	adc	s7, s3

	mov	3*8(ap), %rdx
	mov	s3, 0*8(rp)
	mulx	%rdx, s6, s7		C a3^2
	adc	s6, s4
	adc	s7, s2

	mov	4*8(ap), %rdx
	mov	s4, 1*8(rp)
	mov	s2, 2*8(rp)
	mulx	%rdx, s6, s7		C a4^2
	adc	s6, s5
	adc	s7, s1
	mov	s5, 3*8(rp)
	mov	s1, 4*8(rp)

	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_6)
	C   0 1 2 3 4 5
	C 0
	C 1
	C 2     e
	C 3   h x d
	C 4 h x x x d
	C 5 x x x x x d
	push	%rbx
	push	%rbp
	push	%r12

	mov	3*8(ap), %rdx
	xor	R32(s1), R32(s1)
	mulx	1*8(ap), s0, s0		C a3 a1
	mulx	2*8(ap), s2, s3		C a3 a2
	add	s2, s0
	adc	s1, s3
	C 0, 3

	mov	4*8(ap), %rdx
	mulx	0*8(ap), s4, s4		C a4 a0
	mulx	1*8(ap), s2, s5		C a4 a1
	mulx	2*8(ap), s6, s7		C a4 a2
	mulx	3*8(ap), %rdx, s8	C a4 a3
	add	s4, s0
	adc	s5, s3
	adc	s1, s7
	add	s2, s0
	adc	s6, s3
	adc	%rdx, s7
	adc	s1, s8
	C 0, 3, 7, 8

	mov	5*8(ap), %rdx
	test	%al, %al
	mulx	0*8(ap), s2, s4		C a5 a0
	mulx	1*8(ap), s5, s6		C a5 a1
	adox	s2, s0
	adox	s4, s3
	adcx	s5, s3
	adcx	s6, s7
	mulx	2*8(ap), s2, s4		C a5 a2
	mulx	3*8(ap), s5, s6		C a5 a3
	adox	s2, s7
	adox	s4, s8
	adox	s1, s6
	adc	s5, s8
	mulx	4*8(ap), %rdx, s4	C a5 a4
	adc	%rdx, s6
	adc	s1, s4
	C 0, 3, 7, 8, 6, 4

	add	s0, s0
	mov	2*8(ap), %rdx
	adc	s3, s3
	adc	s7, s7
	adc	s8, s8
	adc	s6, s6
	adc	s4, s4
	adc	R32(s1), R32(s1)
	C 0, 3, 7, 8, 6, 4, 1
	C (2, 5)

	mulx	%rdx, s2, s2		C a2^2
	mov	3*8(ap), %rdx
	add	s2, s0
	mulx	%rdx, s2, s5		C a3^2
	adc	s2, s3
	adc	s5, s7
	mov	4*8(ap), %rdx
	mov	s3, 0*8(rp)
	mulx	%rdx, s2, s5		C a4^2
	mov	s7, 1*8(rp)
	mov	5*8(ap), %rdx
	adc	s2, s8
	adc	s5, s6
	mov	s8, 2*8(rp)
	mulx	%rdx, s2, s5		C a5^2
	mov	s6, 3*8(rp)
	adc	s2, s4
	adc	s5, s1
	mov	s4, 4*8(rp)
	mov	s1, 5*8(rp)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_7)
	C   0 1 2 3 4 5 6
	C 0
	C 1
	C 2
	C 3     h d
	C 4   h x x d
	C 5 h x x x x d
	C 6 x x x x x x d
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	mov	3*8(ap), %rdx
	mulx	2*8(ap), s0, s0		C a3 a2
	C 0

	xor	R32(s1), R32(s1)

	mov	4*8(ap), %rdx
	mulx	1*8(ap), s5, s5		C a4 a1
	mulx	2*8(ap), s2, s3		C a4 a2
	mulx	3*8(ap), s6, s4		C a4 a3
	add	s5, s0
	adc	s1, s3
	add	s2, s0
	adc	s6, s3
	adc	s1, s4
	C 0, 3, 4

	mov	5*8(ap), %rdx
	mulx	0*8(ap), s8, s8		C a5 a0
	mulx	1*8(ap), s2, s6		C a5 a1
	mulx	2*8(ap), s5, s7		C a5 a2
	add	s2, s0
	adc	s5, s3
	mulx	3*8(ap), s9, s2		C a5 a3
	mulx	4*8(ap), %rdx, s5	C a5 a4
	adc	s9, s4
	adc	s1, s2
	add	s8, s0
	adc	s6, s3
	adc	s7, s4
	adc	%rdx, s2
	adc	s1, s5
	C 0, 3, 4, 2, 5

	mov	6*8(ap), %rdx
	test	%al, %al
	mulx	0*8(ap), s6, s7		C a6 a0
	mulx	1*8(ap), s8, s9		C a6 a1
	adcx	s6, s0
	adox	s7, s3
	adcx	s8, s3
	adox	s9, s4
	mulx	2*8(ap), s6, s7		C a6 a2
	mulx	3*8(ap), s8, s9		C a6 a3
	adcx	s6, s4
	adox	s7, s2
	adcx	s8, s2
	adox	s9, s5
	mulx	4*8(ap), s6, s7		C a6 a4
	mulx	5*8(ap), s8, s9		C a6 a5
	adox	s1, s7
	adc	s6, s5
	adc	s8, s7
	adc	s1, s9
	C 0, 3, 4, 2, 5, 7, 9

	add	s0, s0
	adc	s3, s3
	adc	s4, s4
	adc	s2, s2
	adc	s5, s5
	adc	s7, s7
	adc	s9, s9
	adc	R32(s1), R32(s1)
	C 0, 3, 4, 2, 5, 7, 9, 1

	mov	3*8(ap), %rdx
	mulx	%rdx, %rdx, s6		C a3^2
	add	%rdx, s0
	adc	s6, s3
	mov	4*8(ap), %rdx
	mulx	%rdx, %rdx, s6		C a4^2
	mov	s3, 0*8(rp)
	adc	%rdx, s4
	adc	s6, s2
	mov	5*8(ap), %rdx
	mulx	%rdx, s3, s6		C a5^2
	mov	s4, 1*8(rp)
	mov	s2, 2*8(rp)
	adc	s3, s5
	adc	s6, s7
	mov	6*8(ap), %rdx
	mulx	%rdx, s3, s6		C a6^2
	mov	s5, 3*8(rp)
	mov	s7, 4*8(rp)
	adc	s3, s9
	adc	s6, s1
	mov	s9, 5*8(rp)
	mov	s1, 6*8(rp)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_sqrhigh_8)
	C   0 1 2 3 4 5 6 7
	C 0
	C 1
	C 2
	C 3       e
	C 4     h x d
	C 5   h x x x d
	C 6 h x x x x x d
	C 7 x x x x x x x d
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	R32(s1), R32(s1)
	mov	4*8(ap), %rdx
	mulx	2*8(ap), s0, s0		C a4 a2
	mulx	3*8(ap), s4, s3		C a4 a3
	add	s4, s0
	adc	s1, s3
	C 0, 3

	mov	5*8(ap), %rdx
	mulx	1*8(ap), s2, s2		C a5 a1
	mulx	2*8(ap), s4, s5		C a5 a2
	mulx	3*8(ap), s6, s8		C a5 a3
	mulx	4*8(ap), s7, s9		C a5 a4
	add	s2, s0
	adc	s5, s3
	adc	s1, s8
	add	s4, s0
	adc	s6, s3
	adc	s7, s8
	adc	s1, s9
	C 0, 3, 8, 9

	mov	6*8(ap), %rdx
	test	%al, %al
	mulx	0*8(ap), s4, s4		C a6 a0
	mulx	1*8(ap), s5, s2		C a6 a1
	adox	s4, s0
	adcx	s5, s0
	adox	s2, s3
	mulx	2*8(ap), s4, s7		C a6 a2
	mulx	3*8(ap), s5, s2		C a6 a3
	adcx	s4, s3
	adox	s7, s8
	adcx	s5, s8
	adox	s2, s9
	mulx	4*8(ap), s4, s7		C a6 a4
	mulx	5*8(ap), s5, s2		C a6 a5
	adcx	s4, s9
	adox	s1, s7
	adcx	s5, s7
	adcx	s1, s2
	C 0, 3, 8, 9, 7, 2

	mov	7*8(ap), %rdx
	mulx	0*8(ap), s6, s5		C a7 a0
	mulx	1*8(ap), s10, s4	C a7 a1
	adcx	s6, s0
	adox	s5, s3
	adcx	s10, s3
	adox	s4, s8
	mulx	2*8(ap), s6, s5		C a7 a2
	mulx	3*8(ap), s10, s4	C a7 a3
	adcx	s6, s8
	adox	s5, s9
	mulx	4*8(ap), s6, s5		C a7 a4
	adcx	s10, s9
	adox	s4, s7
	adcx	s6, s7
	adox	s5, s2
	mulx	5*8(ap), s10, s4	C a7 a5
	mulx	6*8(ap), s6, s5		C a7 a6
	adox	s1, s4
	adc	s10, s2
	adc	s6, s4
	adc	s1, s5
	C 0, 3, 8, 9, 7, 2, 4, 5

	add	s0, s0
	adc	s3, s3
	adc	s8, s8
	adc	s9, s9
	adc	s7, s7
	adc	s2, s2
	adc	s4, s4
	adc	s5, s5
	adc	R32(s1), R32(s1)
	C 0, 3, 8, 9, 7, 2, 4, 5, 1

	mov	3*8(ap), %rdx
	mulx	%rdx, s6, s6		C a3^2
	add	s6, s0
	mov	4*8(ap), %rdx
	mulx	%rdx, s10, s6		C a4^2
	adc	s10, s3
	adc	s6, s8
	mov	5*8(ap), %rdx
	mulx	%rdx, s10, s6		C a5^2
	mov	s3, 0*8(rp)
	mov	s8, 1*8(rp)
	adc	s10, s9
	adc	s6, s7
	mov	6*8(ap), %rdx
	mulx	%rdx, s10, s6		C a6^2
	mov	s9, 2*8(rp)
	mov	s7, 3*8(rp)
	adc	s10, s2
	adc	s6, s4
	mov	7*8(ap), %rdx
	mulx	%rdx, s10, s6		C a7^2
	mov	s2, 4*8(rp)
	mov	s4, 5*8(rp)
	adc	s10, s5
	adc	s6, s1
	mov	s5, 6*8(rp)
	mov	s1, 7*8(rp)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()
