dnl
dnl Copyright (C) 2024 Albin Ahlbäck
dnl
dnl This file is part of FLINT.
dnl
dnl FLINT is free software: you can redistribute it and/or modify it under
dnl the terms of the GNU Lesser General Public License (LGPL) as published
dnl by the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.  See <https://www.gnu.org/licenses/>.
dnl

include(`config.m4')

define(`rp', `%rdi')
define(`ap', `%rsi')
define(`bp', `%rdx')
define(`n', `%rcx')

define(`s0', `%r8')
define(`s1', `%r9')
define(`s2', `%r10')
define(`s3', `%r11')

define(`sx', `%rax')

dnl NOTE: This function requires n >= 4

dnl NOTE: This function could easily not be inplace without pushing registers,
dnl       but currently I do not know if this function is going to be used on
dnl       other functions than multiplications.

dnl NOTE: Make this function shorter.

	TEXT

	ALIGN(16)
PROLOGUE(flint_mpn_2add_n_inplace)
	mov	R32(n), R32(sx)
	lea	L(tab)(%rip), s0
	shr	$3, R32(n)
	and	$7, R32(sx)	C Also resets OF and CF

ifdef(`PIC',
`	movslq	(s0,sx,4), sx
	lea	(sx,s0), s0
	xor	R32(sx), R32(sx)
	jmp	*s0
',`
	jmp	*(s0,sx,8)
')

L(p1):	mov	0*8(rp), s3
	mov	1*8(rp), s0
	mov	2*8(rp), s1
	adcx	0*8(ap), s3
	adox	0*8(bp), s3
	adcx	1*8(ap), s0
	adox	1*8(bp), s0
	adcx	2*8(ap), s1
	lea	1*8(rp), rp
	lea	1*8(bp), bp
	lea	1*8(ap), ap
	mov	s3, -1*8(rp)
	jmp	L(a0)

L(p0):	mov	0*8(rp), s0
	mov	1*8(rp), s1
	adcx	0*8(ap), s0
	adox	0*8(bp), s0
	adcx	1*8(ap), s1
	jmp	L(a0)

L(p5):	mov	0*8(rp), s3
	mov	1*8(rp), s0
	mov	2*8(rp), s1
	adcx	0*8(ap), s3
	adox	0*8(bp), s3
	adcx	1*8(ap), s0
	adox	1*8(bp), s0
	adcx	2*8(ap), s1
	lea	-3*8(rp), rp
	lea	-3*8(bp), bp
	lea	-3*8(ap), ap
	mov	s3, 3*8(rp)
	jmp	L(a4)

L(p4):	mov	0*8(rp), s0
	mov	1*8(rp), s1
	adcx	0*8(ap), s0
	adox	0*8(bp), s0
	adcx	1*8(ap), s1
	lea	-4*8(rp), rp
	lea	-4*8(bp), bp
	lea	-4*8(ap), ap
	jmp	L(a4)

L(p7):	mov	0*8(rp), s1
	mov	1*8(rp), s2
	mov	2*8(rp), s3
	adcx	0*8(ap), s1
	adox	0*8(bp), s1
	adcx	1*8(ap), s2
	adox	1*8(bp), s2
	adcx	2*8(ap), s3
	lea	-1*8(rp), rp
	lea	-1*8(bp), bp
	lea	-1*8(ap), ap
	mov	s1, 1*8(rp)
	jmp	L(a6)

L(p6):	mov	0*8(rp), s2
	mov	1*8(rp), s3
	adcx	0*8(ap), s2
	adox	0*8(bp), s2
	adcx	1*8(ap), s3
	lea	-2*8(rp), rp
	lea	-2*8(bp), bp
	lea	-2*8(ap), ap
	jmp	L(a6)

L(p3):	mov	0*8(rp), s1
	mov	1*8(rp), s2
	mov	2*8(rp), s3
	adcx	0*8(ap), s1
	adox	0*8(bp), s1
	adcx	1*8(ap), s2
	adox	1*8(bp), s2
	adcx	2*8(ap), s3
	lea	3*8(rp), rp
	lea	3*8(bp), bp
	lea	3*8(ap), ap
	mov	s1, -3*8(rp)
	jmp	L(a2)

L(p2):	mov	0*8(rp), s2
	mov	1*8(rp), s3
	adcx	0*8(ap), s2
	adox	0*8(bp), s2
	adcx	1*8(ap), s3
	lea	2*8(rp), rp
	lea	2*8(bp), bp
	lea	2*8(ap), ap
	C jmp	L(a2)

	C n = 12 -> n = 1, kx = 4
	C 2, 3, 4, 5

	ALIGN(32)
L(a2):	mov	0*8(rp), s0	C 01 start
	mov	1*8(rp), s1
	adox	-1*8(bp), s3
	adcx	0*8(ap), s0
	mov	s2, -2*8(rp)
	mov	s3, -1*8(rp)	C 23 end
	adox	0*8(bp), s0
	adcx	1*8(ap), s1
L(a0):	mov	2*8(rp), s2	C 23 start
	mov	3*8(rp), s3
	adox	1*8(bp), s1
	adcx	2*8(ap), s2
	lea	-1(n), R32(n)
	mov	s0, 0*8(rp)
	mov	s1, 1*8(rp)	C 01 end
	adox	2*8(bp), s2
	adcx	3*8(ap), s3
L(a6):	mov	4*8(rp), s0	C 01 start
	mov	5*8(rp), s1
	adox	3*8(bp), s3
	adcx	4*8(ap), s0
	mov	s2, 2*8(rp)
	mov	s3, 3*8(rp)	C 23 end
	adox	4*8(bp), s0
	adcx	5*8(ap), s1
L(a4):	mov	6*8(rp), s2	C 23 start
	mov	7*8(rp), s3
	adox	5*8(bp), s1
	adcx	6*8(ap), s2
	mov	s0, 4*8(rp)
	mov	s1, 5*8(rp)	C 01 end
	adox	6*8(bp), s2
	adcx	7*8(ap), s3
	jrcxz	L(end)
	lea	8*8(bp), bp
	lea	8*8(ap), ap
	lea	8*8(rp), rp
	jmp	L(a2)

L(end):	adox	7*8(bp), s3
	mov	s2, 6*8(rp)
	mov	s3, 7*8(rp)
	seto	R8(sx)		C sx < 8 prior, so it is contained in 8 bits
	adc	R32(n), R32(sx)	C n = 0

	ret
EPILOGUE()

	JUMPTABSECT
	ALIGN(8)
L(tab):	JMPENT(	L(p0), L(tab))
	JMPENT(	L(p1), L(tab))
	JMPENT(	L(p2), L(tab))
	JMPENT(	L(p3), L(tab))
	JMPENT(	L(p4), L(tab))
	JMPENT(	L(p5), L(tab))
	JMPENT(	L(p6), L(tab))
	JMPENT(	L(p7), L(tab))
	TEXT
