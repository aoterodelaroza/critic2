set(SOURCES blas.f dasum.f daxpy.f dbdsqr.f dgebd2.f dgebrd.f
            dgeco.f dgedi.f dgeev.f dgefa.f dgelq2.f dgelqf.f dgeqr2.f dgeqrf.f
            dgesvd.f disnan.f dlabrd.f dlaisnan.f dlas2.f dlasq1.f dlasq2.f
            dlasq3.f dlasq4.f dlasq5.f dlasq6.f dlasr.f dlasrt.f dlasv2.f dorgbr.f
            dorgl2.f dorglq.f dorm2r.f dormbr.f dorml2.f dormlq.f dormqr.f
            dp1vlu.f dpcoef.f dpolft.f dscal.f dswap.f eispack.f idamax.f lubksb.f
            ludcmp.f tqlgrm.f trace.f tredig.f)

add_library(oldlibs OBJECT ${SOURCES})
