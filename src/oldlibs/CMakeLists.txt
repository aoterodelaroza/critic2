set(CMAKE_BUILD_TYPE Release)

set(SOURCES daxpy.f dbdsqr.f dcopy.f ddot.f dgebak.f dgebal.f dgebd2.f dgebrd.f
  dgeev.f dgehd2.f dgehrd.f dgelq2.f dgelqf.f dgemm.f dgemv.f dgeqr2.f dgeqrf.f
  dger.f dgesvd.f dgetrf2.f dgetrf.f dgetri.f dhseqr.f disnan.f dlabad.f dlabrd.f
  dlacpy.f dladiv.f dlae2.f dlaev2.f dlaexc.f dlahqr.f dlahr2.f dlaisnan.f dlaln2.f
  dlamch.f dlange.f dlanst.f dlansy.f dlanv2.f dlapy2.f dlaqr0.f dlaqr1.f dlaqr2.f
  dlaqr3.f dlaqr4.f dlaqr5.f dlarfb.f dlarf.f dlarfg.f dlarft.f dlarfx.f dlartg.f
  dlas2.f dlascl.f dlaset.f dlasq1.f dlasq2.f dlasq3.f dlasq4.f dlasq5.f dlasq6.f
  dlasr.f dlasrt.f dlassq.f dlasv2.f dlaswp.f dlasy2.f dlasyf.f dlatrd.f dnrm2.f
  dorg2l.f dorg2r.f dorgbr.f dorghr.f dorgl2.f dorglq.f dorgql.f dorgqr.f dorgtr.f
  dorm2r.f dormbr.f dormhr.f dorml2.f dormlq.f dormqr.f drot.f dscal.f dsteqr.f
  dsterf.f dswap.f dsyev.f dsymv.f dsyr2.f dsyr2k.f dsyr.f dsytd2.f dsytf2.f
  dsytrd.f dsytrf.f dsytri.f dtrevc3.f dtrexc.f dtrmm.f dtrmv.f dtrsm.f dtrti2.f
  dtrtri.f idamax.f ieeeck.f iladlc.f iladlr.f ilaenv.f iparmq.f lsame.f xerbla.f
  )
add_library(oldlibs OBJECT ${SOURCES})
