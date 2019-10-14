set(SOURCES geom2.c geom.c global.c io.c libqhull.c mem.c merge.c
    poly2.c poly.c qset.c random.c rboxlib.c stat.c user.c usermem.c
    userprintf.c userprintf_rbox.c)

add_library(qhull OBJECT ${SOURCES})
