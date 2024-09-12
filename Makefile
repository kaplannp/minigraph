CXX=g++
CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -std=c99 -msse4 -O3
CPPFLAGS= -g -Wall -std=c++20
INCLUDES=
OBJS=		kalloc.o kthread.o algo.o sys.o gfa-base.o gfa-io.o gfa-aug.o gfa-bbl.o gfa-ed.o \
            sketch.o misc.o bseq.o options.o shortk.o miniwfa.o \
			index.o lchain.o gchain1.o galign.o gcmisc.o map-algo.o cal_cov.o \
			format.o gmap.o ggsimple.o ggen.o asm-call.o serialization.cpp.o
DUMP_OBJS= kalloc.dump.o kthread.dump.o algo.dump.o sys.dump.o gfa-base.dump.o gfa-io.dump.o gfa-aug.dump.o gfa-bbl.dump.o gfa-ed.dump.o sketch.dump.o misc.dump.o bseq.dump.o options.dump.o shortk.dump.o miniwfa.dump.o index.dump.o lchain.dump.o gchain1.dump.o galign.dump.o gcmisc.dump.o map-algo.dump.o cal_cov.dump.o format.dump.o gmap.dump.o ggsimple.dump.o ggen.dump.o asm-call.dump.o serialization.cpp.dump.o deserialization.cpp.dump.o loadParams.cpp.dump.o

PROG=		minigraph
LIBS=		-lz -lpthread -lm

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

.SUFFIXES:.c .o
.PHONY:all clean depend

all:$(PROG)

serialization.cpp.o: serialization.cpp
	$(CXX) -c $(CPPFLAGS) $(INCLUDES) serialization.cpp -o serialization.cpp.o

%.cpp.dump.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(INCLUDES) $< -o $@ -DWRITE_KERNEL_INS_AND_OUTS -DLOAD_INS_FROM_DUMP_TO_LOAD

.c.o:
		$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@ -DWRITE_KERNEL_INS_AND_OUTS

%.dump.o: %.c
		$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@ -DWRITE_KERNEL_INS_AND_OUTS -DLOAD_INS_FROM_DUMP_TO_LOAD

minigraph:$(OBJS) main.o
		$(CXX) $(CFLAGS) $^ -o $@ $(LIBS) -DWRITE_KERNEL_INS_AND_OUTS

#To use this you name a file in this directory DumpToLoad
minigraph.load:$(DUMP_OBJS) main.o
		$(CXX) $(CFLAGS) $^ -o $@ $(LIBS) -DWRITE_KERNEL_INS_AND_OUTS -DLOAD_INS_FROM_DUMP_TO_LOAD

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM
		rm -f minigraph.load

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

algo.o: kalloc.h algo.h miniwfa.h kvec-km.h ksort.h
asm-call.o: mgpriv.h minigraph.h gfa.h ggen.h bseq.h gfa-priv.h algo.h
bseq.o: bseq.h kvec-km.h kalloc.h kseq.h
cal_cov.o: mgpriv.h minigraph.h gfa.h gfa-priv.h algo.h kalloc.h
format.o: kalloc.h mgpriv.h minigraph.h gfa.h
galign.o: mgpriv.h minigraph.h gfa.h kalloc.h miniwfa.h
gchain1.o: mgpriv.h minigraph.h gfa.h ksort.h khashl.h kalloc.h gfa-priv.h serialization.h
gcmisc.o: mgpriv.h minigraph.h gfa.h kalloc.h
gfa-aug.o: gfa-priv.h gfa.h ksort.h
gfa-base.o: gfa-priv.h gfa.h kstring.h khashl.h kalloc.h ksort.h
gfa-bbl.o: gfa-priv.h gfa.h kalloc.h ksort.h kvec.h
gfa-ed.o: gfa-priv.h gfa.h kalloc.h ksort.h khashl.h kdq.h kvec-km.h
gfa-io.o: kstring.h gfa-priv.h gfa.h kseq.h
ggen.o: kthread.h kalloc.h sys.h bseq.h ggen.h minigraph.h gfa.h mgpriv.h
ggen.o: gfa-priv.h
ggsimple.o: mgpriv.h minigraph.h gfa.h gfa-priv.h kalloc.h bseq.h algo.h
ggsimple.o: sys.h ggen.h kvec-km.h
gmap.o: kthread.h kalloc.h bseq.h sys.h mgpriv.h minigraph.h gfa.h gfa-priv.h
index.o: mgpriv.h minigraph.h gfa.h khashl.h kalloc.h kthread.h kvec-km.h
index.o: sys.h
kalloc.o: kalloc.h
kthread.o: kthread.h
lchain.o: mgpriv.h minigraph.h gfa.h kalloc.h krmq.h
main.o: mgpriv.h minigraph.h gfa.h gfa-priv.h sys.h ketopt.h
map-algo.o: kalloc.h mgpriv.h minigraph.h gfa.h khashl.h ksort.h
miniwfa.o: miniwfa.h kalloc.h
misc.o: mgpriv.h minigraph.h gfa.h ksort.h
options.o: mgpriv.h minigraph.h gfa.h sys.h
shortk.o: mgpriv.h minigraph.h gfa.h ksort.h kavl.h algo.h khashl.h kalloc.h
sketch.o: kvec-km.h kalloc.h mgpriv.h minigraph.h gfa.h
sys.o: sys.h
