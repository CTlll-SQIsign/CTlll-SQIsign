CFLAGS = -DNDEBUG

build/run: src/run.c src/inputs.h src/bkz_constants.h  src/read_data.h build/libbkz.a build/libcthelpers.a build/libctintbig.a build/libinputs.a
	gcc $(CFLAGS) src/run.c -L./build -linputs -lbkz -lcthelpers -lctintbig -lgmp  -o build/run

build/bench: src/bench.c src/inputs.h src/bench.h src/bkz_constants.h  src/read_data.h build/libbkz.a build/liblll.a  build/libcthelpers.a build/libctintbig.a build/libinputs.a
	gcc $(CFLAGS) src/bench.c -L./build -linputs -lbkz -llll -lcthelpers -lctintbig -lgmp  -o build/bench

build/small_test: src/small_test.c build/libtesttesthelpers.a build/libtesthelpers.a build/libnctintbig.a build/libbkz.a build/liblll.a  build/libcthelpers.a build/libctintbig.a build/libtestctintbig.a build/libtestbkz.a build/libtestlll.a build/libinputs.a
	gcc $(CFLAGS) src/small_test.c -L./build -linputs  -ltesttesthelpers -ltesthelpers -ltestctintbig -ltestbkz -ltestlll -lbkz -llll -lcthelpers -lctintbig -lnctintbig -lgmp -o build/small_test

build/test: src/test.c src/inputs.h src/bkz_constants.h build/libtesthelpers.a build/libnctintbig.a build/libbkz.a build/liblll.a build/libcthelpers.a build/libctintbig.a build/libinputs.a src/limbnum.h src/read_data.h
	gcc $(CFLAGS) src/test.c -L./build -ltesthelpers  -lbkz -linputs -llll -lcthelpers -lctintbig -lnctintbig -lgmp -o build/test

build/quality_test: src/test.c src/inputs.h src/bkz_constants.h build/libtesthelpers.a build/libnctintbig.a build/libbkz.a build/liblll.a build/libcthelpers.a build/libctintbig.a build/libinputs.a src/limbnum.h src/read_data.h
	gcc $(CFLAGS) -DPRINT_FLAG=1 src/test.c -L./build  -linputs-ltesthelpers  -lbkz -llll -lcthelpers -lctintbig -lnctintbig -lgmp -o build/quality_test

build/bench_stats: src/bench_stats.c src/inputs.h src/bench.h src/bkz_constants.h build/libtesthelpers.a build/libnctintbig.a build/libbkz.a build/liblll.a build/libcthelpers.a build/libctintbig.a build/libinputs.a src/limbnum.h src/read_data.h
	gcc $(CFLAGS) src/bench_stats.c -L./build -linputs  -ltesthelpers  -lbkz  -llll -lcthelpers -lctintbig -lnctintbig -lgmp -o build/bench_stats

build/bench_stats_lll: src/bench_stats_lll.c src/inputs.h src/bench.h src/bkz_constants.h build/libtesthelpers.a build/libnctintbig.a build/libbkz.a build/liblll.a build/libcthelpers.a build/libctintbig.a build/libinputs.a src/limbnum.h src/read_data.h
	gcc $(CFLAGS) src/bench_stats_lll.c -L./build -linputs  -ltesthelpers  -lbkz  -llll -lcthelpers -lctintbig -lnctintbig -lgmp -o build/bench_stats_lll

build/libctgrind.so: src/ctgrind.c
	gcc -o build/libctgrind.so -shared src/ctgrind.c -Wall -std=c99 -fPIC

build/libctgrind.so.1: build/libctgrind.so
	ln -s build/libctgrind.so build/libctgrind.so.1

build/bench_ctgrind: src/bench_ctgrind.c src/inputs.h src/ctgrind.h src/bkz_constants.h build/libctgrind.so.1 build/libbkz.a build/libcthelpers.a build/libctintbig.a build/libinputs.a src/limbnum.h src/read_data.h
	gcc $(CFLAGS) -Wall -ggdb  -std=c99  -Wextra src/bench_ctgrind.c -L./build -linputs -lbkz -lcthelpers -lctintbig -lgmp -lctgrind -o build/bench_ctgrind

run: build/run
	./build/run

test: build/test
	./build/test

bench: build/bench
	./build/bench

bench_stats: build/bench_stats
	./build/bench_stats

bench_stats_lll: build/bench_stats_lll
	./build/bench_stats_lll

bench_ctgrind: build/libctgrind.so build/libctgrind.so.1 build/bench_ctgrind
	LD_LIBRARY_PATH="$(CURDIR)/build/" valgrind -s --track-origins=yes --leak-check=full --show-leak-kinds=all --verbose --log-file=ctgrind_log.log ./build/bench_ctgrind

small_test: build/small_test
	./build/small_test

quality_test: build/quality_test


build/libinputs.a: build/inputs.o
	ar -rcs build/libinputs.a  build/inputs.o

build/inputs.o: src/inputs.h src/inputs.c src/ct_intbig/ct_intbig.h src/limbnum.h
	gcc -c $(CFLAGS) src/inputs.c -o build/inputs.o

src/inputs.c: 
	touch src/inputs.c


build/libtesttesthelpers.a: build/test_test_helpers.o
	ar -rcs build/libtesttesthelpers.a  build/test_test_helpers.o

build/test_test_helpers.o: src/test_helpers/test_test_helpers.c src/test_helpers/test_test_helpers.h src/test_helpers/test_helpers.h src/test_helpers/nct_helpers.h src/nct_intbig/nct_intbig.h  src/limbnum.h
	gcc -c $(CFLAGS) src/test_helpers/test_test_helpers.c -o build/test_test_helpers.o


build/libtesthelpers.a: build/test_code.o build/test_helpers.o build/nct_helpers.o
	ar -rcs build/libtesthelpers.a  build/test_code.o build/test_helpers.o build/nct_helpers.o

build/test_code.o: src/test_helpers/test_code.c src/test_helpers/test_code.h src/test_helpers/test_helpers.h src/test_helpers/nct_helpers.h src/nct_intbig/nct_intbig.h src/limbnum.h
	gcc -c $(CFLAGS) src/test_helpers/test_code.c -o build/test_code.o

build/test_helpers.o: src/test_helpers/test_helpers.c src/test_helpers/test_helpers.h src/test_helpers/nct_helpers.h src/nct_intbig/nct_intbig.h src/limbnum.h
	gcc -c $(CFLAGS) src/test_helpers/test_helpers.c -o build/test_helpers.o

build/nct_helpers.o: src/test_helpers/nct_helpers.c src/test_helpers/nct_helpers.h src/nct_intbig/nct_intbig.h src/limbnum.h
	gcc -c $(CFLAGS) src/test_helpers/nct_helpers.c -o build/nct_helpers.o


build/libnctintbig.a: build/nct_intbig.o
	ar -rcs build/libnctintbig.a build/nct_intbig.o

build/nct_intbig.o: src/nct_intbig/nct_intbig.c src/nct_intbig/nct_intbig.h src/limbnum.h
	gcc -c $(CFLAGS) src/nct_intbig/nct_intbig.c -o build/nct_intbig.o



build/libtestlll.a: build/test_lll.o
	ar -rcs build/libtestlll.a build/test_lll.o

build/test_lll.o: src/lll/test_lll.c src/lll/test_lll.h src/lll/lll.h src/ct_intbig/ct_intbig.h src/limbnum.h
	gcc -c $(CFLAGS) src/lll/test_lll.c -o build/test_lll.o


build/liblll.a: build/lll.o build/lll_constants.o
	ar -rcs build/liblll.a build/lll.o build/lll_constants.o

build/lll.o: src/lll/lll.c src/lll/lll.h src/ct_intbig/ct_intbig.h src/limbnum.h
	gcc -c $(CFLAGS) src/lll/lll.c -o build/lll.o

build/lll_constants.o: src/lll/lll_constants.c src/lll/lll_constants.h src/ct_intbig/ct_intbig.h src/limbnum.h
	gcc -c $(CFLAGS) src/lll/lll_constants.c -o build/lll_constants.o



build/libtestbkz.a: build/test_bkz.o
	ar -rcs build/libtestbkz.a build/test_bkz.o

build/test_bkz.o: src/bkz/test_bkz.c src/bkz/test_bkz.h src/bkz/bkz.h src/ct_helpers/ct_helpers.h src/ct_intbig/ct_intbig.h src/limbnum.h
	gcc -c $(CFLAGS) src/bkz/test_bkz.c -o build/test_bkz.o


build/libbkz.a: build/bkz.o
	ar -rcs build/libbkz.a build/bkz.o

build/bkz.o: src/bkz/bkz.c src/bkz/bkz.h src/ct_helpers/ct_helpers.h src/ct_intbig/ct_intbig.h src/limbnum.h
	gcc -c $(CFLAGS) src/bkz/bkz.c -o build/bkz.o



build/libcthelpers.a: build/ct_helpers.o
	ar -rcs build/libcthelpers.a build/ct_helpers.o

build/ct_helpers.o: src/ct_helpers/ct_helpers.c src/ct_helpers/ct_helpers.h src/ct_intbig/ct_intbig.h src/limbnum.h
	gcc -c $(CFLAGS) src/ct_helpers/ct_helpers.c -o build/ct_helpers.o



build/libtestctintbig.a: build/test_ibz.o build/test_ibq.o
	ar -rcs build/libtestctintbig.a build/test_ibz.o build/test_ibq.o

build/test_ibz.o: src/ct_intbig/test_ibz.c src/ct_intbig/ct_intbig.h src/ct_intbig/test_ct_intbig.h src/limbnum.h
	gcc -c $(CFLAGS) src/ct_intbig/test_ibz.c -o build/test_ibz.o

build/test_ibq.o: src/ct_intbig/test_ibq.c src/ct_intbig/ct_intbig.h src/ct_intbig/test_ct_intbig.h src/limbnum.h
	gcc -c $(CFLAGS) src/ct_intbig/test_ibq.c -o build/test_ibq.o


build/libctintbig.a: build/ct_intbig.o
	ar -rcs build/libctintbig.a build/ct_intbig.o

build/ct_intbig.o: src/ct_intbig/ct_intbig.c src/ct_intbig/ct_intbig.h src/limbnum.h
	gcc -c $(CFLAGS) src/ct_intbig/ct_intbig.c -o build/ct_intbig.o
