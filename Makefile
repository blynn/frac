.PHONY: test target clean

CF_OBJS:=cf.o mobius.o famous.o bihom.o taylor.o test.o newton.o tee.o
TESTS:=bihom_test cf_test famous_test mobius_test newton_test tee_test
BINS:=pi hakmem

target : $(BINS)

libfrac.a : $(CF_OBJS)
	ar rvs $@ $^

%.o : %.c
	gcc -O3 -std=c99 -Wall -c -o $@ $<
	#gcc -g -std=c99 -Wall -c -o $@ $<

% : %.c libfrac.a
	gcc -O3 -std=c99 -Wall -o $@ $< -lgmp -lpthread -lfrac -L .
	#gcc -g -std=c99 -Wall -o $@ $< -lgmp -lpthread -lfrac -L .

test: $(TESTS)

clean:
	-rm *.o $(TESTS) $(BINS) libfrac.a
