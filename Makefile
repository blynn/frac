CF_FILES:=cf.c cf_mobius.c

target : pi pi2 sqrt2

% : %.c $(CF_FILES)
	#gcc -O3 -std=c99 -Wall -o $@ $^ -lgmp -lpthread
	gcc -g -std=c99 -Wall -o $@ $^ -lgmp -lpthread

test :
	gcc -g -std=c99 -Wall -o cf_test cf_test.c cf.c -lgmp -lpthread
