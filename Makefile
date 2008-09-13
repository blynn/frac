CF_FILES:=cf.c cf_mobius.c cf_famous.c bihom.c

target : pi pi2 sqrt2 epow

% : %.c $(CF_FILES)
	#gcc -O3 -std=c99 -Wall -o $@ $^ -lgmp -lpthread
	gcc -g -std=c99 -Wall -o $@ $^ -lgmp -lpthread

test :
	gcc -g -std=c99 -Wall -o cf_test cf_test.c $(CF_FILES) -lgmp -lpthread
	gcc -g -std=c99 -Wall -o famous_test famous_test.c $(CF_FILES) -lgmp -lpthread
