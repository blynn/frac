CF_FILES:=cf.c mobius.c famous.c bihom.c taylor.c
TEST_FILES:=$(CF_FILES) test.c

target : pi pi2 sqrt2 epow

% : %.c $(CF_FILES)
	#gcc -O3 -std=c99 -Wall -o $@ $^ -lgmp -lpthread
	gcc -g -std=c99 -Wall -o $@ $^ -lgmp -lpthread

test :
	gcc -g -std=c99 -Wall -o cf_test cf_test.c $(TEST_FILES) -lgmp -lpthread
	gcc -g -std=c99 -Wall -o famous_test famous_test.c $(TEST_FILES) -lgmp -lpthread
	gcc -g -std=c99 -Wall -o bihom_test bihom_test.c $(TEST_FILES) -lgmp -lpthread
