// Tee: read input channel and write to two output channels

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "cf.h"

struct tee_s {
  int n;
  cf_t parent;
};
typedef struct tee_s tee_t[1];
typedef struct tee_s *tee_ptr;

static void *branch_loop(cf_t cf) {
  tee_ptr t = cf_data(cf);
  mpz_t z;
  mpz_init(z);
  while(cf_wait(cf)) {
    cf_put_int(t->parent, t->n);  // Causes parent to put to cf.
    cf_signal(t->parent);
  }
  cf_put_int(t->parent, -1 - t->n);  // Notify parent of our destruction.
  mpz_clear(z);
  free(t);
  return NULL;
}

struct parent_data_s {
  cf_t kid[2];
  cf_t in;
};
typedef struct parent_data_s parent_data_t[1];
typedef struct parent_data_s *parent_data_ptr;

// TODO: Destroy after children are destroyed.
static void *parent_loop(cf_t cf) {
  struct backlog_s {
    mpz_t z;
    struct backlog_s *next;
  };
  typedef struct backlog_s *backlog_ptr;

  backlog_ptr head[2];
  backlog_ptr last[2];
  head[0] = NULL;
  head[1] = NULL;
  last[0] = NULL;
  last[1] = NULL;
  mpz_t z;
  mpz_init(z);
  parent_data_ptr pd = cf_data(cf);
  for (;;) {
    cf_wait_special(cf);
    cf_get(z, cf);
    int k = mpz_get_ui(z);
    if (k < 0) {
      k = -k + 1;
      pd->kid[k] = NULL;
      backlog_ptr pnext = head[k], p;
      do {
	p = pnext;
	pnext = p->next;
	mpz_clear(p->z);
	free(p);
      } while(pnext);
      if (!pd->kid[1-k]) break;
    } else {
      if (head[k]) {
	backlog_ptr p = head[k];
	mpz_set(z, p->z);
	mpz_clear(p->z);
	head[k] = p->next;
	free(p);
	if (!head[k]) last[k] = NULL;
      } else {
	cf_get(z, pd->in);
	if (pd->kid[1-k]) {
	  backlog_ptr p = malloc(sizeof(*p));
	  mpz_init(p->z);
	  mpz_set(p->z, z);
	  p->next = NULL;
	  if (last[1-k]) last[1-k]->next = p;
	  last[1-k] = p;
	  if (!head[1-k]) head[1-k] = p;
	}
      }
    }
    cf_put(pd->kid[k], z);
  }
  mpz_clear(z);
  return NULL;
}

cf_t cf_new_parent(cf_t in) {
  parent_data_ptr p = malloc(sizeof(*p));
  p->in = in;
  p->kid[0] = NULL;
  p->kid[1] = NULL;
  return cf_new(parent_loop, p);
}

cf_t cf_new_branch(cf_t parent, int n) {
  tee_ptr t = malloc(sizeof(*t));
  t->n = n;
  t->parent = parent;
  parent_data_ptr p = cf_data(parent);
  cf_t res = cf_new(branch_loop, t);
  p->kid[n] = res;
  return res;
}

void cf_tee(cf_t *out_array, cf_t in) {
  cf_t parent = cf_new_parent(in);
  out_array[0] = cf_new_branch(parent, 0);
  out_array[1] = cf_new_branch(parent, 1);
}
