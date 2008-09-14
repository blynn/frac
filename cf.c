// Demand channels. See squint paper by McIlroy.
//
// TODO: Handle messy thread problems. What happens if a thread quits
// but then another tries to signal and read its channel?
// TODO: What if the continued fraction terminates?
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <gmp.h>

struct channel_s {
  void *data;
  struct channel_s *next;
};
typedef struct channel_s channel_t[1];
typedef struct channel_s *channel_ptr;

struct cf_s {
  // Each continued fraction is a separate thread.
  pthread_t thread;
  // Helgrind prints warnings for these condition variables.
  // Rewrite with semaphores?
  // When queue is empty, and there is demand for the next term.
  sem_t demand_sem;
  // When the queue was empty, and we just added to it.
  pthread_cond_t read_cond;
  pthread_mutex_t chan_mu;
  channel_ptr chan, next;

  int quitflag;
  void *data;
};

typedef struct cf_s *cf_t;

void *cf_data(cf_t cf) {
  return cf->data;
}

// A bit like cooperative multitasking. Continued fractions are expected
// to call this as often as practical, and on a return value of 0,
// to drop everything and stop.
int cf_wait(cf_t cf) {
  for (;;) {
    sem_wait(&cf->demand_sem);
    // The wait is over!
    if (cf->quitflag) {
      return 0;
    }
    pthread_mutex_lock(&cf->chan_mu);
    // ... but we keep waiting unless the channel is empty.
    if (!cf->chan) break;
    pthread_mutex_unlock(&cf->chan_mu);
    // The channel could be emptied in the meantime, but that
    // implies at least one sem_post() call, so we'll notice next iteration.
  }
  pthread_mutex_unlock(&cf->chan_mu);
  return 1;
}

void cf_free(cf_t cf) {
  // These two statements force a thread out of its next/current cf_wait.
  cf->quitflag = 1;
  sem_post(&cf->demand_sem);

  pthread_join(cf->thread, NULL);
  pthread_mutex_lock(&cf->chan_mu);
  channel_ptr c = cf->chan;
  while (c) {
    channel_ptr cnext = c->next;
    free(c->data);
    free(c);
    c = cnext;
  }
  pthread_mutex_unlock(&cf->chan_mu);
  sem_destroy(&cf->demand_sem);
  free(cf);
}

void cf_put(cf_t cf, mpz_t z) {
  // TODO: Block or something if there's a large backlog on the queue.
  channel_ptr cnew = malloc(sizeof(*cnew));
  mpz_ptr znew = malloc(sizeof(*znew));
  mpz_init(znew);
  mpz_set(znew, z);
  cnew->data = znew;
  cnew->next = NULL;
  pthread_mutex_lock(&cf->chan_mu);
  if (cf->chan) {
    cf->next->next = cnew;
  } else {
    // Channel is empty. Now that we're populating it, send signal
    // in case someone is waiting for data.
    cf->chan = cnew;
    pthread_cond_signal(&cf->read_cond);
  }
  cf->next = cnew;
  pthread_mutex_unlock(&cf->chan_mu);
}

void cf_put_int(cf_t cf, int n) {
  mpz_t z;
  mpz_init(z);
  mpz_set_si(z, n);
  cf_put(cf, z);
  mpz_clear(z);
}

void cf_get(mpz_t z, cf_t cf) {
  pthread_mutex_lock(&cf->chan_mu);
  if (!cf->chan) {
    // If channel is empty, send demand signal and wait for read signal.
    sem_post(&cf->demand_sem);
    pthread_cond_wait(&cf->read_cond, &cf->chan_mu);
  }
  channel_ptr c = cf->chan;
  cf->chan = c->next;
  pthread_mutex_unlock(&cf->chan_mu);
  mpz_ptr znew = c->data;
  mpz_set(z, znew);
  mpz_clear(znew);
  free(c->data);
  free(c);
}

cf_t cf_new(void *(*func)(cf_t), void *data) {
  cf_t cf = malloc(sizeof(*cf));
  cf->chan = NULL;
  cf->next = NULL;
  cf->quitflag = 0;
  cf->data = data;
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_mutex_init(&cf->chan_mu, NULL);
  sem_init(&cf->demand_sem, 0, 0);
  pthread_cond_init(&cf->read_cond, NULL);
  pthread_create(&cf->thread, &attr, (void*(*)(void *)) func, cf);
  pthread_attr_destroy(&attr);
  return cf;
}
