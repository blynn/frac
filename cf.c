// Demand channels. See squint paper by McIlroy.
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

struct channel_s {
  void *data;
  struct channel_s *next;
};
typedef struct channel_s channel_t[1];
typedef struct channel_s *channel_ptr;

struct cf_s {
  // Each continued fraction is a separate thread.
  pthread_t thread;

  // Signal 'demand' to compute next term.
  // 'ack' acknowledges that the signal was received,
  // to prevent a new 'demand' signal while processing a previous one.
  pthread_cond_t demand, ack;
  pthread_mutex_t demand_mu, ack_mu;

  // Demand channel.
  pthread_cond_t read_cond;
  pthread_mutex_t chan_mu;
  channel_ptr chan, next;

  int quitflag;
  void *data;
};

typedef struct cf_s *cf_t;

void cf_wait(cf_t cf) {
  // Wait for 'demand' signal.
  pthread_cond_wait(&cf->demand, &cf->demand_mu);
  // Acknowledge it, allowing future 'demand' signals.
  pthread_mutex_lock(&cf->ack_mu);
  pthread_cond_signal(&cf->ack);
  pthread_mutex_unlock(&cf->ack_mu);
  if (cf->quitflag) {
    pthread_mutex_unlock(&cf->demand_mu);
    pthread_exit(NULL);
  }
}

void cf_signal(cf_t cf) {
  pthread_mutex_lock(&cf->demand_mu);
  pthread_cond_signal(&cf->demand);
  pthread_mutex_unlock(&cf->demand_mu);
  pthread_cond_wait(&cf->ack, &cf->ack_mu);
}

void cf_free(cf_t cf) {
  cf->quitflag = 1;
  pthread_mutex_lock(&cf->demand_mu);
  pthread_cond_signal(&cf->demand);
  pthread_mutex_unlock(&cf->demand_mu);
  pthread_cond_wait(&cf->ack, &cf->ack_mu);
  pthread_join(cf->thread, NULL);
  pthread_mutex_destroy(&cf->demand_mu);
  pthread_mutex_destroy(&cf->ack_mu);
  pthread_cond_destroy(&cf->demand);
  pthread_cond_destroy(&cf->ack);
  free(cf);
}

void cf_put(cf_t cf, int i) {
  channel_ptr cnew = malloc(sizeof(*cnew));
  cnew->data = (void *) i;
  cnew->next = NULL;
  pthread_mutex_lock(&cf->chan_mu);
  if (cf->chan) {
    cf->next->next = cnew;
  } else {
    // Channel is empty so send signal in case someone is waiting for it.
    cf->chan = cnew;
    pthread_cond_signal(&cf->read_cond);
  }
  cf->next = cnew;
  pthread_mutex_unlock(&cf->chan_mu);
}

int cf_get(cf_t cf) {
  pthread_mutex_lock(&cf->chan_mu);
  // Wait for signal if channel is empty.
  if (!cf->chan) {
    pthread_cond_wait(&cf->read_cond, &cf->chan_mu);
  }
  channel_ptr c = cf->chan;
  int i = (int) c->data;
  cf->chan = c->next;
  free(c);
  pthread_mutex_unlock(&cf->chan_mu);
  return i;
}

cf_t cf_new(void *(*func)(cf_t)) {
  cf_t cf = malloc(sizeof(*cf));
  pthread_attr_t attr;
  pthread_cond_init(&cf->demand, NULL);
  pthread_mutex_init(&cf->demand_mu, NULL);
  pthread_cond_init(&cf->ack, NULL);
  pthread_mutex_init(&cf->ack_mu, NULL);
  pthread_cond_init(&cf->read_cond, NULL);
  pthread_mutex_init(&cf->chan_mu, NULL);
  pthread_mutex_lock(&cf->demand_mu);
  pthread_mutex_lock(&cf->ack_mu);
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  cf->chan = NULL;
  cf->next = NULL;
  cf->quitflag = 0;

  pthread_create(&cf->thread, &attr, (void*(*)(void *)) func, cf);
  pthread_attr_destroy(&attr);
  return cf;
}

static void *sqrt2(cf_t cf) {
  cf_put(cf, 1);
  for (;;) {
    cf_wait(cf);
    cf_put(cf, 2);
  }
  return NULL;
}

cf_t cf_new_sqrt2() {
  return cf_new(sqrt2);
}

int main() {
  cf_t a;
  a = cf_new_sqrt2();
  for (int i = 0; i < 10; i++) {
    cf_signal(a);
    printf("a: %d\n", cf_get(a));
  }
  cf_free(a);
  return 0;
}
