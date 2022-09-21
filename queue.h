#include "lib.h"

#ifndef QUEUE_H
#define QUEUE_H

typedef struct node{
    double* data;
    struct node* prev;
    struct node* next;
} node;

typedef struct queue{
    int len;
    node* head;
    node* rear;
} queue;

static node* new_node(double* data);
queue* create_queue(void);
void queue_delete(void *self);
bool queue_is_empty(const queue *self);
bool push_queue(queue *self, double* data);
double* pop_queue(queue* self);

#endif
