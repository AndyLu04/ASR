#include "queue.h"

static node* new_node(double* data)
{
    node* my_node = malloc(sizeof(node));
    my_node->data = data;
    my_node->prev = NULL;
    my_node->next = NULL;
    return my_node;
}

queue* create_queue(void)
{
    queue* my_q = malloc(sizeof(queue));
    my_q->len = 0;
    my_q->head = NULL;
    my_q->rear = NULL;
    return my_q;
}

void queue_delete(void *self)
{
    if (!self)
        return;
    node *curr = ((queue *) self)->head;
    while (curr) {
        node *temp = curr;
        curr = curr->next;
        free(temp);
    }
    free(self);
}

bool queue_is_empty(const queue *self)
{
    assert(self);
    return !(self->head) ? true : false;
}

bool push_queue(queue *self, double* data)
{
    assert(self);
    node* my_node = new_node(data);

    if (!(self->head)) {
        self->head = my_node;
        self->rear = my_node;
        return true;
    }

    self->rear->next = my_node;
    my_node->prev = self->rear;
    self->rear = my_node;
    self->len += 1;
    return true;
}

double* pop_queue(queue* self)
{
    assert(!queue_is_empty(self));

    double* popped = self->head->data;

    if (self->head == self->rear) {
        free(self->head);
        self->head = NULL;
        self->rear = NULL;
    }
    else {
        node *curr = self->head;
        self->head = curr->next;
        free(curr);
    }

    self->len -= 1;

    return popped;
}
