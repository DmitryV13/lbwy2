#include <assert.h>
#include <malloc.h>
#include <stdio.h>

typedef enum {
    false,
    true
} bool;

typedef struct {
    int val;
    struct node *next;
} nodeOfList;

typedef struct {
    nodeOfList *node;
} abstract_data_t;

nodeOfList* newCell() {
    nodeOfList *newC = malloc(sizeof(nodeOfList));
    newC->val = 0;
    newC->next = NULL;
    return newC;
};

abstract_data_t make_empty() {
    abstract_data_t new;
    new.node=NULL;
    return new;
};

size_t length(nodeOfList * list){
    if(list== NULL)
        return 0;
    int l=0;
    while(list){
        l++;
        list=list->next;
    }
    return l;
}

abstract_data_t make_copy(const abstract_data_t *initObject) {
    if(length(initObject->node)==0){
        return make_empty();
    }
    abstract_data_t newObject;
    newObject.node=newCell();
    nodeOfList *temp = newObject.node;
    nodeOfList *initNode=initObject->node;
    while (initNode->next) {
        temp->val = initNode->val;
        temp->next = newCell();
        temp = temp->next;
        initNode = initNode->next;
    }
    temp->val = initNode->val;
    return newObject;
};

abstract_data_t make_from_array(const int arr[], size_t length) {
    if(length==0){
        return make_empty();
    }
    abstract_data_t newObject;
    newObject.node=newCell();
    nodeOfList *temp = newObject.node;
    temp->val = arr[0];
    for (int i = 1; i < length; i++) {
        temp->next = newCell();
        temp = temp->next;
        temp->val = arr[i];
    }
    return newObject;
};

abstract_data_t make_of_size(size_t length) {
    if(length==0){
        return make_empty();
    }
    abstract_data_t newObject;
    newObject.node=newCell();
    nodeOfList *temp = newObject.node;
    for (int i = 1; i < length; i++) {
        temp->next = newCell();
        temp = temp->next;
    }
    return newObject;
};

void clearNode(nodeOfList* target){
    if (target->next == NULL) {
        free(target);
    } else {
        clearNode(target->next);
        target->next=NULL;
        free(target);
    }
}
void clear(abstract_data_t *object) {
    if (length(object->node) == 0) return;
    clearNode(object->node);
    object->node=NULL;
};

int normalization(int ind, int n){
    if(ind>=n) return n-1;
    if(ind<0) return n+ind<0?0:n+ind;
    return ind;
}

int at(const abstract_data_t *target, int index){
    index= normalization(index, length(target->node));
    nodeOfList *temp=target->node;
    for (int i = 0; i < index; ++i) {
        temp=temp->next;
    }
    return temp->val;
}

bool is_equal(const abstract_data_t *object, const abstract_data_t *target) {
    if(length(object->node)!= length(target->node)) return false;
    nodeOfList *objectTemp=object->node;
    nodeOfList *targetTemp=target->node;
    while (objectTemp) {
        if (objectTemp->val != targetTemp->val)
            return false;
        objectTemp = objectTemp->next;
        targetTemp = targetTemp->next;
    }
    return true;
};

void output(abstract_data_t *d) {
    nodeOfList *temp=d->node;
    while(temp){
        printf("%d",temp->val);
        temp=temp->next;
    }
}

int main(int argc, char const *argv[]) {
    int array[] = {1, 4, 7, 9};

    abstract_data_t a = make_from_array(array, 4);

    assert(9 == at(&a, 3));
    assert(9 == at(&a, -1));
    assert(9 == at(&a, 13));
    assert(1 == at(&a, -13));

    abstract_data_t b = make_copy(&a);
    assert(is_equal(&a, &b));

    clear(&a);
    clear(&b);
}
