#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
struct list{
  void *value;
  struct list *next;
};
typedef struct list List;

List *addList(List *l, void *t){
  List *ret;
  if(l==NULL){
    ret=(List *)malloc(sizeof(List));
    ret->value=t;
    ret->next=NULL;
  }else{
    ret=(List *)malloc(sizeof(List));
    ret->value=t;
    ret->next=l;
  }
  return ret;
}

int contains(List *l, void *c){
  while(l!=NULL){
    if(l->value==c) return 1;
    l=l->next;
  }
  return 0;
}

List *intersectList(List *l1, List *l2){
  List *curr1=l1;
  List *intersectionList=NULL;
  while(curr1!=NULL){
    if(contains(l2, curr1->value)>0 ){
      intersectionList=addList(intersectionList, curr1->value);
    }
    curr1=curr1->next;
  }
  return intersectionList;
}


void deallocList(List *l){
  List *next;
  while(l!=NULL){
    next=l->next;
    free(l);
    l=next;
  }
}

List *removeList(List *l, void *c){
  List *ret;
  List *head=l, *pred;
  if(l==NULL) return NULL;
  if(l->value==c){
    ret=l->next;
    free(l);
    return ret;
  }else{
    while(l!=NULL && l->value!=c){
      pred=l;
      l=l->next;
    }
    if(l!=NULL){
      pred->next=l->next;
      free(l);
    }
    return head;
  }
}

void printPointerList(List *l){
  int c=0;
  printf("List: ");
  while(l!=NULL){
    c++;
    printf("%p -> ", l);
    l=l->next;
    if(c%10==0) sleep(1);
  }
  printf("\\ \n");
}

void printPointerValue(List *l){
  int c=0;
  printf("List: ");
  while(l!=NULL){
    c++;
    printf("%p -> ", l->value);
    l=l->next;
    if(c%10==0) sleep(1);
  }
  printf("\\ \n");
}

int sizeList(List *l){
  int ret=0;
  while(l!=NULL){
    ret++;
    l=l->next;
  }
  return ret;
}


List *duplicateElem(List *l){
  List *ret=NULL;
  if(l!=NULL){
    ret=addList(NULL, l->value);
    ret->next=l->next;
  }
  return ret;
}

List *unionList(List *l1, List *l2){
  if(l1==NULL) return l2;
  if(l2==NULL) return l1;
  List *ret=NULL;
  List *curr;
  int isHead=1;
  while(l1!=NULL){
    ret=addList(ret, l1->value);
    l1=l1->next;
  }
  List *last=ret;
  while(l2!=NULL){
    if(contains(ret, l2->value)==0){
      ret=addList(ret, l2->value);
    }
      l2=l2->next;
  }
  return last;
}
