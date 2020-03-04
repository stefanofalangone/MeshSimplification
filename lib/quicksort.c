#include <stdio.h>

enum { SIZE = 10 };

void quicksort(float *target, int left, int right) {
  if(left >= right) return;
  int i = left, j = right;
  float tmp, pivot = target[i];
  for(;;) {
    while(target[i] < pivot) i++;
    while(pivot < target[j]) j--;
    if(i >= j) break;
    tmp = target[i]; target[i] = target[j]; target[j] = tmp;
    i++; j--;
  }
  quicksort(target, left, i-1);
  quicksort(target, j+1, right);
}

int main() {
  int i;
  float array[SIZE] = { 7.7 , 2, 6, 3, 8, 5, 4, 1, 9, 7 };

  quicksort(array, 0, SIZE-1);

  for(i=0; i<SIZE; i++)
    printf("%f ", array[i]);
  printf("\n");
}
