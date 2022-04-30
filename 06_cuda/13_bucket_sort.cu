#include <cstdio>
#include <cstdlib>
#include <vector>


__global__ void bucket_sort(int *key, int n, int range) {
  int i = threadIdx.x;
  extern __shared__ int s[];
  int *bucket = s;
  int *offset = (int *)&bucket[range];
  int *b = (int *)&offset[range];

  if (i < range) bucket[i] = 0;
  atomicAdd(&bucket[key[i]], 1);

  if (i < range) {
    offset[i] = bucket[i];
    for (int j=1; j<range; j<<=1) {
      b[i] = offset[i];
      __syncthreads();
      if (i-j>=0) offset[i] += b[i-j];
      __syncthreads();
    }
  }

  if (i < range) {
      int j = 0;
      if(i>0) j = offset[i-1];
      for (; bucket[i]>0; bucket[i]--) {
        key[j++] = i;
      }
  }
}



int main() {
  int n = 50;
  int range = 5;
  // std::vector<int> key(n);

  int *key;
  cudaMallocManaged(&key, n*sizeof(int));

  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  /*
  std::vector<int> bucket(range);
  for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }
  for (int i=0; i<n; i++) {
    bucket[key[i]]++;
  }
  for (int i=0, j=0; i<range; i++) {
    for (; bucket[i]>0; bucket[i]--) {
      key[j++] = i;
    }
  }
  */

  bucket_sort<<<1,n, range*sizeof(int)*3>>>(key, n, range);
  cudaDeviceSynchronize();

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");

  cudaFree(key);

}