#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <immintrin.h>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
  for(int i=0; i<N; i++) {
    /*
    for(int j=0; j<N; j++) {
      if(i != j) {
        float rx = x[i] - x[j];
        float ry = y[i] - y[j];
        float r = std::sqrt(rx * rx + ry * ry);
        fx[i] -= rx * m[j] / (r * r * r);
        fy[i] -= ry * m[j] / (r * r * r);
      }
    }
    */
    float a[N];
    __m256 zerovec = _mm256_setzero_ps();
    __m256 onevec = _mm256_set1_ps(1);

    __m256 xivec = _mm256_set1_ps(x[i]);
    __m256 yivec = _mm256_set1_ps(y[i]);
    __m256 xvec = _mm256_load_ps(x);
    __m256 yvec = _mm256_load_ps(y);
    
    __m256 mask = _mm256_cmp_ps(xivec, xvec, _CMP_EQ_OQ);
    
    // rx ry
    __m256 rxvec = _mm256_sub_ps(xivec, xvec);
    __m256 ryvec = _mm256_sub_ps(yivec, yvec);

    // 1/r
    __m256 onervec1 = _mm256_mul_ps(rxvec, rxvec);
    __m256 onervec2 = _mm256_mul_ps(ryvec, ryvec);
    __m256 onervec = _mm256_add_ps(onervec1, onervec2);
    onervec = _mm256_blendv_ps(onervec, onevec, mask);
    onervec = _mm256_rsqrt_ps(onervec);

    __m256 mvec = _mm256_load_ps(m);
    
    // fx[i]
    __m256 fxivec = _mm256_mul_ps(rxvec, mvec);
    fxivec = _mm256_mul_ps(fxivec, onervec);
    fxivec = _mm256_mul_ps(fxivec, onervec);
    fxivec = _mm256_mul_ps(fxivec, onervec);
    fxivec = _mm256_blendv_ps(fxivec, zerovec, mask);

    // fy[i]
    __m256 fyivec = _mm256_mul_ps(ryvec, mvec);
    fyivec = _mm256_mul_ps(fyivec, onervec);
    fyivec = _mm256_mul_ps(fyivec, onervec);
    fyivec = _mm256_mul_ps(fyivec, onervec);
    fyivec = _mm256_blendv_ps(fyivec, zerovec, mask);

    // sum of fx[i]
    __m256 bvec = _mm256_permute2f128_ps(fxivec, fxivec, 1);
    bvec = _mm256_add_ps(bvec, fxivec);
    bvec = _mm256_hadd_ps(bvec, bvec);
    bvec = _mm256_hadd_ps(bvec, bvec);
    _mm256_store_ps(a, bvec);
    fx[i] -= a[0]

    // sum of fy[i]
    bvec = _mm256_permute2f128_ps(fyivec, fyivec, 1);
    bvec = _mm256_add_ps(bvec, fyivec);
    bvec = _mm256_hadd_ps(bvec, bvec);
    bvec = _mm256_hadd_ps(bvec, bvec);
    _mm256_store_ps(a, bvec);
    fy[i] -= a[0]

    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
