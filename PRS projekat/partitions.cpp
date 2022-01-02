#include <algorithm>
#include <random>
#include <functional>
#include <cassert>
#include <omp.h>
#include <x86intrin.h>
#include <iostream>

#define FORCE_INLINE inline __attribute__((always_inline))

/// Vraća random generisani integer između min i max
int getRandomInt(int min = 1, int max = 1000) {
    // generate thread safe random number
    //std::cout<<"Uso"<<std::endl;
    std::uniform_int_distribution<int> dice_distribution(min, max);
    std::mt19937 random_number_engine; // pseudorandom number generator
    auto dice_roller = std::bind(dice_distribution, random_number_engine);
    int random_roll = dice_roller();  // Generate one of the integers [min,max]
    //std::cout<<"Random roll" <<random_roll;
    return random_roll;
}

int getRandomInt2(int min = 1, int max = 1000){
    return rand() % max;
}

///Lomuto's partition Scheme
int standardPartition(int* array, int low, int high) {
    int pivot = array[high];
    // Index of smaller element
    int i = low - 1;

    for(int j = low; j <= high-1; j++){
        // If current element is smaller than or equal to pivot
        if(array[j] <= pivot){
            // increment index of smaller element
            i++;
            std::swap(array[i], array[j]);
        }
    }
    std::swap(array[i+1], array[high]);
    return i+1;
}

int standardPartition32(uint32_t * array, int low, int high) {
    int pivot = array[high];
    // Index of smaller element
    int i = low - 1;

    for(int j = low; j <= high-1; j++){
        // If current element is smaller than or equal to pivot
        if(array[j] <= pivot){
            // increment index of smaller element
            i++;
            std::swap(array[i], array[j]);
        }
    }
    std::swap(array[i+1], array[high]);
    return i+1;
}


/// Generates Random Pivot, swaps pivot with end element and calls the partition function
/// O(1.386 * n * log n)
int partition_randomPivot(int* array, int low, int high) {
    int random_roll = getRandomInt(low,high);

    std::swap(array[random_roll], array[high]);
    return standardPartition(array, low, high);
}

/// Optimizirana particija, tako da se za pivot uzima Median-of-three O(1.188*n*log n)
///Median-of-three: https://en.wikipedia.org/wiki/Quicksort#:~:text=Median%2Dof%2Dthree%20code%20snippet%20for%20Lomuto%20partition%3A
int partition_medianOfThreePivot(int* array, int low, int high) {
    int mid = (low + high) / 2;
    if (array[mid] < array[low]) {
        std::swap(array[low], array[mid]);
    }
    if (array[high] < array[low]) {
        std::swap(array[low], array[high]);
    }
    if (array[mid] < array[high]) {
        std::swap(array[mid], array[high]);
    }
    return standardPartition(array, low, high);
}

int partition_medianOfThreePivot32(uint32_t * array, int low, int high) {
    int mid = (low + high) / 2;
    if (array[mid] < array[low]) {
        std::swap(array[low], array[mid]);
    }
    if (array[high] < array[low]) {
        std::swap(array[low], array[high]);
    }
    if (array[mid] < array[high]) {
        std::swap(array[mid], array[high]);
    }
    return standardPartition32(array, low, high);
}

void scalar_partition_epi32(uint32_t *array, uint32_t pivot, int &left, int &right) {

    while (left <= right) {

        while (array[left] < pivot) {
            left += 1;
        }

        while (array[right] > pivot) {
            right -= 1;
        }

        if (left <= right) {
            const uint32_t t = array[left];
            array[left]      = array[right];
            array[right]     = t;

            left  += 1;
            right -= 1;
        }
    }
}

__m256i FORCE_INLINE bitmask_to_bytemask_epi32(uint8_t bm) {

const __m256i mask = _mm256_set1_epi32(bm);
const __m256i bits = _mm256_setr_epi32(0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80);
const __m256i tmp  = _mm256_and_si256(mask, bits);

return _mm256_cmpeq_epi32(tmp, bits);
}


void FORCE_INLINE align_masks(uint8_t& a, uint8_t& b, uint8_t& rem_a, uint8_t& rem_b, __m256i& shuffle_a, __m256i& shuffle_b) {

assert(a != 0);
assert(b != 0);

uint8_t tmpA = a;
uint8_t tmpB = b;

uint32_t __attribute__((__aligned__(32))) tmpshufa[8];
uint32_t __attribute__((__aligned__(32))) tmpshufb[8];

while (tmpA != 0 && tmpB != 0) {
int idx_a = __builtin_ctz(tmpA);
int idx_b = __builtin_ctz(tmpB);

tmpA = tmpA & (tmpA - 1);
tmpB = tmpB & (tmpB - 1);

tmpshufb[idx_a] = idx_b;
tmpshufa[idx_b] = idx_a;
}

a = a ^ tmpA;
b = b ^ tmpB;

assert(a != 0);
assert(b != 0);
assert(_mm_popcnt_u64(a) == _mm_popcnt_u64(b));

rem_a = tmpA;
rem_b = tmpB;

shuffle_a = _mm256_load_si256((__m256i*)tmpshufa);
shuffle_b = _mm256_load_si256((__m256i*)tmpshufb);
}


__m256i FORCE_INLINE merge(const __m256i mask, const __m256i a, const __m256i b) {

    return _mm256_or_si256(
            _mm256_and_si256(mask, a),
            _mm256_andnot_si256(mask, b)
    );
}


void FORCE_INLINE swap_epi32(
        __m256i& a, __m256i& b,
        uint8_t mask_a, const __m256i shuffle_a,
        uint8_t mask_b, const __m256i shuffle_b) {

const __m256i to_swap_b = _mm256_permutevar8x32_epi32(a, shuffle_a);
const __m256i to_swap_a = _mm256_permutevar8x32_epi32(b, shuffle_b);
const __m256i ma    = bitmask_to_bytemask_epi32(mask_a);
const __m256i mb    = bitmask_to_bytemask_epi32(mask_b);
//a = merge(ma, to_swap_a, a);
//b = merge(mb, to_swap_b, b);
a = merge(ma, to_swap_a, a);
b = merge(mb, to_swap_b, b);
}



#define _mm256_iszero(vec) (_mm256_testz_si256(vec, vec) != 0)

void FORCE_INLINE partition_epi32(uint32_t* array, const uint32_t pv, int& left, int& right) {

const int N = 8; // the number of items in a register (256/32)

__m256i L;
__m256i R;
uint8_t maskL = 0;
uint8_t maskR = 0;

const __m256i pivot = _mm256_set1_epi32(pv);

int origL = left;
int origR = right;

while (true) {

if (maskL == 0) {
while (true) {
if (right - (left + N) + 1 < 2*N) {
goto end;
}

L = _mm256_loadu_si256((__m256i*)(array + left));
const __m256i bytemask = _mm256_cmpgt_epi32(pivot, L);

if (_mm256_testc_ps((__m256)bytemask, (__m256)_mm256_set1_epi32(-1))) {
left += N;
} else {
maskL = ~_mm256_movemask_ps((__m256)bytemask);
break;
}
}

}

if (maskR == 0) {
while (true) {
if ((right - N) - left + 1 < 2*N) {
goto end;
}

R = _mm256_loadu_si256((__m256i*)(array + right - N + 1));
const __m256i bytemask = _mm256_cmpgt_epi32(pivot, R);
if (_mm256_iszero(bytemask)) {
right -= N;
} else {
maskR = _mm256_movemask_ps((__m256)bytemask);
break;
}
}

}

assert(left <= right);
assert(maskL != 0);
assert(maskR != 0);

uint8_t mL;
uint8_t mR;
__m256i shuffleL;
__m256i shuffleR;

align_masks(maskL, maskR, mL, mR, shuffleL, shuffleR);
swap_epi32(L, R, maskL, shuffleL, maskR, shuffleR);

maskL = mL;
maskR = mR;

if (maskL == 0) {
_mm256_storeu_si256((__m256i*)(array + left), L);
left += N;
}

if (maskR == 0) {
_mm256_storeu_si256((__m256i*)(array + right - N + 1), R);
right -= N;
}

} // while

end:

assert(!(maskL != 0 && maskR != 0));

if (maskL != 0) {
_mm256_storeu_si256((__m256i*)(array + left), L);
} else if (maskR != 0) {
_mm256_storeu_si256((__m256i*)(array + right - N + 1), R);
}

if (left < right) {
int less    = 0;
int greater = 0;
const int all = right - left + 1;

for (int i=left; i <= right; i++) {
less    += int(array[i] < pv);
greater += int(array[i] > pv);
}

if (all == less) {
// all elements in range [left, right] less than pivot
scalar_partition_epi32(array, pv, origL, left);
} else if (all == greater) {
// all elements in range [left, right] greater than pivot
scalar_partition_epi32(array, pv, left, origR);
} else {
scalar_partition_epi32(array, pv, left, right);
}
}
}
