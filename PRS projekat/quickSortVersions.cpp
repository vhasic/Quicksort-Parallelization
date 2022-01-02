#include "partitions.cpp"
#include <omp.h>
#include <x86intrin.h>
#include <immintrin.h>
#include <cstdint>
#include <cassert>
#include <iostream>

#define FORCE_INLINE inline __attribute__((always_inline))

template<class T> inline void Log(const __m256i & value)
{
    const size_t n = sizeof(__m256i) / sizeof(T);
    T buffer[n];
    _mm256_storeu_si256((__m256i*)buffer, value);
    for (int i = 0; i < n; i++)
        std::cout << buffer[i] << " ";
}

void sequentialQuickSort(int* array, int first, int last) {
    if (first < last) {
        int j = standardPartition(array, first, last);
        sequentialQuickSort(array, first, j - 1);
        sequentialQuickSort(array, j + 1, last);
    }
}

void sequentialQuickSortRandomPivot(int* array, int first, int last) {
    if (first < last) {
        int j = partition_randomPivot(array, first, last);
        sequentialQuickSort(array, first, j - 1);
        sequentialQuickSort(array, j + 1, last);
    }
}

void sequentialQuickSortMedianOfThreePivot(int* array, int first, int last) {
    if (first < last) {
        int j = partition_medianOfThreePivot(array, first, last);
        sequentialQuickSort(array, first, j - 1);
        sequentialQuickSort(array, j + 1, last);
    }
}

//TODO: Uraditi particiju sa SIMD implementacijom

/**
 * Funkcija vrši paralelno sortiranje niza uz pomoć QuickSort algoritma
 * Ukoliko je niz premali (< granice) da se izvršava paralelno, jer bi se gubilo previše vremena na paralelnom overhead-u,
 * onda  se takav niz sortira sekvencijalno
 *
 * @param array Pokazivač na niz koji se sortira
 * @param first Donja granica niza koji se sortira
 * @param last Gornja granica niza koji se sortira
 * @param sequentialLimit Granica za sekvencijalno izvršavanje
 */

//**********************************************************************************************************************
//32 bitni integer brojevi

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

void FORCE_INLINE partition_epi32(uint32_t* array, int& left, int& right) {

    uint32_t pv = partition_medianOfThreePivot32(array, left, right);
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
            std::cout<<"MASKR"<<",";
            while (true) {
                std::cout<<"WHILE"<<",";
                if ((right - N) - left + 1 < 2*N) {
                    std::cout<<"END"<<",";
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


//*******************************************************************************************************************

void quickSortTasks(int* array, int first, int last, int sequentialLimit) {
    if (first < last) {                                               // sekvencijalno sortiranje
        if (last - first < sequentialLimit) {
            return sequentialQuickSort(array, first, last);
        }
        else {                                                       // inače paralelno
            int j = partition_medianOfThreePivot(array, first, last);
//            int j = standardPartition(array, first, last);

#pragma omp task default(none) firstprivate(array, first, j, sequentialLimit)
            {
                quickSortTasks(array, first, j - 1, sequentialLimit);
                //std::printf("Hello from thread %d of %d \n", omp_get_thread_num(), omp_get_num_threads());
            }
#pragma omp task default(none) firstprivate(array, last, j, sequentialLimit)
            {
                quickSortTasks(array, j + 1, last, sequentialLimit);
                //std::printf("Hello from thread %d of %d \n", omp_get_thread_num(), omp_get_num_threads());
            }
        }
    }
}


void quickSortSections(int* array, int first, int last, int sequentialLimit, int numThreads=4) {
    if (first < last) {                                               // sekvencijalno sortiranje
        if (last - first < sequentialLimit) {
            return sequentialQuickSort(array, first, last);
        }
        else {                                                       // inače paralelno
            int j = partition_medianOfThreePivot(array, first, last);
//            int j = standardPartition(array, first, last);

#pragma omp parallel sections default(none) firstprivate(array, first, last, j, sequentialLimit) num_threads(numThreads)
            {
#pragma omp section
                {
                    quickSortTasks(array, first, j - 1, sequentialLimit);
                }
#pragma omp section
                {
                    quickSortTasks(array, j + 1, last, sequentialLimit);
                }
            }
        }
    }
}

//Izvršavaju samo 2 niti, i na mom dvojezgrenom procesoru ima najbolje rezultate
void quickSortTasks_v2(int* array, int first, int last, int sequentialLimit) {
    if (first < last) {                                               // sekvencijalno sortiranje
        if (last - first < sequentialLimit) {
            return sequentialQuickSort(array, first, last);
        }
        else {                                                       // inače paralelno
            int j = partition_medianOfThreePivot(array, first, last);
//            int j = standardPartition(array, first, last);

#pragma omp parallel default(none) shared(array) firstprivate(first,last,sequentialLimit,j)
            {
#pragma omp single nowait
                {
#pragma omp task default(none) firstprivate(array, first, j, sequentialLimit)
                    {
                        quickSortTasks(array, first, j - 1, sequentialLimit);
                        //std::printf("Hello from thread %d of %d \n", omp_get_thread_num(), omp_get_num_threads());
                    }
#pragma omp task default(none) firstprivate(array, last, j, sequentialLimit)
                    {
                        quickSortTasks(array, j + 1, last, sequentialLimit);
                        //std::printf("Hello from thread %d of %d \n", omp_get_thread_num(), omp_get_num_threads());
                    }
                }
            }
        }
    }
}