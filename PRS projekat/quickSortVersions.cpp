#include "partitions.cpp"
#include <omp.h>
#include <x86intrin.h>
#include <immintrin.h>
#include <cstdint>
#include <cassert>
#include <iostream>

#define FORCE_INLINE inline __attribute__((always_inline))

template<class T>
inline void Log(const __m256i &value) {
    const size_t n = sizeof(__m256i) / sizeof(T);
    T buffer[n];
    _mm256_storeu_si256((__m256i *) buffer, value);
    for (int i = 0; i < n; i++)
        std::cout << buffer[i] << " ";
}

void sequentialQuickSort(int *array, int first, int last) {
    if (first < last) {
        int j = standardPartition(array, first, last);
        sequentialQuickSort(array, first, j - 1);
        sequentialQuickSort(array, j + 1, last);
    }
}

void sequentialQuickSortRandomPivot(int *array, int first, int last) {
    if (first < last) {
        int j = partition_randomPivot(array, first, last);
        sequentialQuickSort(array, first, j - 1);
        sequentialQuickSort(array, j + 1, last);
    }
}

void sequentialQuickSortMedianOfThreePivot(int *array, int first, int last) {
    if (first < last) {
        int j = partition_medianOfThreePivot(array, first, last);
        sequentialQuickSort(array, first, j - 1);
        sequentialQuickSort(array, j + 1, last);
    }
}

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

void quicksort_32(uint32_t *array, int left, int right, int sequentialLimit) {

    int i = left;
    int j = right;

    const uint32_t pivot = array[(i + j) / 2];
    const int AVX2_REGISTER_SIZE = 8; // in 32-bit words

    if (j - i >= sequentialLimit) {
        partition_epi32(array, pivot, i, j);
    } else {
        scalar_partition_epi32(array, pivot, i, j);
    }

    if (left < j) {
        quicksort_32(array, left, j, sequentialLimit);
    }

    if (i < right) {
        quicksort_32(array, i, right, sequentialLimit);
    }

}


//*******************************************************************************************************************

void quickSortTasks(int *array, int first, int last, int sequentialLimit) {
    if (first < last) {                                               // sekvencijalno sortiranje
        if (last - first < sequentialLimit) {
            return sequentialQuickSort(array, first, last);
        } else {                                                       // inače paralelno
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


void quickSortSections(int *array, int first, int last, int sequentialLimit, int numThreads = 4) {
    if (first < last) {                                               // sekvencijalno sortiranje
        if (last - first < sequentialLimit) {
            return sequentialQuickSort(array, first, last);
        } else {                                                       // inače paralelno
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
void quickSortTasks_v2(int *array, int first, int last, int sequentialLimit) {
    if (first < last) {                                               // sekvencijalno sortiranje
        if (last - first < sequentialLimit) {
            return sequentialQuickSort(array, first, last);
        } else {                                                       // inače paralelno
            int j = partition_medianOfThreePivot(array, first, last);
//            int j = standardPartition(array, first, last);

#pragma omp parallel default(none) shared(array) firstprivate(first, last, sequentialLimit, j)
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