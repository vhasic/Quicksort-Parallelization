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

template<typename Tip>
void sequentialQuickSort(Tip *array, int first, int last) {
    if (first < last) {
        int j = standardPartition(array, first, last);
        sequentialQuickSort(array, first, j - 1);
        sequentialQuickSort(array, j + 1, last);
    }
}

template<typename Tip>
void sequentialQuickSortRandomPivot(Tip *array, int first, int last) {
    if (first < last) {
        int j = partition_randomPivot(array, first, last);
        sequentialQuickSortRandomPivot(array, first, j - 1);
        sequentialQuickSortRandomPivot(array, j + 1, last);
    }
}

template<typename Tip>
void sequentialQuickSortMedianOfThreePivot(Tip *array, int first, int last) {
    if (first < last) {
        int j = partition_medianOfThreePivot(array, first, last);
        sequentialQuickSortMedianOfThreePivot(array, first, j - 1);
        sequentialQuickSortMedianOfThreePivot(array, j + 1, last);
    }
}



//**********************************************************************************************************************
//32 bitni integer brojevi

void quicksort_32(uint32_t *array, int left, int right, int sequentialLimit) {

    int i = left;
    int j = right;

    const uint32_t pivot = array[(i + j) / 2];

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
template<typename Tip>
void quickSortTasks(Tip *array, int first, int last, int sequentialLimit) {
    if (first < last) {                                               // sekvencijalno sortiranje
        if (last - first < sequentialLimit) {
            return sequentialQuickSortMedianOfThreePivot(array, first, last);
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

/*void quickSortTasksAVX(uint32_t *array, int left, int right, int sequentialLimit) {
    // Median of three pivot
    int mid = (left + right) / 2;
    if (array[mid] < array[left]) {
        std::swap(array[left], array[mid]);
    }
    if (array[right] < array[left]) {
        std::swap(array[left], array[right]);
    }
    if (array[mid] < array[right]) {
        std::swap(array[mid], array[right]);
    }
    uint32_t pivot = array[right];

    int i = left;
    int j = right;

    // AVX particija
    if (j - i >= sequentialLimit) {
        partition_epi32(array, pivot, i, j);
    } else {
        scalar_partition_epi32(array, pivot, i, j);
    }

    // paralelizacija sa taskovima
    if (left < j) {
#pragma omp task default(none) firstprivate(array, left, j, sequentialLimit)
        {
            quicksort_32(array, left, j, sequentialLimit);
            //std::printf("Hello from thread %d of %d \n", omp_get_thread_num(), omp_get_num_threads());
        }
    }

    if (i < right) {
#pragma omp task default(none) firstprivate(array, right, i, sequentialLimit)
        {
            quicksort_32(array, i, right, sequentialLimit);
            //std::printf("Hello from thread %d of %d \n", omp_get_thread_num(), omp_get_num_threads());
        }
    }
}*/

void quickSortTasksAVX(uint32_t *array, int left, int right, int sequentialLimit) {
    // ovo je bolje za average case ali lošije za worst case
    // Median of three pivot
    int mid = (left + right) / 2;
    if (array[right] < array[left]) {
        std::swap(array[left], array[right]);
    }
    if (array[mid] < array[left]) {
        std::swap(array[left], array[mid]);
    }
    if (array[right] < array[mid]) {
        std::swap(array[mid], array[right]);
    }
    uint32_t pivot = array[mid];

    int i = left;
    int j = right;
//    ovo je bolje za worst case ali je lošije za average case
//    const uint32_t pivot = array[(i + j) / 2];

    // AVX particija
    if (j - i >= sequentialLimit) {
        partition_epi32(array, pivot, i, j);
    } else {
        scalar_partition_epi32(array, pivot, i, j);
    }

    // paralelizacija sa taskovima
    if (left < j) {
#pragma omp task default(none) firstprivate(array, left, j, sequentialLimit)
        {
            quickSortTasksAVX(array, left, j, sequentialLimit);
            //std::printf("Hello from thread %d of %d \n", omp_get_thread_num(), omp_get_num_threads());
        }
    }

    if (i < right) {
#pragma omp task default(none) firstprivate(array, right, i, sequentialLimit)
        {
            quickSortTasksAVX(array, i, right, sequentialLimit);
            //std::printf("Hello from thread %d of %d \n", omp_get_thread_num(), omp_get_num_threads());
        }
    }
}

template<typename Tip>
void quickSortSections(Tip *array, int first, int last, int sequentialLimit, int numThreads = 4) {
    if (first < last) {                                               // sekvencijalno sortiranje
        if (last - first < sequentialLimit) {
            return sequentialQuickSortMedianOfThreePivot(array, first, last);
        } else {                                                       // inače paralelno
            int j = partition_medianOfThreePivot(array, first, last);

#pragma omp parallel sections default(none) firstprivate(array, first, last, j, sequentialLimit) num_threads(numThreads)
            {
#pragma omp section
                {
                    quickSortSections(array, first, j - 1, sequentialLimit);
                }
#pragma omp section
                {
                    quickSortSections(array, j + 1, last, sequentialLimit);
                }
            }
        }
    }
}