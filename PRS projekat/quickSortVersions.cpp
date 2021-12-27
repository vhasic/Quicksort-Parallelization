#include "partitions.cpp"


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
void quickSortTasks(int* array, int first, int last, int sequentialLimit) {
    if (first < last) {                                               // sekvencijalno sortiranje
        if (last - first < sequentialLimit) {
            return sequentialQuickSort(array, first, last);
        }
        else {                                                       // inače paralelno
//            int j = partition_medianOfThreePivot(array, first, last);
            int j = standardPartition(array, first, last);

#pragma omp task default(none) firstprivate(array, first, j, sequentialLimit)
            {
                quickSortTasks(array, first, j - 1, sequentialLimit);
                //                std::printf("Hello from thread %d of %d \n", omp_get_thread_num(), omp_get_num_threads());
            }
#pragma omp task default(none) firstprivate(array, last, j, sequentialLimit)
            {
                quickSortTasks(array, j + 1, last, sequentialLimit);
                //                std::printf("Hello from thread %d of %d \n", omp_get_thread_num(), omp_get_num_threads());
            }
        }
    }
}


void quickSortSections(int* array, int first, int last, int sequentialLimit) {
    if (first < last) {                                               // sekvencijalno sortiranje
        if (last - first < sequentialLimit) {
            return sequentialQuickSort(array, first, last);
        }
        else {                                                       // inače paralelno
            int j = standardPartition(array, first, last);

#pragma omp parallel sections default(none) firstprivate(array, first, last, j, sequentialLimit) num_threads(4)
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
