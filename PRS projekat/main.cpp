#include "partitions.cpp"
#include <ctime>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <omp.h>

/// Upis razultata u datoteku filename u csv formatu broj niti, vrijeme izvršavanja
void writeResultInFile(std::string filename, int brojNiti, int trajanje) {
    std::ofstream izlazni_tok(filename, std::ios_base::app);
    izlazni_tok << brojNiti << "," << trajanje << std::endl;
}

/// Provjerava da li je niz sortiran pravilno u rastućem poretku
bool isSorted(const int* a, int n) {
    for (int i = 1; i < n; i++) {
        if (a[i - 1] > a[i]) {
            return false;
        }
    }
    return true;
}


void sequentialQuickSort(int* array, int first, int last) {
    if (first < last) {
        int j = standardPartition(array, first, last);
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

int main() {
    //todo: neka ovakva inicijalizacija bude samo privremeno dok ne skontamo kako radi openMP
    // poslije bi trebali napraviti klasu sa različitim metodama za inicijalizaciju:
    // double, int, najgori slučaj (opadajuće sortiran), prosječan slučaj (random elementi)...
    int n = 100000;
    int* array1 = new int[n];
    int* array2 = new int[n];
        //for (int i = 0; i < n; i++) {
        //    array1[i] = n - i;
        //    array2[i] = n - i;
        //}
    for (int i = 0; i < n; i++) {
        int randomInt = getRandomInt();
        array1[i] = randomInt;
        array2[i] = randomInt;
    }
    // todo riješiti problem sa stack overflow-om

    clock_t timeBefore = std::clock();
    sequentialQuickSort(array1, 0, n - 1);
    clock_t timeAfter = std::clock();

    int executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    //std::cout << "Vrijeme izvrsenja: " << executionTime << " ms." << std::endl;
    std::printf("Vrijeme izvrsenja: %d ms\n", executionTime);
    writeResultInFile("rezultati.csv", 1, executionTime);

    // paralelizacija sa taskovima
    timeBefore = std::clock();
#pragma omp parallel default(none) shared(array2,n) num_threads(4)
    {
#pragma omp single nowait
        {
            quickSortTasks(array2, 0, n - 1, 100);
        }
    }
    timeAfter = std::clock();

    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    //std::cout << "Vrijeme izvrsenja: " << executionTime << " ms." << std::endl;
    std::printf("Vrijeme izvrsenja: %d ms\n", executionTime);
    writeResultInFile("rezultati.csv", 4, executionTime);

    // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
    if (isSorted(array2, n)) {
        //std::cout << "OK" << std::endl;
        std::printf("OK\n");
    }
    else {
        //std::cout << "Nije OK" << std::endl;
        std::printf("Nije OK\n");
    }

    // način da se dobije broj hardverski podržanih niti na računaru uz pomoć open mp
    //int numThreads = omp_get_num_procs();
    // std::cout << "Broj niti: " << numThreads << std::endl;


    delete[] array1;
    delete[] array2;
    return 0;
}