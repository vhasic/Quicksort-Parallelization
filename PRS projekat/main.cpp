#include "quickSortVersions.cpp"
#include <ctime>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <vector>
#include <omp.h>


/// Upis razultata u datoteku filename u csv formatu broj niti, vrijeme izvršavanja
void writeResultInFile(int brojNiti, int trajanje, std::string filename = "rezultati.csv") {
    std::ofstream izlazni_tok(filename, std::ios_base::app);
    izlazni_tok << brojNiti << "," << trajanje << std::endl;
}

/// Provjerava da li je niz pravilno sortiran u rastućem poretku
bool isSorted(const int *a, int n) {
    for (int i = 1; i < n; i++) {
        if (a[i - 1] > a[i]) {
            return false;
        }
    }
    return true;
}

void showArray(int* array, int first, int last){
    for (int i = 0; i<last-first;i++) {
        std::cout<<array[i]<<",";
    }
    std::cout<<std::endl;
}

/// Inicijalizira niz na najgori slučaj za QuickSort algoritam, niz već sortiran
int* initializeArrayWorstCase(int *array, int n){
    for (int i = 0; i < n; i++) {
//        array[i] = i+1;
        array[i] = n-i;
    }
    return array;
}

/// Inicijalizira niz na prosječan slučaj, generisanjem radnom elemenata između 1 i n
int* initializeArrayAverageCase(int *array, int n){
    for (int i = 0; i < n; i++) {
        array[i] = getRandomInt(1,n);
    }
    return array;
}

int main() {
    clock_t timeBefore, timeAfter;
    int executionTime;
    int n = 10000;
    int *array = new int[n];

//    array= initializeArrayWorstCase(array,n);
//    array= initializeArrayAverageCase(),n);

    //sekvencijalna verzija sa standardnom particijom
    array= initializeArrayWorstCase(array,n);
    timeBefore = std::clock();
    sequentialQuickSort(array, 0, n - 1);
    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::printf("Vrijeme izvrsenja: %d ms\n", executionTime);
    if (isSorted(array, n)) {
        std::printf("OK\n");
        writeResultInFile(1, executionTime);
    }
    else {
        std::printf("Nije OK\n");
    }

    //sekvencijalna verzija sa random pivotom u particiji
    array= initializeArrayWorstCase(array,n);
    timeBefore = std::clock();
    sequentialQuickSortRandomPivot(array, 0, n - 1);
    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::printf("Vrijeme izvrsenja: %d ms\n", executionTime);
    // provjera da li je niz dobro sortiran
    if (isSorted(array, n)) {
        std::printf("OK\n");
        writeResultInFile(1, executionTime);
    }
    else {
        std::printf("Nije OK\n");
    }

    //sekvencijalna verzija sa median of three pivotom u particiji
    array= initializeArrayWorstCase(array,n);
    timeBefore = std::clock();
    sequentialQuickSortMedianOfThreePivot(array, 0, n - 1);
    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::printf("Vrijeme izvrsenja: %d ms\n", executionTime);
    // provjera da li je niz dobro sortiran
    if (isSorted(array, n)) {
        std::printf("OK\n");
        writeResultInFile(1, executionTime);
    }
    else {
        std::printf("Nije OK\n");
    }


    // paralelizacija sa sekcijama
    array= initializeArrayWorstCase(array,n);
    timeBefore = std::clock();
    quickSortSections(array, 0, n - 1,100);
    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::printf("Vrijeme izvrsenja: %d ms\n", executionTime);
    writeResultInFile(4, executionTime);
    if (isSorted(array, n)) {
        std::printf("OK\n");
    }
    else {
        std::printf("Nije OK\n");
    }


    // paralelizacija sa taskovima
    array= initializeArrayWorstCase(array,n);
    std::vector<int> arrayNumberOfThreads({4});
    for (int numThreads: arrayNumberOfThreads) {
        timeBefore = std::clock();
#pragma omp parallel default(none) shared(array,n) num_threads(numThreads)
        {
#pragma omp single nowait
            {
                quickSortTasks(array, 0, n - 1, 100);
            }
        }
        timeAfter = std::clock();

        executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
        std::printf("Vrijeme izvrsenja: %d ms\n", executionTime);
        writeResultInFile(numThreads, executionTime);

        // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
        if (isSorted(array, n)) {
            std::printf("OK\n");
        }
        else {
            std::printf("Nije OK\n");
        }
    }

    // način da se dobije broj hardverski podržanih niti na računaru uz pomoć open mp
    //int numThreads = omp_get_num_procs();
    // std::cout << "Broj niti: " << numThreads << std::endl;


    delete[] array;
    return 0;
}