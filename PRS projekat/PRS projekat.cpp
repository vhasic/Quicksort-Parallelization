#include <algorithm>
#include <random>
#include <ctime>
#include <iostream>
#include <fstream>
#include <omp.h>

/// Vraća random generisani integer između min i max
int getRandomInt(int min = 1, int max = 1000) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> randomInt(min, max);

    return randomInt(rng);
}

/// Upis razultata u datoteku filename u csv formatu broj niti, vrijeme izvršavanja
void writeResultInFile(std::string filename, int brojNiti, int trajanje) {
    std::ofstream izlazni_tok(filename, std::ios_base::app);
    izlazni_tok << brojNiti << "," << trajanje << std::endl;
}

/// Provjerava da li je niz sortiran pravilno u rastućem poretku
bool isSorted(const int* a, int n)
{
    for (int i = 1; i < n; i++)
    {
        if (a[i - 1] > a[i])
        {
            return false;
        }
    }
    return true;
}

int standardPartition(int* array, int first, int last) {
    int pivot = array[first];
    int p = first + 1;
    while (p <= last) {
        if (array[p] <= pivot) p++;
        else break;
    }
    for (int i = p + 1; i <= last; i++) {
        if (array[i] < pivot) {
            std::swap(array[i], array[p]);
            p++;
        }
    }
    std::swap(array[first], array[p - 1]);
    return p - 1;
}

void sequentialQuickSort(int* array, int first, int last) {
    if (first < last) {
        int j = standardPartition(array, first, last);
        sequentialQuickSort(array, first, j - 1);
        sequentialQuickSort(array, j + 1, last);
    }
}

/// Optimizirana particija, tako da se za pivot uzima random element
int optimizedPartition(int* array, int first, int last) {
    int randomPosition = getRandomInt(first, last);   // random pivot element
    int pivot = array[randomPosition];
    int p = first + 1;
    while (p <= last) {
        if (array[p] <= pivot) p++;
        else break;
    }
    for (int i = p + 1; i <= last; i++) {
        if (array[i] < pivot) {
            std::swap(array[i], array[p]);
            p++;
        }
    }
    std::swap(array[first], array[p - 1]);
    return p - 1;
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
            int j = optimizedPartition(array, first, last);

#pragma omp single nowait
            {
#pragma omp task default(none) shared(array) private(first,j,sequentialLimit)
                quickSortTasks(array, first, j - 1, sequentialLimit);
#pragma omp task default(none) shared(array) private(last,j,sequentialLimit)
                quickSortTasks(array, j + 1, last, sequentialLimit);
            }
#pragma omp taskwait
        }
    }
}

int main() {
    //todo: neka ovakva inicijalizacija bude samo privremeno dok ne skontamo kako radi openMP
    // poslije bi trebali napraviti klasu sa različitim metodama za inicijalizaciju:
    // double, int, najgori slučaj (opadajuće sortiran), prosječan slučaj (random elementi)...
    int n = 1000000;
    int* array1 = new int[n];
    int* array2 = new int[n];
    //    for (int i = 0; i < n; i++) {
    //        array1[i] = n - i;
    //        array2[i] = n - i;
    //    }
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
    std::cout << "Vrijeme izvrsenja: " << executionTime << " ms." << std::endl;
    writeResultInFile("rezultati.csv", 1, executionTime);

    // paralelizacija sa taskovima
    timeBefore = std::clock();
    quickSortTasks(array2, 0, n - 1, 1000);
    timeAfter = std::clock();

    // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
    if (isSorted(array2, n)) {
        std::cout << "OK" << std::endl;
    }
    else {
        std::cout << "Nije OK" << std::endl;
    }

    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::cout << "Vrijeme izvrsenja: " << executionTime << " ms." << std::endl;
    writeResultInFile("rezultati.csv", 4, executionTime);

    // način da se dobije broj hardverski podržanih niti na računaru uz pomoć open mp
    //int numThreads = omp_get_num_procs();
    // std::cout << "Broj niti: " << numThreads << std::endl;


    delete[] array1;
    delete[] array2;
    return 0;
}

// za n=10000000 preko 2 min i 20 sek