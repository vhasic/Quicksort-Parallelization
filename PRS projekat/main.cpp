#include "quickSortVersions.cpp"
#include <ctime>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <vector>
#include <iterator>

/// Vraća true ako datoteka sa nazivom name postoji, inače vraća false
bool fileExists(const std::string &name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    }
    return false;
}

/// Upis razultata u datoteku filename u csv formatu broj niti, vrijeme izvršavanja
void writeResultInFile(int brojNiti, int trajanje, int velicinaNiza, std::string nacin,
                       std::string filename = "rezultati.csv") {
    bool datotekaNePostoji = false;
    if (!fileExists(filename)) {
        datotekaNePostoji = true;
    }
    std::ofstream izlazni_tok(filename, std::ios_base::app);
    // upisivanje prvog reda (nazivi kkolona), ako je datoteka tek kreirana, kako bi se poslije moglo pristupati po kolonama kroz csv
    if (datotekaNePostoji) {
        izlazni_tok << "nacin" << "," << "velicinaNiza" << "," << "brojNiti" << "," << "trajanje" << std::endl;
    }
//    izlazni_tok << brojNiti << "," << velicinaNiza << "," << trajanje << "ms," << nacin << std::endl;
    izlazni_tok << nacin << "," << velicinaNiza << "," << brojNiti << "," << trajanje << std::endl;
}

/// Provjerava da li je niz pravilno sortiran u rastućem poretku
template<typename Tip>
bool isSorted(const Tip *a, int n) {
    for (int i = 1; i < n; i++) {
        if (a[i - 1] > a[i]) {
            return false;
        }
    }
    return true;
}

template<typename Tip>
void showArray(Tip *array, int first, int last) {
    for (int i = 0; i <= last - first; i++) {
        std::cout << array[i] << ",";
    }
    std::cout << std::endl;
}

/// Inicijalizira niz na najgori slučaj za QuickSort algoritam, niz već sortiran
template<typename Tip>
Tip *initializeArrayWorstCase(Tip *array, int n) {
    for (int i = 0; i < n; i++) {
        array[i] = n - i;
    }
    return array;
}

/// Inicijalizira niz na prosječan slučaj, generisanjem radnom elemenata između 1 i n
template<typename Tip>
Tip *initializeArrayAverageCase(Tip *array, int n) {
    for (int i = 0; i < n; i++) {
        array[i] = getRandomInt(1, n);
    }
    return array;
}

int main() {
    clock_t timeBefore, timeAfter;
    int executionTime;
    std::string versionString="";
    int n = 1000000000;
    int *initializedArray = new int[n];
    initializedArray = initializeArrayAverageCase(initializedArray, n);
//    showArray(initializedArray,0,n-1);

    int *array = new int[n];
    uint32_t *array2 = new uint32_t[n];
    const int numHarwareThreads = omp_get_num_procs();     // način da se dobije broj hardverski podržanih niti na računaru uz pomoć openMp
    const int optimalNumberOfThreads = numHarwareThreads*4;
    const int sequentialLimit = n / optimalNumberOfThreads;

/*    //sekvencijalna verzija sa standardnom particijom
    versionString="Sekvencijalno-standardnaParticija-worstCase: ";
    array = initializeArrayWorstCase(array, n);
    timeBefore = std::clock();
    sequentialQuickSort(array, 0, n - 1);
    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
    if (isSorted(array, n)) {
        std::printf("OK\n");
//        writeResultInFile(1, executionTime, n, versionString);
    } else {
        std::printf("Nije OK\n");
    }

    //sekvencijalna verzija sa random pivotom u particiji
    versionString="Sekvencijalno-randomPivot-worstCase: ";
    array = initializeArrayWorstCase(array, n);
    timeBefore = std::clock();
    sequentialQuickSortRandomPivot(array, 0, n - 1);
    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
    // provjera da li je niz dobro sortiran
    if (isSorted(array, n)) {
        std::printf("OK\n");
//        writeResultInFile(1, executionTime, n, versionString);
    } else {
        std::printf("Nije OK\n");
    }

    //Najbolja sekvencijalna verzija sa median of three pivotom u particiji najgori slučaj
    versionString="Sekvencijalno-medianOfThreePivot-worstCase: ";
    array = initializeArrayWorstCase(array, n);
    timeBefore = std::clock();
    sequentialQuickSortMedianOfThreePivot(array, 0, n - 1);
    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
    // provjera da li je niz dobro sortiran
    if (isSorted(array, n)) {
        std::printf("OK\n");
        writeResultInFile(1, executionTime, n, versionString);
    } else {
        std::printf("Nije OK\n");
    }

    // Najbolja sekvencijalna verzija sa median of three pivotom u particiji prosječan slučaj
    versionString="Sekvencijalno-medianOfThreePivot-averageCase: ";
    std::copy(initializedArray, initializedArray + n, array);
    timeBefore = std::clock();
    sequentialQuickSortMedianOfThreePivot(array, 0, n - 1);
    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
    // provjera da li je niz dobro sortiran
    if (isSorted(array, n)) {
        std::printf("OK\n");
        writeResultInFile(1, executionTime, n, versionString);
    } else {
        std::printf("Nije OK\n");
    }

    //****************************************
    // SIMD worst case
    versionString="AVX2-worstCase: ";
    array2 = initializeArrayWorstCase(array2, n);

    timeBefore = std::clock();
    quicksort_32(array2, 0, n - 1, sequentialLimit);
    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
    // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
    if (isSorted(array2, n)) {
        std::printf("OK\n");
        writeResultInFile(1, executionTime, n, versionString);
    } else {
        std::printf("Nije OK\n");
    }

    //SIMD average case
    versionString="AVX2-averageCase: ";
    std::copy(initializedArray, initializedArray + n, array2);
    timeBefore = std::clock();
    quicksort_32(array2, 0, n - 1, sequentialLimit);
    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
    // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
    if (isSorted(array2, n)) {
        std::printf("OK\n");
        writeResultInFile(1, executionTime, n, versionString);
    } else {
        std::printf("Nije OK\n");
    }

    //*********************************************************************************
    // paralelizacija

    // paralelizacija sa sekcijama
    versionString="Paralelno-Sekcije-medianOfThreePivot-worstCase: ";
    array = initializeArrayWorstCase(array, n);
    timeBefore = std::clock();
    quickSortSections(array, 0, n - 1, sequentialLimit);
    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
    if (isSorted(array, n)) {
        std::printf("OK\n");
        writeResultInFile(4, executionTime, n, versionString);
    } else {
        std::printf("Nije OK\n");
    }

    std::vector<int> arrayNumberOfThreads({4});
    for (int numThreads: arrayNumberOfThreads) {
        // paralelizacija sa taskovima
        versionString= "Paralelno-Taskovi-medianOfThreePivot-worstCase: ";
        array = initializeArrayWorstCase(array, n);
        timeBefore = std::clock();
#pragma omp parallel default(none) firstprivate(sequentialLimit) shared(array, n) num_threads(numThreads)
        {
#pragma omp single nowait
            {
                quickSortTasks(array, 0, n - 1, sequentialLimit);
            }
        }
        timeAfter = std::clock();
        executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
        std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
        // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
        if (isSorted(array, n)) {
            std::printf("OK\n");
            writeResultInFile(numThreads, executionTime, n, versionString);
        } else {
            std::printf("Nije OK\n");
        }


        // paralelizacija sa taskovima v2
        versionString="Paralelno-Taskovi_v2-medianOfThreePivot-worstCase: ";
        array = initializeArrayWorstCase(array, n);
        timeBefore = std::clock();
        quickSortTasks_v2(array, 0, n - 1, sequentialLimit);
        timeAfter = std::clock();

        executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
        std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;

        // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
        if (isSorted(array, n)) {
            std::printf("OK\n");
            writeResultInFile(numThreads, executionTime, n, versionString);
        } else {
            std::printf("Nije OK\n");
        }

        // paralelizacija sa taskovima prosječan slučaj
        versionString="Paralelno-Taskovi_v2-medianOfThreePivot-averageCase: ";
        std::copy(initializedArray, initializedArray + n, array);
        timeBefore = std::clock();
        quickSortTasks_v2(array, 0, n - 1, sequentialLimit);
        timeAfter = std::clock();

        executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
        std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;

        // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
        if (isSorted(array, n)) {
            std::printf("OK\n");
            writeResultInFile(numThreads, executionTime, n, versionString);
        } else {
            std::printf("Nije OK\n");
        }

        // Najbolja verzija: paralelizacija sa taskovima + AVX najgori slučaj
        versionString= "Paralelno-Taskovi-medianOfThreePivot-AVX-worstCase: ";
        array2 = initializeArrayWorstCase(array2, n);
        timeBefore = std::clock();
#pragma omp parallel default(none) firstprivate(sequentialLimit) shared(array2, n) num_threads(numThreads)
        {
#pragma omp single nowait
            {
                quickSortTasksAVX(array2, 0, n - 1, sequentialLimit);
            }
        }
        timeAfter = std::clock();
        executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
        std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
        // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
        if (isSorted(array, n)) {
            std::printf("OK\n");
            writeResultInFile(numThreads, executionTime, n, versionString);
        } else {
            std::printf("Nije OK\n");
        }

        // Najbolja verzija: paralelizacija sa taskovima + AVX prosječan slučaj
        versionString= "Paralelno-Taskovi-medianOfThreePivot-AVX-averageCase: ";
        std::copy(initializedArray, initializedArray + n, array2);
        timeBefore = std::clock();
#pragma omp parallel default(none) firstprivate(sequentialLimit) shared(array2, n) num_threads(numThreads)
        {
#pragma omp single nowait
            {
                quickSortTasksAVX(array2, 0, n - 1, sequentialLimit);
            }
        }
        timeAfter = std::clock();
        executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
        std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
        // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
        if (isSorted(array, n)) {
            std::printf("OK\n");
            writeResultInFile(numThreads, executionTime, n, versionString);
        } else {
            std::printf("Nije OK\n");
        }
    }*/

    //Biliotečna funkcija najgori slučaj
    auto compare =[] (const void * a, const void * b)
    {
        return ( *(int*)a - *(int*)b );
    };
    versionString="Bibliotecna-worstCase: ";
    array2 = initializeArrayWorstCase(array2, n);
    timeBefore = std::clock();
    std::qsort(array2,n,sizeof(uint32_t), compare);
    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
    writeResultInFile(1, executionTime, n, versionString);

    // Bibliotečna funkcija prosječan slučaj
    versionString="Bibliotecna-averageCase: ";
    std::copy(initializedArray, initializedArray + n, array2);
    timeBefore = std::clock();
    std::qsort(array2,n,sizeof(uint32_t), compare);
    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
    writeResultInFile(1, executionTime, n, versionString);


    //Najbolja sekvencijalna verzija sa median of three pivotom u particiji najgori slučaj
    versionString="Sekvencijalno-medianOfThreePivot-worstCase: ";
    array2 = initializeArrayWorstCase(array2, n);
    timeBefore = std::clock();
    sequentialQuickSortMedianOfThreePivot(array2, 0, n - 1);
    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
    // provjera da li je niz dobro sortiran
    if (isSorted(array2, n)) {
        std::printf("OK\n");
        writeResultInFile(1, executionTime, n, versionString);
    } else {
        std::printf("Nije OK\n");
    }

    // Najbolja sekvencijalna verzija sa median of three pivotom u particiji prosječan slučaj
    versionString="Sekvencijalno-medianOfThreePivot-averageCase: ";
    std::copy(initializedArray, initializedArray + n, array2);
    timeBefore = std::clock();
    sequentialQuickSortMedianOfThreePivot(array2, 0, n - 1);
    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
    // provjera da li je niz dobro sortiran
    if (isSorted(array2, n)) {
        std::printf("OK\n");
        writeResultInFile(1, executionTime, n, versionString);
    } else {
        std::printf("Nije OK\n");
    }

    //Najbolja paralelna verzija
    std::vector<int> arrayNumberOfThreads({numHarwareThreads, optimalNumberOfThreads});
    for (int numThreads: arrayNumberOfThreads) {
        // Najbolja verzija: paralelizacija sa taskovima + AVX najgori slučaj
        versionString= "Paralelno-Taskovi-medianOfThreePivot-AVX-worstCase: ";
        array2 = initializeArrayWorstCase(array2, n);
        timeBefore = std::clock();
#pragma omp parallel default(none) firstprivate(sequentialLimit) shared(array2, n) num_threads(numThreads)
        {
#pragma omp single nowait
            {
                quickSortTasksAVX(array2, 0, n - 1, sequentialLimit);
            }
        }
        timeAfter = std::clock();
        executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
        std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
        // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
        if (isSorted(array2, n)) {
            std::printf("OK\n");
            writeResultInFile(numThreads, executionTime, n, versionString);
        } else {
            std::printf("Nije OK\n");
        }

        // Najbolja verzija: paralelizacija sa taskovima + AVX prosječan slučaj
        versionString= "Paralelno-Taskovi-medianOfThreePivot-AVX-averageCase: ";
        std::copy(initializedArray, initializedArray + n, array2);
        timeBefore = std::clock();
#pragma omp parallel default(none) firstprivate(sequentialLimit) shared(array2, n) num_threads(numThreads)
        {
#pragma omp single nowait
            {
                quickSortTasksAVX(array2, 0, n - 1, sequentialLimit);
            }
        }
        timeAfter = std::clock();
        executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
        std::cout<<"Vrijeme izvrsenja za "<< versionString<< executionTime<<std::endl;
        // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
        if (isSorted(array2, n)) {
            std::printf("OK\n");
            writeResultInFile(numThreads, executionTime, n, versionString);
        } else {
            std::printf("Nije OK\n");
        }
    }

    delete[] array;
    delete[] initializedArray;
    delete[] array2;
    return 0;
}