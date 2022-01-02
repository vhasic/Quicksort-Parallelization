#include "quickSortVersions.cpp"
#include <ctime>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <vector>


/// Vraća true ako datoteka sa nazivom name postoji, inače vraća false
bool fileExists (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    }
    return false;
}

/// Upis razultata u datoteku filename u csv formatu broj niti, vrijeme izvršavanja
void writeResultInFile(int brojNiti, int trajanje, int velicinaNiza, std::string nacin, std::string filename = "rezultati.csv") {
    bool datotekaNePostoji = false;
    if(!fileExists(filename)){
        datotekaNePostoji = true;
    }
    std::ofstream izlazni_tok(filename, std::ios_base::app);
    // upisivanje prvog reda (nazivi kkolona), ako je datoteka tek kreirana, kako bi se poslije moglo pristupati po kolonama kroz csv
    if(datotekaNePostoji){
        izlazni_tok << "nacin" << "," << "velicinaNiza" << "," << "brojNiti" << "," << "trajanje" << std::endl;
    }
//    izlazni_tok << brojNiti << "," << velicinaNiza << "," << trajanje << "ms," << nacin << std::endl;
    izlazni_tok << nacin << "," << velicinaNiza << "," << brojNiti << "," << trajanje << std::endl;
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

bool isSorted32(const uint32_t *a, int n) {
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

void showArray32(uint32_t * array, int first, int last){
    for (int i = 0; i<last-first+1;i++) {
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

uint32_t* initializeArrayWorstCase32(uint32_t *array, int n){
    for (int i = 0; i < n; i++) {
//        array[i] = i+1;
        array[i] = n-i;
    }
    return array;
}

// ovakav način inicijalizacije ima problem, jer se svakog puta inicijalizira drugačije niz
// a zbog izbora pivota to ima velikog utjecaja na vrijeme izvršavanja
// pa se različita vremena izvršavanja ne mogu porediti;
// jedino da se random inicijalizira jedan globalni niz i da se on uvijek kopira u druge nizove
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
    int n = 200000;
    int *array = new int[n];
    const int numHarwareThreads = omp_get_num_procs();     // način da se dobije broj hardverski podržanih niti na računaru uz pomoć open mp
    const int sequentialLimit = n/(numHarwareThreads*4);

//    array= initializeArrayWorstCase(array,n);
//    array= initializeArrayAverageCase(),n);

    //sekvencijalna verzija sa standardnom particijom
//    array= initializeArrayWorstCase(array,n);
//    timeBefore = std::clock();
//    sequentialQuickSort(array, 0, n - 1);
//    timeAfter = std::clock();
//    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
//    std::printf("Vrijeme izvrsenja: %d ms\n", executionTime);
//    if (isSorted(array, n)) {
//        std::printf("OK\n");
//        writeResultInFile(1, executionTime,n,"Sekvencijalno-standardnaParticija");
//    }
//    else {
//        std::printf("Nije OK\n");
//    }

    //sekvencijalna verzija sa random pivotom u particiji
//    array= initializeArrayWorstCase(array,n);
//    timeBefore = std::clock();
//    sequentialQuickSortRandomPivot(array, 0, n - 1);
//    timeAfter = std::clock();
//    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
//    std::printf("Vrijeme izvrsenja: %d ms\n", executionTime);
//    // provjera da li je niz dobro sortiran
//    if (isSorted(array, n)) {
//        std::printf("OK\n");
//        writeResultInFile(1, executionTime,n,"Sekvencijalno-randomPivot");
//    }
//    else {
//        std::printf("Nije OK\n");
//    }

    //sekvencijalna verzija sa median of three pivotom u particiji
//    array= initializeArrayWorstCase(array,n);
//    timeBefore = std::clock();
//    sequentialQuickSortMedianOfThreePivot(array, 0, n - 1);
//    timeAfter = std::clock();
//    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
//    std::printf("Vrijeme izvrsenja: %d ms\n", executionTime);
//    // provjera da li je niz dobro sortiran
//    if (isSorted(array, n)) {
//        std::printf("OK\n");
//        writeResultInFile(1, executionTime,n,"Sekvencijalno-medianOfThreePivot");
//    }
//    else {
//        std::printf("Nije OK\n");
//    }


    // paralelizacija sa sekcijama
//    array= initializeArrayWorstCase(array,n);
//    timeBefore = std::clock();
//    quickSortSections(array, 0, n - 1,sequentialLimit);
//    timeAfter = std::clock();
//    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
//    std::printf("Vrijeme izvrsenja: %d ms\n", executionTime);
//    writeResultInFile(4, executionTime,n,"Paralelno-Sekcije-medianOfThreePivot");
//    if (isSorted(array, n)) {
//        std::printf("OK\n");
//    }
//    else {
//        std::printf("Nije OK\n");
//    }


    // paralelizacija sa taskovima
    std::vector<int> arrayNumberOfThreads({4});
    for (int numThreads: arrayNumberOfThreads) {
        array= initializeArrayWorstCase(array,n);
        timeBefore = std::clock();
#pragma omp parallel default(none) firstprivate(sequentialLimit) shared(array,n) num_threads(numThreads)
        {
#pragma omp single nowait
            {
                quickSortTasks(array, 0, n - 1, sequentialLimit);
            }
        }
        timeAfter = std::clock();

        executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
        std::printf("Vrijeme izvrsenja: %d ms\n", executionTime);
        writeResultInFile(numThreads, executionTime,n,"Paralelno-Taskovi-medianOfThreePivot");

        // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
        if (isSorted(array, n)) {
            std::printf("OK\n");
        }
        else {
            std::printf("Nije OK\n");
        }

//        verzija 2

        array= initializeArrayWorstCase(array,n);
        timeBefore = std::clock();
        quickSortTasks_v2(array,0,n-1,sequentialLimit);
        timeAfter = std::clock();

        executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
        std::printf("Vrijeme izvrsenja: %d ms\n", executionTime);
        writeResultInFile(numThreads, executionTime,n, "Paralelno-Taskovi_v2-medianOfThreePivot");

        // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
        if (isSorted(array, n)) {
            std::printf("OK\n");
        }
        else {
            std::printf("Nije OK\n");
        }


        //        new shit 2
        int n2= 75;
        uint32_t *array2 = new uint32_t[n2];
        array2= initializeArrayWorstCase32(array2, n2);
        int i = 0;
        int j = n2-1;
        int j1 = j+1;

        const uint32_t pivot = array2[(i + j)/2];

        timeBefore = std::clock();
        quicksort_32(array2,i, j);
        timeAfter = std::clock();

        executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
        std::printf("Vrijeme izvrsenja asim: %d ms\n", executionTime);
        writeResultInFile(numThreads, executionTime,n2, "AVX2");

        // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
        if (isSorted32(array2, n2)) {
            std::printf("OK\n");
            showArray32(array2, 0, n2-1);
        }
        else {
            std::printf("Nije OK\n");

        }
        delete[] array2;

    }

    delete[] array;
    return 0;
}