#include "helperFunctions.cpp"

/**
 *
 * Funkcija izvršava verziju Quick Sorta zadanu sa parametrom versionString
 * @tparam Tip Tip elemenata niza koji se sortira
 * @param array Niz koji se sortira
 * @param n Veličina array
 * @param sequentialLimit Kada će se izvršavati sekvencijalno
 * @param numThreads Broj niti za paralelno izvršavanje
 * @param versionString Verzija koja se izvršava
 * @param worstCase true ako se treba izvršiti najgori slučaj, false ako se želi izvršiti prosječan slučaj
 * @param initializedArray Niz kojim se inicijalizira u slučaju da je worstCase = false
 */
template<typename Tip>
void executeVersionQuickSort(Tip* array, int n, int sequentialLimit, int numThreads, std::string versionString, bool worstCase, int* initializedArray){
    clock_t timeBefore, timeAfter;
    int executionTime;
    if(worstCase) {
        array = initializeArrayWorstCase(array, n);
    }
    else{
        std::copy(initializedArray, initializedArray + n, array);
    }
    timeBefore = std::clock();

    if (versionString=="Sekvencijalno-standardnaParticija-averageCase: " || versionString=="Sekvencijalno-standardnaParticija-worstCase: "){
        sequentialQuickSort(array, 0, n - 1);
//        numThreads=1; //ako se greškom pošalje neki drugi broj
    }
    else if (versionString=="Sekvencijalno-randomPivot-averageCase: " || versionString=="Sekvencijalno-randomPivot-worstCase: "){
        sequentialQuickSortRandomPivot(array, 0, n - 1);
    }
    else if (versionString=="Sekvencijalno-medianOfThreePivot-averageCase: " || versionString=="Sekvencijalno-medianOfThreePivot-worstCase: "){
        sequentialQuickSortMedianOfThreePivot(array, 0, n - 1);
    }
    else if (versionString=="Sekvencijalno-AVX2-averageCase: " || versionString=="Sekvencijalno-AVX2-worstCase: "){
        quicksort_32(array, 0, n - 1, sequentialLimit);
    }
    else if (versionString=="Sekvencijalno-bibliotecna-averageCase: " || versionString=="Sekvencijalno-bibliotecna-worstCase: "){
        auto compare = [](const void *a, const void *b) {
            return (*(int *) a - *(int *) b);
        };
        std::qsort(array, n, sizeof(Tip), compare);
    }
    else if (versionString=="Paralelno-Taskovi-medianOfThreePivot-averageCase: " || versionString=="Paralelno-Taskovi-medianOfThreePivot-worstCase: "){
#pragma omp parallel default(none) firstprivate(sequentialLimit) shared(array, n) num_threads(numThreads)
        {
#pragma omp single nowait
            {
                quickSortTasks(array, 0, n - 1, sequentialLimit);
            }
        }
    }
    else if (versionString=="Paralelno-Taskovi-medianOfThreePivot-AVX-averageCase: " || versionString=="Paralelno-Taskovi-medianOfThreePivot-AVX-worstCase: "){
#pragma omp parallel default(none) firstprivate(sequentialLimit) shared(array, n) num_threads(numThreads)
        {
#pragma omp single nowait
            {
                quickSortTasksAVX(array, 0, n - 1, sequentialLimit);
            }
        }
    }
    else if (versionString=="Paralelno-Sekcije-medianOfThreePivot-averageCase: " || versionString=="Paralelno-Sekcije-medianOfThreePivot-worstCase: "){
        quickSortSections(array, 0, n - 1, sequentialLimit, numThreads);
    }

    timeAfter = std::clock();
    executionTime = (timeAfter - timeBefore) / (CLOCKS_PER_SEC / 1000);
    std::cout<<"Vrijeme izvrsenja za n= "<<n<<" i "<< versionString<< executionTime<<std::endl;
    // provjera da li je niz dobro sortiran, tj. da li je paralelizacija korektna
    if (isSorted(array, n)) {
        std::printf("OK\n");
        writeResultInFile(numThreads, executionTime, n, versionString);
    } else {
        std::printf("Nije OK\n");
    }
}