#include "executionVersions.cpp"


int main() {

    std::vector<int> arrayNumberOfElements({100'000, 500'000, 1'000'000, 5'000'000, 10'000'000, 50'000'000, 100'000'000});
//    std::vector<int> arrayNumberOfElements({100'000});
    int br=0;
    for (int n: arrayNumberOfElements) {
        int *initializedArray = new int[n];
        initializedArray = initializeArrayAverageCase(initializedArray, n);
//    showArray(initializedArray,0,n-1);

        int *array = new int[n];
        uint32_t *array2 = new uint32_t[n];

        // Neka se ove sekvencijalne verzije samo jednom izvrše da se uporedi da je ona Median of three najbolja
        // i nju dalje koristiti kako nebismo čekali satima da se izvrši
        if(br==0){
            executeVersionQuickSort(array2,n,n,1,"Sekvencijalno-standardnaParticija-worstCase: ",true,initializedArray);
            executeVersionQuickSort(array2,n,n,1,"Sekvencijalno-standardnaParticija-averageCase: ",false,initializedArray);

            executeVersionQuickSort(array2,n,n,1,"Sekvencijalno-randomPivot-worstCase: ",true,initializedArray);
            executeVersionQuickSort(array2,n,n,1,"Sekvencijalno-randomPivot-averageCase: ",false,initializedArray);

            executeVersionQuickSort(array2,n,n,1,"Sekvencijalno-medianOfThreePivot-worstCase: ",true,initializedArray);
            executeVersionQuickSort(array2,n,n,1,"Sekvencijalno-medianOfThreePivot-averageCase: ",false,initializedArray);
        }
        br++;

        executeVersionQuickSort(array2,n,n,1,"Sekvencijalno-medianOfThreePivot-worstCase: ",true,initializedArray);
        executeVersionQuickSort(array2,n,n,1,"Sekvencijalno-medianOfThreePivot-averageCase: ",false,initializedArray);

        executeVersionQuickSort(array2,n,n,1,"Sekvencijalno-AVX2-worstCase: ",true,initializedArray);
        executeVersionQuickSort(array2,n,n,1,"Sekvencijalno-AVX2-averageCase: ",false,initializedArray);

        executeVersionQuickSort(array2,n,n,1,"Sekvencijalno-bibliotecna-worstCase: ",true,initializedArray);
        executeVersionQuickSort(array2,n,n,1,"Sekvencijalno-bibliotecna-averageCase: ",false,initializedArray);


        //************************** Paralelne verzije***************************

        const int numHarwareThreads = omp_get_num_procs();     // način da se dobije broj hardverski podržanih niti na računaru uz pomoć openMp
        std::vector<int> arrayNumberOfThreads({ numHarwareThreads, numHarwareThreads*8, numHarwareThreads*16, numHarwareThreads*32});
        for (int numThreads: arrayNumberOfThreads) {
            int sequentialLimit = n / numThreads;

            executeVersionQuickSort(array2,n,sequentialLimit,numThreads,"Paralelno-Taskovi-medianOfThreePivot-worstCase: ",true,initializedArray);
            executeVersionQuickSort(array2,n,sequentialLimit,numThreads,"Paralelno-Taskovi-medianOfThreePivot-averageCase: ",false,initializedArray);

            executeVersionQuickSort(array2,n,sequentialLimit,numThreads,"Paralelno-Taskovi-AVX-worstCase: ",true,initializedArray);
            executeVersionQuickSort(array2,n,sequentialLimit,numThreads,"Paralelno-Taskovi-AVX-averageCase: ",false,initializedArray);


            executeVersionQuickSort(array2,n,sequentialLimit,numThreads,"Paralelno-Sekcije-medianOfThreePivot-worstCase: ",true,initializedArray);
            executeVersionQuickSort(array2,n,sequentialLimit,numThreads,"Paralelno-Sekcije-medianOfThreePivot-averageCase: ",false,initializedArray);
        }

        delete[] array;
        delete[] initializedArray;
        delete[] array2;
    }
    return 0;
}