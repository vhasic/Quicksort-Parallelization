// PRS projekat.cpp : Defines the entry point for the application.
//

#include "PRS projekat.h"

#include <algorithm>
#include <random>
#include <ctime>
#include <iostream>
#include <fstream>
#include <omp.h>


int dajRandomInt(int min = 1, int max = 1000) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> randomInt(min, max); // distribution in range [1, 1000]

    return randomInt(rng);
}

void zapisiRezultatUDatoteku(std::string filename, int brojNiti, int trajanje) {
    std::ofstream izlazni_tok(filename, std::ios_base::app);
    izlazni_tok << "Broj niti: " << brojNiti << ", vrijeme izvršavanja: " << trajanje << std::endl;
}

int particijaStandardna(int* niz, int prvi, int zadnji) {
    int pivot = niz[prvi];
    int p = prvi + 1;
    while (p <= zadnji) {
        if (niz[p] <= pivot) p++;
        else break;
    }
    for (int i = p + 1; i <= zadnji; i++) {
        if (niz[i] < pivot) {
            std::swap(niz[i], niz[p]);
            p++;
        }
    }
    std::swap(niz[prvi], niz[p - 1]);
    return p - 1;
}

void quickSortSekvencijalni(int* niz, int prvi, int zadnji) {
    if (prvi < zadnji) {
        int j = particijaStandardna(niz, prvi, zadnji);
        quickSortSekvencijalni(niz, prvi, j - 1);
        quickSortSekvencijalni(niz, j + 1, zadnji);
    }
}

int particijaOptimizirana(int* niz, int prvi, int zadnji) {
    int broj = dajRandomInt(prvi, zadnji);   // random pivot element
    int pivot = niz[broj];
    int p = prvi + 1;
    while (p <= zadnji) {
        if (niz[p] <= pivot) p++;
        else break;
    }
    for (int i = p + 1; i <= zadnji; i++) {
        if (niz[i] < pivot) {
            std::swap(niz[i], niz[p]);
            p++;
        }
    }
    std::swap(niz[prvi], niz[p - 1]);
    return p - 1;
}

void quickSortTasks(int* niz, int prvi, int zadnji, int sekvencijalnaGranica) {
    if (prvi < zadnji) {
        if (zadnji - prvi < sekvencijalnaGranica)
        {
            // small list => sequential.
            return quickSortSekvencijalni(niz, prvi, zadnji);
        }
        else {

            int j = particijaOptimizirana(niz, prvi, zadnji);
            // create two tasks
#pragma omp task shared(niz)
            quickSortTasks(niz, prvi, j - 1, sekvencijalnaGranica);
#pragma omp task shared(niz)
            quickSortTasks(niz, j + 1, zadnji, sekvencijalnaGranica);
        }
    }
}

void paralellQuickSortTasks(int* niz, int prvi, int zadnji) {
    constexpr int TASK_LIMIT = 1000;
#pragma omp parallel omp_thread_limit(4)
    {
#pragma omp single
        quickSortTasks(niz, prvi, zadnji, TASK_LIMIT - 1);
#pragma omp taskwait
    }
}

// ovdje je problem jer svaka nit izvršava zasebno paralelni dio
//void quickSortParalelniNNiti(int *niz, int prvi, int zadnji) {
//    if (prvi < zadnji) {
//        int j = particijaOptimizirana(niz, prvi, zadnji);
//#pragma omp parallel sections
//        {
//#pragma omp section
//            {
//                quickSortParalelniNNiti(niz, prvi, j - 1);
//            }
//#pragma omp section
//            {
//                quickSortParalelniNNiti(niz, j + 1, zadnji);
//            }
//        }
//    }
//}


int main() {
    int n = 1000000;
    int* niz1 = new int[n];
    int* niz2 = new int[n];
    for (int i = 0; i < n; i++) {
        int broj = dajRandomInt();
        niz1[i] = broj;
        niz2[i] = broj;
    }
    clock_t vrijeme1 = std::clock();
    quickSortSekvencijalni(niz1, 0, n - 1);
    clock_t vrijeme2 = std::clock();

    int ukvrijeme = (vrijeme2 - vrijeme1) / (CLOCKS_PER_SEC / 1000);
    std::cout << "Vrijeme izvrsenja: " << ukvrijeme << " ms." << std::endl;
    zapisiRezultatUDatoteku("rezultati.txt", 1, ukvrijeme);


    //sa n niti i random pivotom
/*    vrijeme1 = std::clock();
    quickSortParalelniNNiti(niz2, 0, n - 1);
    vrijeme2 = std::clock();

    ukvrijeme = (vrijeme2 - vrijeme1) / (CLOCKS_PER_SEC / 1000);
    std::cout << "Vrijeme izvrsenja: " << ukvrijeme << " ms." << std::endl;
    zapisiRezultatUDatoteku("rezultati.txt", n, ukvrijeme);*/


    //sa taskovima
    vrijeme1 = std::clock();
    paralellQuickSortTasks(niz2, 0, n - 1);
    vrijeme2 = std::clock();

    ukvrijeme = (vrijeme2 - vrijeme1) / (CLOCKS_PER_SEC / 1000);
    std::cout << "Vrijeme izvrsenja: " << ukvrijeme << " ms." << std::endl;
    zapisiRezultatUDatoteku("rezultati.txt", n, ukvrijeme);

    int numThreads = omp_get_num_procs();
    std::cout << "Broj niti: " << numThreads << std::endl;


    delete[] niz1;
    delete[] niz2;
    return 0;
}

// za n=10000000 preko 2 min i 20 sek