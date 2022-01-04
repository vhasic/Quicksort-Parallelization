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
