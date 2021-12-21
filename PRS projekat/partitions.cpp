#include <algorithm>
#include <random>

/// Vraća random generisani integer između min i max
int getRandomInt(int min = 1, int max = 1000) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> randomInt(min, max);

    return randomInt(rng);
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

/// Optimizirana particija, tako da se za pivot uzima random element O(1.386 * n * log n)
int partition_randomPivot(int* array, int first, int last) {
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

/// Optimizirana particija, tako da se za pivot uzima Median-of-three O(1.188*n*log n)
///Median-of-three: https://en.wikipedia.org/wiki/Quicksort#:~:text=Median%2Dof%2Dthree%20code%20snippet%20for%20Lomuto%20partition%3A
int partition_medianOfThreePivot(int* array, int first, int last) {
    int mid = (first + last) / 2;
    if (array[mid] < array[first]) {
        std::swap(array[first], array[mid]);
    }
    if (array[last] < array[first]) {
        std::swap(array[first], array[last]);
    }
    if (array[mid] < array[last]) {
        std::swap(array[mid], array[last]);
    }
    int pivot = array[last];

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
