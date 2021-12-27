#include <algorithm>
#include <random>
#include <functional>

/// Vraća random generisani integer između min i max
int getRandomInt(int min = 1, int max = 1000) {
    // generate thread safe random number
    std::uniform_int_distribution<int> dice_distribution(min, max);
    std::mt19937 random_number_engine; // pseudorandom number generator
    auto dice_roller = std::bind(dice_distribution, random_number_engine);
    int random_roll = dice_roller();  // Generate one of the integers [min,max]
    return random_roll;
}

///Lomuto's partition Scheme
int standardPartition(int* array, int low, int high) {
    int pivot = array[high];
    // Index of smaller element
    int i = low - 1;

    for(int j = low; j <= high-1; j++){
        // If current element is smaller than or equal to pivot
        if(array[j] <= pivot){
            // increment index of smaller element
            i++;
            std::swap(array[i], array[j]);
        }
    }
    std::swap(array[i+1], array[high]);
    return i+1;
}

/// Generates Random Pivot, swaps pivot with end element and calls the partition function
/// O(1.386 * n * log n)
int partition_randomPivot(int* array, int low, int high) {
/*    // generate thread safe random number
    std::uniform_int_distribution<int> dice_distribution(low, high);
    std::mt19937 random_number_engine; // pseudorandom number generator
    auto dice_roller = std::bind(dice_distribution, random_number_engine);
    int random_roll = dice_roller();  // Generate one of the integers [low,high]*/
    int random_roll = getRandomInt(low,high);

    std::swap(array[random_roll], array[high]);
    return standardPartition(array, low, high);
}

/// Optimizirana particija, tako da se za pivot uzima Median-of-three O(1.188*n*log n)
///Median-of-three: https://en.wikipedia.org/wiki/Quicksort#:~:text=Median%2Dof%2Dthree%20code%20snippet%20for%20Lomuto%20partition%3A
int partition_medianOfThreePivot(int* array, int low, int high) {
    int mid = (low + high) / 2;
    if (array[mid] < array[low]) {
        std::swap(array[low], array[mid]);
    }
    if (array[high] < array[low]) {
        std::swap(array[low], array[high]);
    }
    if (array[mid] < array[high]) {
        std::swap(array[mid], array[high]);
    }
    return standardPartition(array, low, high);
}
