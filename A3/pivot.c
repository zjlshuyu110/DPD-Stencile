#include "pivot.h"
#include <stdlib.h>

// Compares two integers for qsort
int compare(const void *v1, const void *v2) {
    int a = *(const int *)v1;
    int b = *(const int *)v2;
    return a - b;
}

// Returns the index of the first element larger than val, or n if none are found
int get_larger_index(int *elements, int n, int val) {
    for (int i = 0; i < n; ++i) {
        if (elements[i] > val)
            return i;
    }
    return n;
}

// Returns the median value from the sorted elements
int get_median(int *elements, int n) {
    if (n == 0) return 0; // or handle error
    if (n % 2 == 1)
        return elements[n / 2];
    else
        return (elements[n / 2 - 1] + elements[n / 2]) / 2;
}

// Selects a pivot based on the strategy and elements provided
int select_pivot(int pivot_strategy, int *elements, int n, MPI_Comm communicator) {
    switch (pivot_strategy) {
        case MEDIAN_ROOT:
            return select_pivot_median_root(elements, n, communicator);
        case MEAN_MEDIAN:
            return select_pivot_mean_median(elements, n, communicator);
        case MEDIAN_MEDIAN:
            return select_pivot_median_median(elements, n, communicator);
        default:
            return select_pivot_smallest_root(elements, n, communicator);
    }
}

// Selects the pivot using the median-of-medians strategy
int select_pivot_median_root(int *elements, int n, MPI_Comm communicator) {
    int median = get_median(elements, n);
    return get_larger_index(elements, n, median);
}

// Selects the pivot using the mean-of-medians strategy (placeholder)
int select_pivot_mean_median(int *elements, int n, MPI_Comm communicator) {
    // Placeholder: implement mean of medians logic as needed
    int median = get_median(elements, n);
    return get_larger_index(elements, n, median);
}

// Selects the pivot using the median-of-medians strategy (placeholder)
int select_pivot_median_median(int *elements, int n, MPI_Comm communicator) {
    // Placeholder: implement median of medians logic as needed
    int median = get_median(elements, n);
    return get_larger_index(elements, n, median);
}

// Selects the pivot as the smallest element's index (fallback strategy)
int select_pivot_smallest_root(int *elements, int n, MPI_Comm communicator) {
    if (n == 0) return 0;
    int smallest = elements[0];
    return get_larger_index(elements, n, smallest);
}