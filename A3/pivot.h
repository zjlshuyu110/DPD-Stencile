#ifndef PIVOT_H_
#define PIVOT_H_

#include <mpi.h>

// Pivot strategies
#define MEDIAN_ROOT 1
#define MEAN_MEDIAN 2
#define MEDIAN_MEDIAN 3

// Root node
#define ROOT 0

/**
 * @param v1 Pointer to first value (integer) to compare
 * @param v2 Pointer to second value (integer) to compare
 * @return 0 if *v1==*v2, a positive and negative number of *v1>*v2 and *v1<*v2 respectively
 */
int compare(const void *v1, const void *v2);

/**
 * Find the index of the first value in elements that is larger than val
 * @param elements Array to search in
 * @param n Length of elements
 * @param val Value to search for
 * @return index of first value that is larger than val, or length of array if all values are smaller
 */
int get_larger_index(int *elements, int n, int val);

/**
 * Find the median in an array. Note that this function assumes that the array
 * is sorted!
 * @param elements Sorted array
 * @param n Length of elements
 * @return median of elements
 */
int get_median(int *elements, int n);

/**
 * Select a pivot element for parallel quick sort. Return the index of the first
 * element that is larger than the pivot. Note that this function assumes that
 * elements is sorted!
 * @param pivot_strategy 0=>smallest on root (Not recommended!) 1=>median on root 2=>mean of medians 3=>median of medians
 * @param elements Elements stored by the current process (sorted!)
 * @param n Length of elements
 * @param communicator Communicator for processes in current group
 * @return The index of the first element after pivot
 */
int select_pivot(int pivot_strategy, int *elements, int n, MPI_Comm communicator);

/**
 * See select pivot!
 */
int select_pivot_median_root(int *elements, int n, MPI_Comm communicator);

/**
 * See select pivot!
 */
int select_pivot_mean_median(int *elements, int n, MPI_Comm communicator);

/**
 * See select pivot!
 */
int select_pivot_median_median(int *elements, int n, MPI_Comm communicator);

/**
 * See select pivot!
 */
int select_pivot_smallest_root(int *elements, int n, MPI_Comm communicator);

#endif
