#ifndef QUICKSORT_H_
#define QUICKSORT_H_

#include <mpi.h>

/**
 * Verify that elements are sorted in ascending order. If not, write an error
 * message to stdout. Thereafter, print all elements (in the order they are
 * stored) to a file named according to the last argument.
 * @param elements Elements to check and print
 * @param n Number of elements
 * @param file_name Name of output file
 * @return 0 on success, -2 on I/O error
 */
int check_and_print(int *elements, int n, char *file_name);

/**
 * Distribute all elements from root to the other processes as evenly as
 * possible. Note that this method allocates memory for my_elements. This must
 * be freed by the caller!
 * @param all_elements Elements to distribute (Not significant in other processes)
 * @param n Number of elements in all_elements
 * @param my_elements Pointer to buffer where the local elements will be stored
 * @return Number of elements received by the current process
 */
int distribute_from_root(int *all_elements, int n, int **my_elements);

/**
 * Gather elements from all processes on root. Put root's elements first and
 * thereafter elements from the other nodes in the order of their ranks (so that
 * elements from process i come after the elements from process i-1).
 * @param all_elements Buffer on root where the elements will be stored
 * @param my_elements Elements to be gathered from the current process
 * @param local_n Number of elements in my_elements
 */
void gather_on_root(int *all_elements, int *my_elements, int local_n);

/**
 * Perform the global part of parallel quick sort. This function assumes that
 * the elements is sorted within each node. When the function returns, all
 * elements owned by process i are smaller than or equal to all elements owned
 * by process i+1, and the elements are sorted within each node.
 * @param elements Pointer to the array of sorted values on the current node. Will point to a(n) (new) array with the sorted elements when the function returns.
 * @param n Length of *elements
 * @param MPI_Comm Communicator containing all processes participating in the global sort
 * @param pivot_strategy Tells how to select the pivot element. See documentation of select_pivot in pivot.h.
 * @return New length of *elements
 */
int global_sort(int **elements, int n, MPI_Comm, int pivot_strategy);

/**
 * Merge v1 and v2 to one array, sorted in ascending order, and store the result
 * in result.
 * @param v1 Array to merge
 * @param n1 Length of v1
 * @param v2 Array to merge
 * @param n2 Length of v2
 * @param result Array for merged result (must be allocated before!)
 */
void merge_ascending(int *v1, int n1, int *v2, int n2, int *result);

/**
 * Read problem size and elements from the file whose name is given as an
 * argument to the function and populate elements accordingly. Note that this
 * method allocates memory for elements. This must be freed by the caller!
 * @param file_name Name of input file
 * @param elements Pointer to array to be created and populated
 * @return Number of elements read and stored
 */
int read_input(char *file_name, int **elements);

/**
 * Check if a number of elements are sorted in ascending order. If they aren't,
 * print an error message specifying the first two elements that are in wrong
 * order.
 * @param elements Array to check
 * @param n Length of elements
 * @return 1 if elements is sorted in ascending order, 0 otherwise
 */
int sorted_ascending(int *elements, int n);

/**
 * Swap the values pointed at by e1 and e2.
 */
void swap(int *e1, int *e2);

#endif
