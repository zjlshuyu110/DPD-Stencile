#include "stencil.h"
#include <stdio.h>
#include <stdlib.h>	
#include <mpi.h>

//Each stencil step must use only the values from the previous step — no partial updates.

//That’s why we have two arrays: input[] and output[].
// therefor we need previous vlues to calculate the next values.even the neighbor's stenciel also basied on the previous values.
// like center pointer we alredy stenciel updated, but the neighbor's stenciel still use previosr center value.
int main(int argc, char **argv) {
	// stage 1. distrubuted the data
	// Initialize MPI
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	if (4 != argc) {
		printf("Usage: stencil input_file output_file number_of_applications\n");
		return 1;
	}
	char *input_name = argv[1]; 
	char *output_name = argv[2];
	int num_steps = atoi(argv[3]); 

	// Read input file
	double *input;
	int num_values;
	// if (0 > (num_values = read_input(input_name, &input))) {
	// 	return 2;
	// }
	if (rank == 0) {
		if (0 > (num_values = read_input(input_name, &input))) {
			MPI_Abort(MPI_COMM_WORLD, 2); // abort all processes
			//MPI_Finalize(); return 2;  also can be used to finalize the MPI environment

		}
	}
	// Broadcast num values to every other processes. so they can allocate memory
	MPI_Bcast(&num_values, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// then we know the size of the input array, so we can allocate memory for it.
    // Calculate the local size for each process
	// for the unblanced unseen data, use the formualr from assignment 1:"local_N = basic + (rank < remainder ? 1 : 0);" and start_index = rank * basic + (rank < remainder ? rank : remainder);


	int basic = num_values / size;
	int remainder = num_values % size;
	
	// Calculate local size and starting index
	int local_N = basic + (rank < remainder ? 1 : 0);
	int start_index = rank * basic + (rank < remainder ? rank : remainder);
	
	// Include ghost cells in the allocation
	const int EXTENT = 2; // STENCIL_WIDTH = 5, there for the largest ghost cell is 2
	const int local_N_with_ghosts = local_N + 2 * EXTENT;


	// Allocate memory for local input and output buffers, since the allocat is expensive at the system level, so first allocat goest memory also
	double* local_input = (double*)malloc(local_N_with_ghosts * sizeof(double));
	if (local_input == NULL) {
		perror("malloc failed for local_input");
		MPI_Finalize();
		return 4;
	}

	double* local_output = (double*)malloc(local_N_with_ghosts * sizeof(double));
	if (local_output == NULL) {
		perror("malloc failed for local_output");
		free(local_input); // Free previously allocated memory
		MPI_Finalize();
		return 4;
	}
	
	// Initialize ghost cells to 0 (optional, for safety)
	for (int i = 0; i < EXTENT; i++) {
		local_input[i] = 0.0; // Left ghost cells
		local_input[local_N + EXTENT + i] = 0.0; // Right ghost cells
	}

	//MPI_Scatter assumes equal pieces for all processes, but now your pieces are different sizes.
	//we need to use MPI_Scatterv instead of MPI_Scatter. (learn from ChatGPT)
	int* sendcounts = NULL;
	int* displacements = NULL;
	if (rank == 0) {
		sendcounts = (int*)malloc(size * sizeof(int));
		displacements = (int*)malloc(size * sizeof(int));
        if (sendcounts == NULL || displacements == NULL) {
            perror("malloc failed for sendcounts or displacements");
            MPI_Abort(MPI_COMM_WORLD, 4);
        }

		int current_disp = 0;
		for (int r = 0; r < size; r++) {
			sendcounts[r] = basic + (r < remainder ? 1 : 0);
			displacements[r] = current_disp;
			current_disp += sendcounts[r];
		}
	}
	

	MPI_Scatterv(
			input, sendcounts, displacements, MPI_DOUBLE,
			local_data + EXTENT, local_N, MPI_DOUBLE,
			0, MPI_COMM_WORLD
		);

	// Free helper arrays on rank 0
    if (rank == 0) {	
        free(sendcounts);
        free(displacements);
        free(input);
    }


	// Stage 2: Ghost Cell Exchange and Boundary Communication // Exchange ghost cells  ->  Apply stencil  ->  Swap arrays  ->  Repeat
	int prev_rank = (rank - 1 + size) % size; // Left neighbor
	int next_rank = (rank + 1) % size;        // Right neighbor

	for (int s = 0; s < num_steps; s++) {
		// Exchange ghost cells
		MPI_Sendrecv(
			local_input + EXTENT, EXTENT, MPI_DOUBLE, prev_rank, 0, // Send leftmost data
			local_input + local_N + EXTENT, EXTENT, MPI_DOUBLE, next_rank, 0, // Receive right ghost cells
			MPI_COMM_WORLD, MPI_STATUS_IGNORE
		);

		MPI_Sendrecv(
			local_input + local_N, EXTENT, MPI_DOUBLE, next_rank, 1, // Send rightmost data
			local_input, EXTENT, MPI_DOUBLE, prev_rank, 1, // Receive left ghost cells
			MPI_COMM_WORLD, MPI_STATUS_IGNORE
		);

		// Apply stencil
		for (int i = EXTENT; i < local_N + EXTENT; i++) {
			double result = 0;
			for (int j = 0; j < STENCIL_WIDTH; j++) {
				result += STENCIL[j] * local_input[i - EXTENT + j];
			}
			local_output[i] = result;
		}

		// Swap input and output buffers
		if (s < num_steps - 1) {
			double *tmp = local_input;
			local_input = local_output;
			local_output = tmp;
		}
}

	// Stencil values
	double h = 2.0*PI/num_values;
	const int STENCIL_WIDTH = 5;
	const int EXTENT = STENCIL_WIDTH/2;
	const double STENCIL[] = {1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)};

	// Start timer
	double start = MPI_Wtime();

	// Allocate data for result
	// so this output is not a ghost cell
	double *output;
	if (NULL == (output = malloc(num_values * sizeof(double)))) {
		perror("Couldn't allocate memory for output");
		return 2;
	}
	// Repeatedly apply stencil--stencil is applied to the whole array, loop.
	//apply a stencil, we're calculating a new value for each point based on itself and its neighbors. but just update the value of the point itself.
	/* use : index = i - EXTENT + j : (index + num_values) % num_values which called : modulo arithmetic __ periodic boundary
	j = 0 → index = 3 - 2 + 0 = 1

	j = 1 → index = 3 - 2 + 1 = 2

	j = 2 → index = 3 - 2 + 2 = 3

	j = 3 → index = 3 - 2 + 3 = 4

	j = 4 → index = 3 - 2 + 4 = 5
	*/

	/*
	i → the current point in input[] we are updating.
	j → the position inside the stencil window.
	EXTENT → the number of neighbors we are considering on each side of the current point.

	Suppose I'm at at i = 3, and input[] = {10, 20, 30, 40, 50, 60, 70, 80}.

	for the output [3]

	Inner loop j moves like this :

	j	Meaning	Formula (index = i - EXTENT + j)	Index	input value
	0	2 steps left	3 - 2 + 0 = 1	1	20
	1	1 step left	3 - 2 + 1 = 2	2	30
	2	center	3 - 2 + 2 = 3	3	40
	3	1 step right	3 - 2 + 3 = 4	4	50
	4	2 steps right	3 - 2 + 4 = 5	5	60


	*/
	for (int s=0; s<num_steps; s++) {
		// Apply stencil-- this loop have 3 parts, the first and last part are the same, but the middle part is different.
		for (int i=0; i<EXTENT; i++) {  // this is left edge
			double result = 0;
			for (int j=0; j<STENCIL_WIDTH; j++) { // j is the offset for setencil window
				int index = (i - EXTENT + j + num_values) % num_values; // i is the current index, EXTENT is the number of neighbors to consider, j is the offset in the stencil.
				result += STENCIL[j] * input[index]; // applied stencile and sum (like weighted average)
			}
			output[i] = result; // output[i] is the new value for the current index stencil result
		}
		for (int i=EXTENT; i<num_values-EXTENT; i++) { // this is the middle part
			// this is the middle part, we just apply the stencil to the current index
			double result = 0;
			for (int j=0; j<STENCIL_WIDTH; j++) {
				int index = i - EXTENT + j;
				result += STENCIL[j] * input[index];
			}
			output[i] = result;
		}
		for (int i=num_values-EXTENT; i<num_values; i++) { // this is right edge
			double result = 0;
			for (int j=0; j<STENCIL_WIDTH; j++) {
				int index = (i - EXTENT + j) % num_values;
				result += STENCIL[j] * input[index];
			}
			output[i] = result;
		}
		// Swap input and output
		//apply the stencil to compute new values: compare to  copying data (slow), we just swap pointers.
		if (s < num_steps-1) {   //num_steps-1 is the last step, so we don't need to swap the pointers. and the numner is the stencile times

			double *tmp = input;
			input = output;
			output = tmp;
			// why there is output =tem
			//	We cannot read and write from the same array during next stencil computation.
			//No memory leak.
	Remember: we want read from input, write to output.
		}
		/*
	stencil uses periodic (looping) boundaries, the array behaves like a circle. So:
	The left edge can safely “look left” and access the end of the array
	The right edge can safely “look right” and access the start of the array
*/
	}
	free(input);
	// Stop timer
	double my_execution_time = MPI_Wtime() - start;

	// Write result
	printf("%f\n", my_execution_time);
#ifdef PRODUCE_OUTPUT_FILE
	if (0 != write_output(output_name, output, num_values)) {
		return 2;
	}
#endif

	// Clean up
	free(output);

	return 0;
}


int read_input(const char *file_name, double **values) { 
	FILE *file;
	if (NULL == (file = fopen(file_name, "r"))) {
		perror("Couldn't open input file");
		return -1;
	}
	int num_values;
	if (EOF == fscanf(file, "%d", &num_values)) {
		perror("Couldn't read element count from input file");
		return -1;
	}
	if (NULL == (*values = malloc(num_values * sizeof(double)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i=0; i<num_values; i++) {
		if (EOF == fscanf(file, "%lf", &((*values)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return num_values;
}


int write_output(char *file_name, const double *output, int num_values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_values; i++) {
		if (0 > fprintf(file, "%.4f ", output[i])) {
			perror("Couldn't write to output file");
		}
	}
	if (0 > fprintf(file, "\n")) {
		perror("Couldn't write to output file");
	}
	if (0 != fclose(file)) {
		perror("Warning: couldn't close output file");
	}
	return 0;
}
