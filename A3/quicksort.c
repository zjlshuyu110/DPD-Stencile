#include "quicksort.h"


int main(int argc, char const *argv[])
{
	return 0;
}



int check_and_print(int *elements, int n, char *file_name)
{
	if (!sorted_ascending(elements, n)) {
		printf("Error! Array failed to sort!");
	}
	FILE *file = fopen(file_name, "wb");
	fwrite(elements, sizeof(double), n, file);
	return 0;
}
