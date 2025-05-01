#!/bin/bash

#SBATCH -A uppmax2025-2-247
#SBATCH -p core -n 4
#SBATCH -t 00:10:00
#SBATCH -J A3-check

################################################################################
# Unpack the file A3.tar.gz (which should be on the format specified in
# the instructions of Assignment 3 in Parallel and Distributed
# Programming spring 2024). Check that the required files exist and can be built
# respectively. Run the binary with different input on different number of
# processes. Exit with exit code 1 if the tar file structure is not the required
# one or the program computes the wrong result. Then also print a message saying
# what's wrong. If everything seems OK, print a message saying that the file is
# ready for submission, and exit with exit code 0.
################################################################################

# Extract tar, enter directory, build (clean first and check that files were removed)
echo "Extracting and entering A3 directory"
tar -xzf A3.tar.gz || exit 1
cd A3 || exit 1
report=`find . -maxdepth 1 -name "A3_Report.pdf" | wc -l`
if [ 1 != "$report" ]; then
    echo "A3_Report.pdf is missing!"
    exit 1
fi
echo "Building"
make clean || exit 1
binary=`find . -maxdepth 1 -name "quicksort" | wc -l`
if [ 0 != "$binary" ]; then
    echo "make clean does not work as expected!"
    exit 1
fi
make || exit 1
binary=`find . -maxdepth 1 -name "quicksort" | wc -l`
if [ 1 != "$binary" ]; then
    echo "quicksort not built by 'make all'!"
    exit 1
fi
echo "OK"

pwd

# Check the output of the program
data_dir=/proj/uppmax2025-2-247/nobackup/A3/inputs
output_file="output_to_remove.txt"
commands=(
	"./quicksort ${data_dir}/input10.txt ${output_file} 1"
	"./quicksort ${data_dir}/input10.txt ${output_file} 2"
	"./quicksort ${data_dir}/input10.txt ${output_file} 3"
	"./quicksort ${data_dir}/input3.txt ${output_file} 3"
)
expected_output=(
	"^\s*35005211\s+294702567\s+304089172\s+521595368\s+1303455736\s+1365180540\s+1540383426\s+1726956429\s+1967513926\s+2044897763\s*$"
	"^\s*35005211\s+294702567\s+304089172\s+521595368\s+1303455736\s+1365180540\s+1540383426\s+1726956429\s+1967513926\s+2044897763\s*$"
	"^\s*35005211\s+294702567\s+304089172\s+521595368\s+1303455736\s+1365180540\s+1540383426\s+1726956429\s+1967513926\s+2044897763\s*$"
	"^\s*846930886\s+1681692777\s+1804289383\s*$"
)

echo "Checking result of serial runs"
for (( i = 0; i < ${#commands[@]}; i++ )); do
	cmd=${commands[$i]}
	$cmd > /dev/null
	cat ${output_file} | tr '\n' ' ' > tmp_check_output_file
	mv tmp_check_output_file ${output_file}
	egrep -q ${expected_output[$i]} ${output_file}
    if [ 0 != $? ] ; then
        echo "Wrong result for command $cmd"
        exit 1
    fi  
    rm ${output_file}
done
echo "OK"

echo "Checking format of output"
output=`${commands[0]}`
if ! [[ $output =~ ^[[:space:]]*[[:digit:]]*\.[[:digit:]]+[[:space:]]*$ ]]; then
     >&2 echo "Wrong format on output: \"$output\" should be one floating point number."
     exit 1
fi
rm ${output_file}
echo "OK"

echo "Checking result of parallel runs"
p=4
for (( i = 0; i < ${#commands[@]}; i++ )); do
	cmd=${commands[$i]}
	output_lines=`mpirun --bind-to none -np $p $cmd | wc -l`
	if [ 1 -lt ${output_lines} ]; then
		echo "Too many output lines! Is your program parallelized?"
		exit 1
	fi
	cat ${output_file} | tr '\n' ' ' > tmp_check_output_file
	mv tmp_check_output_file ${output_file}
	egrep -q ${expected_output[$i]} ${output_file}
    if [ 0 != $? ] ; then
        echo "Wrong result for command $cmd ($p processes)"
        exit 1
    fi  
    rm ${output_file}
done
echo "OK"

make clean || exit 1
echo "Your file is ready for submission. Well done!"

