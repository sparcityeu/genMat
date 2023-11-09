# Matrix-Generator


Use the command "./compile.sh" to compile the program.

Usage: 

./genten sizes[] [options]                                                                                                                               
	-d density : nonzero ratio                                                                                                                                         
	-c cv : coefficient of variation for row (column) degrees    
	-n min : min degree value
	-i imbalance : (max-avg) / avg
	-s is_column : whether the degree properties are for column (or row)
	-b bandwidth : matrix badwidth
	-p profile : matrix profile
	-r random_seed : seed for randomness                                                                                                                          
	-o outfile : to print out the generated tensor
	-h print_header : to print the header names for the output values 
	-p print_debug : to print at some main steps for debugging
	-w write_matrix : to write the generated matrix into a file

An example run command will be like the following:

./genmat 100 100 -d 0.01 -c 0.5 -n 1 -i 1.5 -o data/generated_100.mtx