# Matrix-Generator
The **genMat** program generates sparse matrices with several given features.

To compile the program, use the command
```
./compile.sh
```

USAGE: 

```
./genten sizes[] [options]                                                                                                                               
	-d density : nonzero ratio                              
	-s is_symmetric : whether the matrix is symmetric (1:yes)
	-c is_column : whether the degree properties are for column (1:yes)
	-m min : minimum degree value
	-i imbalance : (max-avg) / avg
	-v cv : coefficient of variation for degrees    
	-u up_bandwidth : upper matrix bandwidth (valid if matrix is square)
	-l low_bandwidth : lower matrix bandwidth (valid if matrix is square)
	-r random_seed : seed for pseudo-randomness                                                                                                                          
	-o outfile : to print out the generated matrix
	-h print_header : to print the header names for the output values 
	-b print_debug : to print at some main steps for debugging
	-w write_matrix : to write the generated matrix into a file
```

An example run command will be like the following:
```
./genmat 100 100 -d 0.01 -c 0.5 -m 1 -i 1.5 -o data/generated_100.mtx
```
