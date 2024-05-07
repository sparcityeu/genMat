# GenMat: A Sparse Matrix Generator
The **genMat** program generates sparse matrices with several given features.
It adopts a combination of normal distribution and log-normal distribution.

To compile the program, use the command
```
./compile.sh
```

USAGE: 

```
./genmat sizes[] [options]                                                                                                                               
	-d density : nonzero ratio                              
	-s is_symmetric : whether the matrix is symmetric (1:symmetric, 0:unsymmetric [default])
	-c is_column : whether the degree properties are for column (1:column, 0:row [default])
	-m min : minimum degree value [default =  1]
	-i imbalance : (max-avg) / avg (to determine max value in a size-independent way) [if not provided, max = size]
	-v cv : coefficient of variation for degrees    [default = 0.5]
	-l low_bandwidth : lower matrix bandwidth (valid if matrix is square) [if not provided, default = size-1]
	-u up_bandwidth : upper matrix bandwidth (valid if matrix is square. enter only if unsymmetric) [if not provided, default = size-1]
	-r random_seed : seed for pseudo-randomness     [default =1]                                                                                                                     
	-o outfile : to print out the generated matrix [default =1]
	-h print_header : to print the header names for the output values  [default = 0]
	-b print_debug : to print at some main steps for debugging [default = 0]
	-w write_matrix : to write the generated matrix into a file
```

Some example run commands will be like the following:
```
./genmat 1000 2000 -d 0.02 -v 0.5 -i 1.5 -o ../sample_data/generated_1000_2000.mtx
```
```
./genmat 1000 1000 -d 0.02 -v 0.5 -i 1.5 -l 200 -u 400 -o ../sample_data/generated_1000_square.mtx
```
```
./genmat 1000 1000 -d 0.02 -v 0.5 -i 1.5 -s 1 -l 200 -o ../sample_data/generated_1000_sym.mtx
``` 
