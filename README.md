# sorting
my version sorting algorithms visualization program inspired by a ton of youtube videos



# usage
call syntax:
	sort [sorting method] [number of elements] [time between frames in ms] [shape of the input]

sorting metods:
	O(n^2) < C:
		bogo
		stooge
		slow

	C = O(n^2):
		selection
		bubble	
		insertion
		cycle
		coctail
		oddEven
		gnome
		heapInsertion

	O(n*log(n)) < C < O(n^2):
		shell
		comb
	
	C = O(n*log(n)):
		partition
		merge
		heap
		partitionInsertion
		partition2

	C = O(n):
		radixLSD{base(default 256)}


input data shapes:
	U - uniformly distributed numbers
	P - perlin noise
		P [steepness]
	B - bell curve made of uniform distributions
		B [number of uniform distributions] [smoothness]
	S - sorted
	I - sorted in reversed order
	O - organ pipes shape
	V - reverse organ pipes shape
