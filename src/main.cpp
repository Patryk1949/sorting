#include <SFML/Graphics.hpp>
#include "sp_Utils.hpp"
#include <unistd.h>
#include <time.h>
#include <algorithm>
#include <random>

const uint sWidth = 1000;
const uint sHeight = 600;

ulong sleepTime = 0;

std::vector<float> arr;
sf::RenderWindow *window;


float thickness;
float highness;
sf::RectangleShape rect;

void printValues(const float *I){
	float heigth;
	window->clear();
	for (int i=0; i<(int)std::size(arr); ++i){
		heigth = arr[i]*highness;
		rect.setPosition(sf::Vector2f{(float)i*thickness, highness-heigth});
		rect.setSize(sf::Vector2f{thickness, heigth});
		window->draw(rect);
	}
	rect.setFillColor(sf::Color::Red);
	heigth = *I*highness;
	rect.setPosition(sf::Vector2f{(float)(I-&*std::begin(arr))*thickness, highness-heigth});
	rect.setSize(sf::Vector2f{thickness, heigth});
	window->draw(rect);
	rect.setFillColor(sf::Color::Yellow);
	
	window->display();
	usleep(sleepTime);
}

void printValues(const float *I, const float *J){
	float heigth;
	window->clear();
	for (int i=0; i<(int)std::size(arr); ++i){
		heigth = arr[i]*highness;
		rect.setPosition(sf::Vector2f{(float)i*thickness, highness-heigth});
		rect.setSize(sf::Vector2f{thickness, heigth});
		window->draw(rect);
	}
	rect.setFillColor(sf::Color::Blue);
	heigth = *I * highness;
	rect.setPosition(sf::Vector2f{(float)(I-&*std::begin(arr))*thickness, highness-heigth});
	rect.setSize(sf::Vector2f{thickness, heigth});
	window->draw(rect);
	rect.setFillColor(sf::Color::Red);
	heigth = *J * highness;
	rect.setPosition(sf::Vector2f{(float)(J-&*std::begin(arr))*thickness, highness-heigth});
	rect.setSize(sf::Vector2f{thickness, heigth});
	window->draw(rect);
	rect.setFillColor(sf::Color::Yellow);
	
	window->display();
	usleep(sleepTime);
}

void printValues(const float *I, const float *J, const float *K){
	float heigth;
	window->clear();
	for (int i=0; i<(int)std::size(arr); ++i){
		heigth = arr[i]*highness;
		rect.setPosition(sf::Vector2f{(float)i*thickness, highness-heigth});
		rect.setSize(sf::Vector2f{thickness, heigth});
		window->draw(rect);
	}
	rect.setFillColor(sf::Color::Blue);
	heigth = *I * highness;
	rect.setPosition(sf::Vector2f{(float)(I-&*std::begin(arr))*thickness, highness-heigth});
	rect.setSize(sf::Vector2f{thickness, heigth});
	window->draw(rect);
	rect.setFillColor(sf::Color::Green);
	heigth = *J * highness;
	rect.setPosition(sf::Vector2f{(float)(J-&*std::begin(arr))*thickness, highness-heigth});
	rect.setSize(sf::Vector2f{thickness, heigth});
	window->draw(rect);
	rect.setFillColor(sf::Color::Red);
	heigth = *K * highness;
	rect.setPosition(sf::Vector2f{(float)(K-&*std::begin(arr))*thickness, highness-heigth});
	rect.setSize(sf::Vector2f{thickness, heigth});
	window->draw(rect);
	rect.setFillColor(sf::Color::Yellow);
	
	window->display();
	usleep(sleepTime);
}

void printValues(const float *I, const float *J, const float *K, const float *L){
	float heigth;
	window->clear();
	for (int i=0; i<(int)std::size(arr); ++i){
		heigth = arr[i]*highness;
		rect.setPosition(sf::Vector2f{(float)i*thickness, highness-heigth});
		rect.setSize(sf::Vector2f{thickness, heigth});
		window->draw(rect);
	}
	rect.setFillColor(sf::Color::Blue);
	heigth = *I * highness;
	rect.setPosition(sf::Vector2f{(float)(I-&*std::begin(arr))*thickness, highness-heigth});
	rect.setSize(sf::Vector2f{thickness, heigth});
	window->draw(rect);
	heigth = *J * highness;
	rect.setPosition(sf::Vector2f{(float)(J-&*std::begin(arr))*thickness, highness-heigth});
	rect.setSize(sf::Vector2f{thickness, heigth});
	window->draw(rect);
	rect.setFillColor(sf::Color::Green);
	heigth = *K * highness;
	rect.setPosition(sf::Vector2f{(float)(K-&*std::begin(arr))*thickness, highness-heigth});
	rect.setSize(sf::Vector2f{thickness, heigth});
	window->draw(rect);
	rect.setFillColor(sf::Color::Red);
	heigth = *L * highness;
	rect.setPosition(sf::Vector2f{(float)(L-&*std::begin(arr))*thickness, highness-heigth});
	rect.setSize(sf::Vector2f{thickness, heigth});
	window->draw(rect);
	rect.setFillColor(sf::Color::Yellow);
	
	window->display();
	usleep(sleepTime);
}

void printValues(){
	float heigth;
	window->clear();
	for (int i=0; i<(int)std::size(arr); ++i){
		heigth = arr[i]*highness;
		rect.setPosition(sf::Vector2f{(float)i*thickness, highness-heigth});
		rect.setSize(sf::Vector2f{thickness, heigth});
		window->draw(rect);
	}
	
	window->display();
	usleep(sleepTime);
}


void selectionSort(float *begin, float *end){
	float *max;
	for (float *I=begin, *J; I!=end; ++I){
		for (J=I+1, max=I; J!=end; ++J){
			max = *max<*J ? max : J;
			printValues(I, max, J);
		}
		sp::swap(I, max);
	}
}

void bubbleSort(float *begin, float *end){
	for (float *I; begin!=end; --end){
		for (I=begin+1; I<end; ++I){
			if (*I < *(I-1))
				sp::swap(I, I-1);

			printValues(end, I);
		}
	}
}

void insertionSort(float *begin, float *end){
	float *J = begin;
	float temp;
	for (float *I=begin+1; I!=end; ++I){
		J = *I<*J ? I : J;
		printValues(begin, J, I);
	}
	temp = std::move(*begin);
	*begin = std::move(*J);
	*J = std::move(temp);

	for (float *I=begin+2; I!=end; ++I){
		temp = std::move(*I);
		J = I-1;
		for (; temp<*J; --J){
			*(J+1) = std::move(*J);
			printValues(I, J);
		}
		*(J+1) = std::move(temp);
	}
}

void gnomeSort(float *begin, float *end){
	for (float *I=begin+1, *J; I!=end; ++I)
		for(J=I; J!=begin; --J){
			if (*(J-1) < *J)
				break;
			sp::swap(J, J-1);

			printValues(I, J-1);
		}
}

void bogoSort(float *begin, float *end){
REPEAT:
	std::random_shuffle(begin, end);

	printValues();

	for (float *I=begin+1; I!=end; ++I)
		if (*I < *(I-1))
			goto REPEAT;
	return;
}

struct FloatView{
	float *begin, *end;
};

void partitionSort(float *begin, float *end){
	std::vector<FloatView> stak;
	stak.reserve(16);
	stak.push_back({begin, end});

	while (!stak.empty()){
		float *L = stak.back().begin;
		float *R = stak.back().end - 1;
		if (L < R){
			const float piv = *L;
			while (L < R){
				printValues(stak.back().begin, stak.back().end-1, L, R);
				while (piv<*R && L<R){
					--R;
					printValues(stak.back().begin, stak.back().end-1, L, R);
				}
				if (L < R){
					*L = *R;
					++L;
				}
				while (*L<piv && L<R){
					++L;
					printValues(stak.back().begin, stak.back().end-1, L, R);
				}
				if (L < R){
					*R = *L;
					--R;
				}
			}
			*L = piv;
			stak.push_back({L+1, stak.back().end});
			stak[std::size(stak)-2].end = L;
		} else stak.pop_back();
	}
}


float *partition(float *const first, float *last, float *const pivot){
	if (first >= last) return first;
	sp::swap(first, pivot);
	
	float lastVal = std::move(*--last);
	*last = *first;
	
	const float *const end = last;

	float *It = first;
	for (;;){
		do{ ++It;
			printValues(first, end, It, last);
		} while (*It < *first);
		*last = std::move(*It);
		do{ --last;
			printValues(first, end, It, last);
		} while (*first < *last);
		if (It >= last) break;
		*It = std::move(*last);
	}
	if (It == last+2){
		*It = std::move(*(last+1));
		--It;
	}
	{
		float *const partitionPoint = It - (*first < lastVal);
		*It = std::move(lastVal);
		sp::swap(first, partitionPoint);
		return partitionPoint;
	}
}

void partitionSort2(float *begin, float *end){
	if (begin >= end-1) return;
	
	float *const pp = partition(begin, end, begin+(end-begin)/3);

	partitionSort2(begin, pp);
	partitionSort2(pp+1, end);
}

void radixLSD(float *const first, float *const last){
	constexpr size_t baseBits = 8;
	constexpr size_t base = 1 << baseBits;
	constexpr size_t bitMask = (base-1);

	uint32_t counts[base];

	std::vector<float> buffer(last-first);

	for (uint32_t step=0; step!=sizeof(float)*8/baseBits; ++step){
		std::copy(first, last, std::begin(buffer));
		
		memset(counts, 0, sizeof(counts));
		for (float *It=&*std::begin(buffer); It!=&*std::end(buffer); ++It)
			++counts[(*(uint32_t *)It >> step*baseBits) & bitMask];

		for (uint32_t *C=std::begin(counts)+1; C!=std::end(counts); ++C) *C += *(C-1);
		
		for (float *It=&*std::end(buffer)-1; It!=&*std::begin(buffer)-1; --It){
			first[--counts[(*(uint32_t *)It >> step*baseBits) & bitMask]] = *It;
			printValues(first + counts[(*(uint32_t *)It >> step*baseBits) & bitMask]);
		}
	}
}


size_t arbaseRadixLSD_Base = 8;
void arbaseRadixLSD(float *const first, float *const last){
	const size_t base = arbaseRadixLSD_Base;
	std::vector<uint32_t> counts(base);

	std::vector<float> buffer(last-first);

	for (uint32_t divisor=1; sp::intLog2(divisor)<sizeof(float)*8; divisor*=base){
		std::fill(std::begin(counts), std::end(counts), 0);
		std::copy(first, last, std::begin(buffer));
		
		for (float *It=&*std::begin(buffer); It!=&*std::end(buffer); ++It)
			++counts[(*(uint32_t *)It / divisor) % base];

		for (uint32_t *C=&*std::begin(counts)+1; C!=&*std::end(counts); ++C) *C += *(C-1);
		
		for (float *It=&*std::end(buffer)-1; It!=&*std::begin(buffer)-1; --It){
			first[--counts[(*(uint32_t *)It / divisor) % base]] = *It;
			printValues(first + counts[(*(uint32_t *)It / divisor) % base]);
		}
	}
}

void americanFlag(float *const first, float *const last){
	constexpr size_t baseBits = 8;
	constexpr size_t base = 1 << baseBits;
	constexpr size_t bitMask = (base-1);


// FIRST ITERATION
	uint32_t counts3[base];
	uint32_t countsBackup3[base];

	memset(countsBackup3, 0, sizeof(countsBackup3));
	for (float *I=first; I!=last; ++I)
		++countsBackup3[((*(uint32_t *)I >> 3*baseBits) & bitMask)];

	for (uint32_t *C=std::begin(countsBackup3)+1; C!=std::end(countsBackup3); ++C) *C += *(C-1);

	counts3[0] = 0;
	memcpy(counts3+1, countsBackup3, sizeof(countsBackup3)-sizeof(uint32_t));

	for (float *I=first; I!=last;){
		const uint32_t index = (*(uint32_t *)I >> 3*baseBits) & bitMask;
		float *const P = first + counts3[index];

		if (I == P){
			++counts3[index];
			++I;
			continue;
		}
		if (counts3[index] == countsBackup3[index]){
			++I;
			continue;
		}

		printValues(I, P);
		std::iter_swap(I, P);
		++counts3[index];
	}
	uint32_t *c3 = counts3;
	for (float *I=first; I<last;){
		const uint32_t val3 = *c3;

// SECOND ITERATION
		uint32_t counts2[base];
		uint32_t countsBackup2[base];

		memset(countsBackup2, 0, sizeof(countsBackup2));
		for (float *J=I; J!=I+val3; ++J)
			++countsBackup2[((*(uint32_t *)J >> 2*baseBits) & bitMask)];

		for (uint32_t *C=std::begin(countsBackup2)+1; C!=std::end(countsBackup2); ++C) *C += *(C-1);

		counts2[0] = 0;
		memcpy(counts2+1, countsBackup2, sizeof(countsBackup2)-sizeof(uint32_t));
		
		for (float *J=I; J!=I+val3;){
			const uint32_t index = (*(uint32_t *)J >> 2*baseBits) & bitMask;
			float *const P = I + counts2[index];

			if (J == P){
				++counts2[index];
				++J;
				continue;
			}
			if (counts2[index] == countsBackup2[index]){
				++J;
				continue;
			}

			printValues(J, P);
			std::iter_swap(J, P);
			++counts2[index];
		}
		uint32_t *c2 = counts2;
		for (float *J=I; J<I+val3;){
			const uint32_t val2 = *c2;

// THIRD ITERATION
			uint32_t counts1[base];
			uint32_t countsBackup1[base];

			memset(countsBackup1, 0, sizeof(countsBackup1));
			for (float *K=J; K!=J+val2; ++K)
				++countsBackup2[((*(uint32_t *)K >> 1*baseBits) & bitMask)];

			for (uint32_t *C=std::begin(countsBackup1)+1; C!=std::end(countsBackup1); ++C) *C += *(C-1);

			counts1[0] = 0;
			memcpy(counts1+1, countsBackup1, sizeof(countsBackup1)-sizeof(uint32_t));
			for (float *K=J; K!=J+val2;){
				const uint32_t index = (*(uint32_t *)K >> 1*baseBits) & bitMask;
				float *const P = J + counts1[index];

				if (K == P){
					++counts1[index];
					++K;
					continue;
				}
				if (counts1[index] == countsBackup1[index]){
					++K;
					continue;
				}

				printValues(K, P);
				std::iter_swap(K, P);
				++counts1[index];
			}
			uint32_t *c1 = counts1;
			for (float *K=J; K<J+val2;){
				const uint32_t val1 = *c1;

// FOURTH ITERATION
				uint32_t counts0[base];
				uint32_t countsBackup0[base];

				memset(countsBackup0, 0, sizeof(countsBackup0));
				for (float *L=K; L!=K+val1; ++L)
					++countsBackup2[((*(uint32_t *)L >> 0*baseBits) & bitMask)];

				for (uint32_t *C=std::begin(countsBackup0)+1; C!=std::end(countsBackup0); ++C) *C += *(C-1);

				counts0[0] = 0;
				memcpy(counts0+1, countsBackup0, sizeof(countsBackup0)-sizeof(uint32_t));

				for (float *L=K; L!=K+val1;){
					const uint32_t index = (*(uint32_t *)L >> 0*baseBits) & bitMask;
					float *const P = K + counts0[index];

					if (L == P){
						++counts0[index];
						++L;
						continue;
					}
					if (counts0[index] == countsBackup0[index]){
						++L;
						continue;
					}

					printValues(L, P);
					std::iter_swap(L, P);
					++counts0[index];
				}

				K += val1;
				while (*++c1 == val1);
			}

			J += val2;
			while (*++c2 == val2);
		}

		I += val3;
		while (*++c3 == val3);
	}
}


void cycleSort(float *begin, float *end){
	for (float *I=begin; I!=end-1; ++I){
	//	float value = *I;
		float *pos;
		do{
			pos = I;
			for (float *J=I+1; J!=end; ++J){
				pos += *J < *I;

				printValues(I, pos, J);
			}
			if (pos == I)
				break;

			for (; *pos==*I; ++pos){
				printValues(I, pos);
			}
			printValues(I, pos);
			sp::swap(pos, I);
			

		} while (pos != I);

	}
}

void ccycleSort(float *begin, float *end){
	// traverse array elements and put it to on
	// the right place
	for (float *I=begin; I!=end-1; ++I){
		// Find position where we put the item. We basically
		// count all smaller elements on right side of item.
		float *pos = I;
		for (float *J=I+1; J!=end; ++J){
			pos += *J < *I;
			
			printValues(J, pos);
		}

		if (pos == I)
			continue;

		for (; *I == *pos; ++pos);
		if (pos != I)
			*pos = *I;


		while (pos != I){
			pos = I;

			// Find position where we put the element
			for (float *J=I+1; J!=end; ++J){
				pos += *J < *I;
				
				printValues(J, pos);
			}

			// ignore all duplicate  elements
			for (; *I == *pos; ++pos);

			// put the item to it's right position
			if (*pos != *I)
				sp::swap(I, pos);
		}
	}
}


namespace priv__{

void mergeRecursion(float *begin1, float *begin2, const int length, const bool swapFlag, const std::vector<float> &source){
	if (length == 1){
		*begin2 = *begin1;
		return;
	}

	const int halfLength = length>>1;
	mergeRecursion(begin1, begin2, halfLength, !swapFlag, source);
	mergeRecursion(begin1+halfLength, begin2+halfLength, length-halfLength, !swapFlag, source);

	if (swapFlag)
		sp::swap(&begin1, &begin2);

	float *I = begin1;
	float *L = begin2;
	begin1 = begin2 + length; // begin1 becomes end2
	begin2 += halfLength; // begin2 becomes mid2
	float *R = begin2;

	const float *lastDigit = I + length - 1;
	for (; L!=begin2 && R!=begin1; ++I){
		*I = *L<*R ? *(L++) : *(R++);
		printValues(I, lastDigit);
	}
	for (; L!=begin2; ++I){
		*I = *(L++);
		printValues(I, lastDigit);
	}
	for (; R!=begin1; ++I){
		*I = *(R++);
		printValues(I, lastDigit);
	}
}

}

inline void mergeSort(float *const begin, float *const end){
	const int length = end - begin;
	std::vector<float> tempArray(length);
	priv__::mergeRecursion(begin, &*std::begin(tempArray), length, false, tempArray);
}

void heapSort(float *const begin, float *const end){
	float *child;
	for (int i=(end-begin-1)>>1; i>=0; --i)
		for(int j=i; (child=begin+(j<<1)+1) < end; ){
			const bool rightChildExists = child + 1 != end;
			if (child[0] <= begin[j] && (child[1] <= begin[j] || !rightChildExists))
				break;

			child += child[1] > child[0] && rightChildExists;

			sp::swap(begin+j, child);
			j = child - begin;

			printValues(begin+i, begin+j);
		}

	for (float *I=end-1; I!=begin; --I){
		sp::swap(begin, I);
		for(int j=0; (child=begin+(j<<1)+1) < I; ){
			const bool rightChildExists = child + 1 != I;
			if (child[0] <= begin[j] && (child[1] <= begin[j] || !rightChildExists))
				break;

			child += child[1] > child[0] && rightChildExists;

			sp::swap(begin+j, child);
			j = child - begin;

			printValues(I, begin + j);
		}
	}
}

void heapInsertionSort(float *const begin, float *const end){
	float *child;
	for (int i=(end-begin-1)>>1; i>=0; --i)
		for(int j=i; (child=begin+(j<<1)+1) < end; ){
			const bool rightChildExists = child + 1 != end;
			if (child[0] >= begin[j] && (child[1] >= begin[j] || !rightChildExists))
				break;

			child += child[1] < child[0] && rightChildExists;

			sp::swap(begin+j, child);
			j = child - begin;

			printValues(begin+j, begin+i);
		}

	for (float *I=begin+2; I!=end; ++I){
		const float temp = std::move(*I);
		float *J = I-1;
		for (; temp<*J; --J){
			*(J+1) = std::move(*J);
			printValues(J, I);
		}
		*(J+1) = std::move(temp);
	}
}



void coctailShakerSort(float *begin, float *end){
	for (float *I; begin<end; ++begin, --end){
		for (I=begin+1; I<end; ++I)
			if (*I < *(I-1)){
				sp::swap(I, I-1);

				printValues(begin, end-1, I);
			}
		for (I=end-2; I>begin; --I)
			if (*I < *(I-1)){
				sp::swap(I, I-1);

				printValues(begin, end-1, I);
			}
	}
}

void oddEvenSort(float *begin, float *end){
	bool isUnsorted  = true;
	while (isUnsorted){
		isUnsorted = false;
		float *I;
		for (I=begin+1; I<end; I+=2){
			if (*I < *(I-1)){
				sp::swap(I, I-1);
				isUnsorted = true;
			}
			printValues(I);
		}
		for (float *I=begin+2; I<end; I+=2){
			if (*I < *(I-1)){
				sp::swap(I, I-1);
				isUnsorted = true;
			}
			printValues(I);
		}
	}
}


void shellSort(float *begin, float *end){
	const size_t length = end - begin;

	for (size_t gap=length>>1; gap; gap>>=1){
	
		for (float *I=begin+gap, *J, temp; I!=end; ++I){
			temp = *I;
			for(J=I-gap; J>=begin && temp<*J; J-=gap){
				*(J+gap) = *J;

				printValues(I, J, J+gap);
			}
			*(J+gap) = temp;
			printValues(I, J, J+gap);
		}
	
	}
}


void combSort(float *begin, float *end){
	for (int gap=((end-begin)*10)/13; gap; gap=(gap*10)/13)
		for (float *I=begin; I!=end-gap; ++I){
			if (*(I+gap) < *I) sp::swap(I, I+gap);
			printValues(I, I+gap);
		}
	for (float *I=begin+1; I!=end; ++I){
		if (*I < *(I-1)) sp::swap(I, I-1);
		printValues(I-1, I);
	}
}

void slowSort(float *begin, float *end){
	if (begin >= end-1)
		return;
	float *const mid = begin + ((end-begin)>>1);                            
	slowSort(begin, mid);
	slowSort(mid, end);
	
	float *const lastElPtr = end - 1;
	if (*lastElPtr < *(mid-1))
		sp::swap(lastElPtr, mid-1);
	slowSort(begin, lastElPtr);
	
	printValues(mid, end-1);
}

void stoogeSort(float *begin, float *end){
	printValues(begin, end-1);
	if (end - begin > 2){
		const int t = (end - begin) / 3;
		stoogeSort(begin  ,end-t);
		stoogeSort(begin+t, end);
		stoogeSort(begin  , end-t);
	} else if (*begin > *(end-1))
		sp::swap(begin, end-1);
}

void partitionInsertionSort(float *begin, float *end){
	if (end - begin < 100){
		shellSort(begin, end);
		return;
	}
	--end;

	float *border = begin;
	
	for (float *I=begin; I<end; ++I)
		if (*I < *end){
			sp::swap(I, border++);

			printValues(border);
		}
	
	sp::swap(border, end);

	printValues(border);

	partitionInsertionSort(begin, border);
	partitionInsertionSort(border+1, end+1);
}





int main(const int argc, const char **argv){
	if (argc == 2)
		if (!strcmp(argv[1], "help")){
			puts(
R"(
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
)"			);
			return 0;
		}

	if (argc < 3){
		puts("Too few arguments!\n");
		return 0;
	}
	char letter = 1;
	if (argc > 3){
		if (*argv[3] >= '0' && *argv[3] <= '9')
			sleepTime = 1000.f*strtof(argv[3], nullptr);
		else
			letter = *argv[3];
		if (argc > 4)
			letter = *argv[4];
	}

	arr.resize(strtol(argv[2], nullptr, 10));

	thickness = (float)sWidth/(float)std::size(arr);
	highness = (float)sHeight;

	rect.setFillColor(sf::Color::Yellow);

	void (*sortingFunction)(float *, float *);

	std::string title;

	if (!strcmp(argv[1], "selection")){
		sortingFunction = selectionSort;
		title = "selection sort";
	} else if (!strcmp(argv[1], "americanFlag")){
		sortingFunction = americanFlag;
		title = "american flag";
	} else  if (!strcmp(argv[1], "bubble")){
		sortingFunction = bubbleSort;
		title = "bubble sort";
	} else if (!strcmp(argv[1], "insertion")){
		sortingFunction = insertionSort;
		title = "insertion sort";
	} else if (!strcmp(argv[1], "bogo")){
		sortingFunction = bogoSort;
		title = "bogo sort";
	} else if (!strcmp(argv[1], "partition")){
		sortingFunction = partitionSort;
		title = "partition sort";
	} else if (!strcmp(argv[1], "partitionInsertion")){
		sortingFunction = partitionInsertionSort;
		title = "partition insertion sort";
	} else if (!strcmp(argv[1], "merge")){
		sortingFunction = mergeSort;
		title = "merge sort";
	} else if (!strcmp(argv[1], "gnome")){
		sortingFunction = gnomeSort;
		title = "gnome sort";
	} else if (!strcmp(argv[1], "heap")){
		sortingFunction = heapSort;
		title = "heap sort";
	} else if (!strcmp(argv[1], "heapInsertion")){
		sortingFunction = heapInsertionSort;
		title = "heap insertion sort";
	} else if (!strcmp(argv[1], "coctail")){
		sortingFunction = coctailShakerSort;
		title = "coctail shaker sort";
	} else if (!strcmp(argv[1], "oddEven")){
		sortingFunction = oddEvenSort;
		title = "odd even sort";
	} else if (!strcmp(argv[1], "shell")){
		sortingFunction = shellSort;
		title = "shell sort";
	} else if (!strcmp(argv[1], "comb")){
		sortingFunction = combSort;
		title = "comb sort";
	} else if (!strcmp(argv[1], "slow")){
		sortingFunction = slowSort;
		title = "slow sort";
	} else if (!strcmp(argv[1], "stooge")){
		sortingFunction = stoogeSort;
		title = "stooge sort";
	} else if (!strcmp(argv[1], "cycle")){
		sortingFunction = cycleSort;
		title = "cycle sort";
	} else if (!strcmp(argv[1], "partition2")){
		sortingFunction = partitionSort2;
		title = "other partition sort";
	} else if (!strncmp(argv[1], "radixLSD", 8)){
		sortingFunction = radixLSD;
		title = "LSD radix sort";
		if (argv[1][8]>='1' && argv[1][8]<='9'){
			const char *numberEnd;
			const size_t base = strtol(argv[1]+8, (char **)&numberEnd, 10);
			title += " base:";
			title.append(argv[1]+8, numberEnd);
			arbaseRadixLSD_Base = base;
			sortingFunction = arbaseRadixLSD;
		}
	} else{
		puts("wrong syntax\n");
		return 0;
	}

	title += ' ';
	title += std::to_string(std::size(arr));
	title += ' ';

	switch (letter){
	case 'U':{
			const float step = 1.f/std::size(arr);
			float sum = step;
			for (int i=0; i<(int)std::size(arr); ++i){
				arr[i] = sum;
				sum += step;
			}
			sp::Rand32 rng(clock());
			std::shuffle(std::begin(arr), std::end(arr), rng);
			title += "shuffled uniformly distrbuted values";
		} break;
	case 'P':{
			size_t steepness = 10;
			if (argc == 6){
				steepness = strtol(argv[5], nullptr, 10);
				if (!steepness){
					puts("steepnes of perlin noise cannot be 0");
					return 1;
				}
			}
			std::fill(std::begin(arr), std::end(arr), 0.f);
			sp::Rand32 rng(clock());
			float maxStep = 0.5f;
			for (size_t partitonSize=std::size(arr)/steepness+1; maxStep>=0.01f; maxStep/=2.f, partitonSize=partitonSize/2+1){
				std::uniform_real_distribution dist(0.f, maxStep);
				float p1 = dist(rng);
				for (size_t i=0; i<std::size(arr); i+=partitonSize){
					const float p2 = dist(rng);
					const float step = (p2 - p1) / partitonSize;
					for (size_t j=i; j<i+partitonSize && j<std::size(arr); ++j, p1+=step)
						arr[j] += p1;
					p1 = p2;
				}
			}
			title += "perlin noise of steepness ";
			title += std::to_string(steepness);
		} break;
	case 'S':{
			const float step = 1.f/std::size(arr);
			float sum = step;
			for (int i=0; i<(int)std::size(arr); ++i){
				arr[i] = sum;
				sum += step;
			}
			title += "sorted values";
		} break;
	case 'R':{
			const float step = 1.f/std::size(arr);
			float sum = 1.f;
			for (int i=0; i<(int)std::size(arr); ++i){
				arr[i] = sum;
				sum -= step;
			}
			title += "values sorted in reversed order";
		} break;
	case 'O':{
			const float step = 2.f / (std::size(arr) & ~1);
			float sum = step;
			int i = 0;
			const int halfSize = std::size(arr)/2;
			for (; i<halfSize; ++i){
				arr[i] = sum;
				sum += step;
			}
			if (std::size(arr) & 1){
				arr[i] = 1.f;
				++i;
			}
			for (; i<(int)std::size(arr); ++i){
				sum -= step;
				arr[i] = sum;
			}
			title += "organ pipes shaped values";
		} break;
	case 'V':{
			const float step = 2.f / (std::size(arr) & ~1);
			float sum = 1.f;
			int i = 0;
			const int halfSize = std::size(arr)/2;
			for (; i<halfSize; ++i){
				arr[i] = sum;
				sum -= step;
			}
			if (std::size(arr) & 1){
				arr[i] = step;
				++i;
			}
			for (; i<(int)std::size(arr); ++i){
				sum += step;
				arr[i] = sum;
			}
			title += "V shaped values";
		} break;
	case 'B':{
			size_t components = 5;
			float smoothness = 64.f;
			if (argc >= 6) components = strtol(argv[5], nullptr, 10);
			if (argc == 7) smoothness = strtof(argv[6], nullptr);
			if (smoothness<=0.f || components<=0){
				puts("Number of composite distributions or smoothness when specyfying the bell shape cannot be 0 or negative");
				return 1;
			}
			sp::Rand32 rng(clock());
			std::uniform_real_distribution<float> dist(0, (float)std::size(arr)/(float)components);
			std::fill(std::begin(arr), std::end(arr), 0.f);
			const uint32_t iterations = (smoothness + 0.5f) * std::size(arr);
			for (uint32_t i=0; i<iterations; ++i){
				float randomSum = 0.f;
				for (uint32_t j=0; j<components; ++j) randomSum += dist(rng);
				arr[(size_t)floor(randomSum)] += 1.f;
			}
			const float scalingVal = 1.f / *std::max_element(std::begin(arr), std::end(arr));
			for (auto &I : arr) I *= scalingVal;
			title += "values forming bell shape made of ";
			title += std::to_string(components);
			title += " uniform distrubutions";
		} break;
	case 1:{
			sp::Rand32 rng(clock());
			std::uniform_real_distribution dist(0.f, 1.f);
			for (auto &I : arr) I = dist(rng);
			title += "random values";
		} break;
	default:
		puts("Wrong argument specifying the input shape");
		return 1;
	}

	sf::RenderWindow windowObject(sf::VideoMode{sWidth, sHeight}, "", sf::Style::Titlebar|sf::Style::Close);
	window = &windowObject;
	window->setTitle(title);

	printValues();
	sortingFunction(arr.data(), arr.data()+std::size(arr));
	printValues();
	{
		sf::Event event;
		while (true){
			while (window->pollEvent(event)){
				if (event.type == sf::Event::Closed){
					window->close();
					return 0;
				}
				if (event.type == sf::Event::GainedFocus){
					printValues();
				}
			}
			usleep(100000);
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Escape) && window->hasFocus()){
				window->close();
				return 0;
			}
		}
	}
	return 0;
}