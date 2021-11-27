#pragma once

#include "sp_Utils.hpp"
#include <numeric>

namespace sp{ // BEGINING OF NAMESPACE //////////////////////////////////////////////////////////////////

template<class T>
void selectionSort(T *begin, T *end);

template<class T>
void bubbleSort(T *begin, T *end);

template<class T>
void insertionSort(T *begin, T *end);

template<class T>
void bogoSort(T *begin, T *end);






template<class T>
void selectionSort(T *begin, T *end){
	T *max;
	for (T *I=begin, *J; I!=end-1; ++I){
		for (J=I+1, max=I; J!=end; ++J)
			max = *max<*J ? max : J;
		sp::swap(I, max);
	}
}

template<class T>
void bubbleSort(T *begin, T *end){
	for (T *I; begin!=end; --end)
		for (I=begin+1; I<end; ++I)
			if (*I < *(I-1))
				sp::swap(I, I-1);
}

template<class T>
void insertionSort(T *begin, T *end){
	for (T *I=begin+1, *J; I!=end; ++I)
		for(J=I; J!=begin; --J){
			if (*(J-1) < *J)
				break;
			sp::swap(J, J-1);
		}
}

template<class T>
void bogoSort(T *begin, T *end){
	while (true){
		std::random_shuffle(begin, end);
		for (T *I=begin+1; I!=end; ++I)
			if (*I < *(I-1))
				return;
	}
}

template<class T>
void quickSort(T *begin, T *end){
	--end;
	if (begin >= end)
		return;
	T *border = begin;
	for (T *I=begin; I<end; ++I){
		if (*I < *end)
			swap(I, border++);
	}
	swap(border, end);
	quickSort(begin, border);
	quickSort(border+1, end+1);
}


}	// END OF NAMESPACE	///////////////////////////////////////////////////////////////////