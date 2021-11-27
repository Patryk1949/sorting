#pragma once

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <array>
#include <functional>
#include <numeric>


namespace sp{ // BEGINING OF NAMESPACE //////////////////////////////////////////////////////////////////

// CONVERTING UNIONS
union ConvF{
	int32_t i;
	float f;
};
union ConvD{
	int64_t i;
	double d;
};
template<class T>
union Conv8{
	static_assert(sizeof(T) == 1, "Type has wrong number of bytes !");
	int8_t i;
	T t;
};
template<class T>
union Conv16{
	static_assert(sizeof(T) == 2, "Type has wrong number of bytes !");
	int16_t i;
	T t;
};
template<class T>
union Conv32{
	static_assert(sizeof(T) == 4, "Type has wrong number of bytes !");
	int32_t i;
	T t;
};
template<class T>
union Conv64{
	static_assert(sizeof(T) == 8, "Type has wrong number of bytes !");
	int64_t i;
	T t;
};

template<class T>
inline constexpr void swap(T *const x, T *const y) noexcept;

template<class T>
inline constexpr auto sign(const T x);

template<class T>
inline constexpr int signInt(const T x);

inline constexpr int factorial(int n);

template<class T>
inline constexpr void clamp(const T &min, T *const x, const T &max);

template<class T>
inline constexpr T clamp(const T &min, const T &x, const T &max);

template<class T>
inline constexpr void clamp(const T &min, T *const x, const T &max, bool (*compare)(T, T));

template<class T>
inline constexpr T clamp(const T &min, T const &x, const T &max, bool (*compare)(T, T));


// TYPE ERASURES
template<class Signature>
class FunctionRef;

template<class Res, class... Args>
class FunctionRef<Res(Args...)>{
public:
	inline constexpr FunctionRef() noexcept : functionPtr{nullptr}, objectPtr{nullptr} {}
	FunctionRef(const FunctionRef &) = default;
	FunctionRef(FunctionRef &&) = default;
	FunctionRef &operator =(const FunctionRef &) = default;
	FunctionRef &operator =(FunctionRef &&) = default;
	~FunctionRef() = default;

	template<class Func>
	inline constexpr FunctionRef(Func *const func) noexcept;

	template<class Func>
	inline constexpr FunctionRef &operator =(Func *const func) noexcept;

	inline constexpr bool isStateless() const noexcept{ return !objectPtr; }

	inline constexpr Res operator ()(Args... args) noexcept;

private:
	void *functionPtr;
	void *objectPtr;
};


// FAST MATH
class Rand32{
public:
	inline constexpr Rand32() noexcept;
	inline constexpr Rand32(const uint32_t seed) noexcept;
	inline constexpr uint32_t min() const;
	inline constexpr uint32_t max() const;
	uint32_t operator ()() noexcept;

	typedef uint32_t result_type;
	uint32_t seed;
};

inline constexpr float qLog(const float x);

inline constexpr float qExp(const float x);

inline constexpr float qSqrt(const float x);

inline constexpr float qInvSqrt(const float x);

template<class T>
inline constexpr T intLog2(T x);


// HEAP
template<class It>
inline void repairHeap(It begin, const It end, const int startingIndex = 0);

template<class It, class Compare>
inline void repairHeap(It begin, const It end, Compare comapre, const int startingIndex = 0);

template<class It>
inline void makeHeap(It begin, const It end);

template<class It, class Compare>
inline void makeHeap(It begin, const It end, Compare compare);


// RANGE
namespace priv__{
	template<class T>
	class RangeClass;

	template<class T>
	class RangeRClass;
} // END OF NAMESPACE PRIV

template<class T>
priv__::RangeClass<T> Range(const T first, const T last, const T step = (T)1);

template<class T>
priv__::RangeRClass<T> RangeR(const T first, const T last, const T step = (T)1);






template<class T>
inline constexpr void swap(T *const x1, T *const x2) noexcept{
	const T temp = std::move(*x1);
	*x1 = std::move(*x2);
	*x2 = std::move(temp);
}

template<class T>
inline constexpr T sign(const T x){
	return (1-(signbit(x)<<1))*(x!=(T)0);
}

template<class T>
inline constexpr int signInt(const T x){
	return (1-(signbit(x)<<1))*(x!=(T)0);
}

inline constexpr int factorial(int n){
	for (int i=n-1; i>1; --i)
		n *= i;
	return n;
}

template<class T>
inline constexpr void clamp(const T &min, T *const x, const T &max){
	*x = *x<max ? (*x>min ? *x : min) : max;
}

template<class T>
inline constexpr T clamp(const T &min, const T &x, const T &max){
	return x<max ? (x>min ? x : min) : max;
}

template<class T>
inline constexpr void clamp(const T &min, T *const x, const T &max, bool (*compare)(T, T)){
	*x = compare(*x, max) ? (compare(*x, min) ? *x : min) : max;
}

template<class T>
inline constexpr T clamp(const T &min, T const &x, const T &max, bool (*compare)(T, T)){
	return compare(x, max) ? (compare(x, min) ? x : min) : max;	
}


// TYPE ERASURES
#define TPL template<class Res, class... Args>
#define CLS FunctionRef<Res(Args...)>

TPL template<class Func>
inline constexpr CLS::FunctionRef(Func *const func) noexcept{
	static_assert(std::is_invocable<Func, Args...>(), "Function's arguments don't match with thoose declared in template parameters.");
	if constexpr (std::is_function<Func>()){
		functionPtr = (void *)func;
		objectPtr = nullptr;
	
	} else{

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpmf-conversions"
		functionPtr = (void *)(Res (Func::*)(Args...))&Func::operator();
#pragma GCC diagnostic pop

		objectPtr = func;
	}
}

TPL template<class Func>
inline constexpr CLS &CLS::operator =(Func *const func) noexcept{
	if constexpr (std::is_function<Func>()){
		objectPtr = nullptr;
		functionPtr = (void *)func;
	}
	if constexpr (std::is_class<Func>()){

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpmf-conversions"
		functionPtr = (void *)(Res (Func::*)(Args...))&Func::operator();
#pragma GCC diagnostic pop

		objectPtr = func;
	}
}

TPL inline constexpr Res CLS::operator ()(Args... args) noexcept{
	if (objectPtr) 
		return ((Res (*)(void *, Args...))functionPtr)(objectPtr, args...);
	else
		return ((Res (*)(Args...))functionPtr)(args...);
}

#undef CLS
#undef TPL


// FAST MATH
inline constexpr Rand32::Rand32() noexcept : seed{0} {}

inline constexpr Rand32::Rand32(const uint32_t firstSeed) noexcept : seed{firstSeed} {}

inline constexpr uint32_t Rand32::min() const{
	return 0;
}

inline constexpr uint32_t Rand32::max() const{
	return std::numeric_limits<uint32_t>::max();
}

uint32_t Rand32::operator ()() noexcept{
	seed += 0xe120fc15;
	uint64_t temp = (uint64_t)seed * 0x4a39b70d;
	temp = (uint64_t)((temp >> 32) ^ temp) * 0x12fad5c9;
	return (temp >> 32) ^ temp;
}


inline constexpr float qLog(const float x){
	constexpr float scaleDown = 1.f/(float)0x00800000;
	constexpr ConvF one{.f = 1.f};
	return (float)(*(int32_t *)&x - one.i)*scaleDown;
}

inline constexpr float qExp(const float x){
	constexpr float scaleUp = 0x00800000;
	constexpr ConvF one{.f = 1.f};
	const ConvF res{.i = (int32_t)(x*scaleUp) + one.i};
	return res.f;
}

inline constexpr float qPow(const float x, const float y){
	constexpr ConvF one{.f = 1.f};
	ConvF res{.i = (int32_t)((float)(*(int32_t *)&x - one.i)*y) + one.i};
	return res.f;
}

inline constexpr float qSqrt(const float x){
	ConvF res{.f = x};
	res.i -= 1<<23;
	res.i >>= 1;
	res.i += 1<<29;
	return (res.f + x/res.f)*0.5f;
}

inline constexpr float qInvSqrt(const float x){
	ConvF res{.i = (int32_t)0x5f3759df - (*(int32_t *)&x>>1)};
	return res.f * (1.5f - 0.5f*x*res.f*res.f);
}


template<class T>
inline constexpr T intLog2(const T x){
	static_assert(std::is_integral<T>(), "Argument to intLog2 must be a integer.");
	typedef std::make_unsigned_t<T> UnsignedT;

	UnsignedT valueToShift = x - 1;
	T result = 0;

	for (;valueToShift & ((UnsignedT)std::numeric_limits<T>::max()<<8); result+=8, valueToShift>>=8);
	for (;valueToShift; ++result, valueToShift>>=1);
	return result;
}


// HEAP
template<class It>
inline void repairHeap(It begin, const It end, const int startingIndex){
	It child;
	for(int i=startingIndex; (child=begin+((i<<1)|1)) < end; i=child-begin){
		child += child[0] < child[child+1 != end];
		
		if (child[0] < begin[i])
			return;

		sp::swap(begin+i, child);
	}
}

template<class It, class Compare>
inline void repairHeap(It begin, const It end, Compare comapre, const int startingIndex){
	It child;
	for(int i=startingIndex; (child=begin+((i<<1)|1)) < end; i=child-begin){
		child += child[0] < child[child+1 != end];
		
		if (child[0] < begin[i])
			return;

		sp::swap(begin+i, child);
	}
}

template<class It>
inline void makeHeap(It begin, const It end){
	for (int i=(end-begin-1)>>1; i>=0; --i)
		repairHeap(begin, end, i);
}

template<class It, class Compare>
inline void makeHeap(It begin, const It end, Compare compare){
	for (int i=(end-begin-1)>>1; i>=0; --i)
		repairHeap(begin, end, compare, i);
}


// RANGES
namespace priv__{

template<class T>
class RangeClass{
public:
	class It : public std::iterator<std::forward_iterator_tag, T>{
	public:
		inline constexpr T operator *() noexcept{ return i; }
		inline constexpr const T &operator ++() { i += step; return i; }
		inline constexpr bool operator !=(const T rhs){ return i <= rhs; }

		T i;
		T step;
	};

	inline constexpr It begin(){
		return It{{}, firstElement, stepSize};
	}
	inline constexpr T end(){ 
		return lastElement;
	}

	T firstElement;
	T lastElement;
	T stepSize;
};

template<class T>
class RangeRClass{
public:
	class It : public std::iterator<std::forward_iterator_tag, T>{
	public:
		inline constexpr T operator *() noexcept{ return i; }
		inline constexpr const T &operator ++() { i -= step; return i; }
		inline constexpr bool operator !=(const T rhs){ return i >= rhs; }

		T i;
		T step;
	};

	inline constexpr It begin(){
		return It{{}, firstElement, stepSize};
	}
	inline constexpr T end(){ 
		return lastElement;
	}

	T firstElement;
	T lastElement;
	T stepSize;
};

} // END OF NAMESPACE PRIV

template<class T>
priv__::RangeClass<T> Range(const T first, const T last, const T step){
	static_assert(std::is_integral<T>() || std::is_floating_point<T>(), "Range must be made of numbers.");
	return priv__::RangeClass<T>{first, last, step};
};

template<class T>
priv__::RangeRClass<T> RangeR(const T first, const T last, const T step){
	static_assert(std::is_integral<T>() || std::is_floating_point<T>(), "Range must be made of numbers.");
	return priv__::RangeRClass<T>{first, last, step};
};

}	// END OF NAMESPACE	///////////////////////////////////////////////////////////////////