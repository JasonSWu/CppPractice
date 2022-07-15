#include <algorithm>
#include <cstdio>
#include <vector>
#include <queue>
#include <map>
#include <limits>

using namespace std;


///////////////////////////////////////////////////////////////////////////////////////////


#define INF (1e9)
#define EPS (1e-9)
#define PI (acos(-1.0)) 

#define SMALLEST_FLOAT				(numeric_limits<float>::lowest())
#define SMALLEST_POSITIVE_FLOAT		(numeric_limits<float>::min())
#define LARGEST_FLOAT				(numeric_limits<float>::max())
#define SMALLEST_DOUBLE				(numeric_limits<double>::lowest())
#define SMALLEST_POSITIVE_DOUBLE	(numeric_limits<double>::min())
#define LARGEST_DOUBLE				(numeric_limits<double>::max())

#define SMALLEST_INT32				(numeric_limits<int>::min())
#define LARGEST_INT32				(numeric_limits<int>::max())
#define SMALLEST_INT64				(numeric_limits<long long>::min())
#define LARGEST_INT64				(numeric_limits<long long>::max())


#define FLOAT_EQUAL(a, b)				(fabs(a-b)<=EPS)
#define FLOAT_EQUAL_ZERO(a)				(fabs(a) <= EPS)
#define FLOAT_GREATER_EQUAL_ZERO(a)		(a > -EPS)
#define FLOAT_LESS_EQUAL_ZERO(a)		(a < EPS)

#define DEG_to_RAD(d)  (d * PI / 180.0)
#define RAD_to_DEG(r)  (r * 180.0 / PI)

#ifndef NULL
	#define NULL 0
#endif


#define LEAST_SIGNIFICANT_BIT(x) (x & (-x))

///////////////////////////////////////////////////////////////////////////////////////////

//
// This is the "Predicate" by way of operator "()" to check if l is greater than r
//
struct IntGreaterThanCompare
{
	bool operator()(const int& l, const int& r)
	{
		return l > r;
	}
};

//
// GreaterThanCompare is used so that priority_queue.top() returns the smallest instead
// 
typedef priority_queue<int> IntMaxHeap;
typedef priority_queue<int, vector<int>, IntGreaterThanCompare> IntMinHeap;


//
//
//
struct StringLessThanCompare
{
	char* T;

	StringLessThanCompare(char* t)
	{
		T = t;
	}

	bool operator()(int a, int b)
	{
		return strcmp(T + a, T + b) < 0;
	}

	bool IsLess(int a, int b)
	{
		return strcmp(T + a, T + b) < 0;
	}
};

//
//	char T[100];
//	int SA[100];
//	
//	StringLessThanCompare SLTC(T);
//	
//	sort(SA, SA + 99, GTC); 
//	
//	or 
//	
//	sort(SA, SA + 99, GTC.IsLess);
//
