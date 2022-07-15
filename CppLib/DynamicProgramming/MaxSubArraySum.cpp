#include <algorithm>
#include <cstdio>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <list>
#include <limits>

using namespace std;

#define EPS (1e-9)

#define FLOAT_EQUAL(a, b)				(fabs(a-b)<=EPS)
#define FLOAT_EQUAL_ZERO(a)				(fabs(a) <= EPS)
#define FLOAT_GREATER_EQUAL_ZERO(a)		(a > -EPS)
#define FLOAT_LESS_EQUAL_ZERO(a)		(a < EPS)


///////////////////////////////////////////////////////////////////////////////////////////


class MaxSubArraySum
{
public:

	//
	// Kadane's Algorithm, enhanced to handle all negatives and zeros
	//
	static int KadaneAlgorithm(int* p, int size, int* pStart = 0, int* pEnd = 0)
	{
		if (p == 0 || size <= 0)
		{
			return 0;
		}

		int _start, _end;
		int& start = (pStart != 0) ? *pStart : _start;
		int& end = (pEnd != 0) ? *pEnd : _end;

		start = end = 0;

		int max_so_far = -numeric_limits<int>::max();
		int max_ending_here = 0;

		int left = 0, right = 0;
		
		bool allNegative = true;
		int largestNegative = -numeric_limits<int>::max();
		int largestNegativePos = 0;

		for (int i = 0; i < size; i++)
		{
			int v = p[i];

			if (v <= 0)
			{
				left = right = i + 1;

				if (largestNegative < v)
				{
					largestNegative = v;
					largestNegativePos = i;
				}
			}
			else
			{
				allNegative = false;
			}

			max_ending_here += v;
			
			if (max_ending_here < 0)
			{
				max_ending_here = 0;
			}

			if (max_so_far < max_ending_here)
			{
				max_so_far = max_ending_here;

				right = i;

				start = left;
				end = right;
			}
		}

		if (allNegative)
		{
			max_so_far = largestNegative;
			start = end = largestNegativePos;
		}

		return max_so_far;
	}

	static void Test()
	{
		int p[] = { -2, -3, 4, -1, -2, 1, 5, -3 };
		int n = sizeof(p) / sizeof(p[0]);

		int start = 0, end = 0;

		int sum = KadaneAlgorithm(p, n, &start, &end);

		printf(" MaxSubArraySum of { -2, -3, 4, -1, -2, 1, 5, -3 } is %d [%d-%d]\n", sum, start, end);

		p[2] = -4;
		p[5] = -1;
		p[6] = -5;

		sum = KadaneAlgorithm(p, n, &start, &end);

		printf(" MaxSubArraySum of { -2, -3, -4, -1, -2, -1, -5, -3 } is %d [%d-%d]\n", sum, start, end);

		p[2] = 0;

		sum = KadaneAlgorithm(p, n, &start, &end);

		printf(" MaxSubArraySum of { -2, -3, 0, -1, -2, -1, -5, -3 } is %d [%d-%d]\n", sum, start, end);

	}
};


