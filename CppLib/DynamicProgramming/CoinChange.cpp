#include <algorithm>
#include <cstdio>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <list>
#include <limits>

using namespace std;



///////////////////////////////////////////////////////////////////////////////////////////

//
// Minimum Number of Coins:
// Given N types of coin, each type with a value and infinite quantity, determine what coins to include 
// so that the total value amounts to TOTAL_VALUE using a minimum total number of coins.
//
// Total Ways of Using Coins:
// Given N types of coin, each type with a value and infinite quantity, determine how many ways to choose coins 
// so that the total value amounts to TOTAL_VALUE.
//


class CoinChange
{
public:

	static int MinimumNumberOfCoins_Backtrack_1D(int TOTAL_VALUE, int* VALUE, int N, vector<int>* pSolution = 0)
	{
		if (TOTAL_VALUE <= 0 || N <= 0)
		{
			return 0;
		}

		//
		// min_num[v] is the minimum number of coins to make up value v.
		//
		// if coin n is used, we would have 
		//
		// min_num[v] = 1 + min_num[v - VALUE[n]]
		//
		// the same logic applies for VALUE[1, 2, ..., n-1], min_num[v] should be the smallest of all
		// 
		vector<int> min_num(TOTAL_VALUE + 1, INT_MAX);

		//
		// memos[v] stores the coin for the solution
		//
		vector<int> memos(TOTAL_VALUE + 1, -1);

		min_num[0] = 0;

		for (int v = 1; v < (TOTAL_VALUE + 1); v++)
		{
			for (int n = 0; n < N; n++)
			{
				if (v >= VALUE[n] && (1 + min_num[v - VALUE[n]]) < min_num[v])
				{
					min_num[v] = 1 + min_num[v - VALUE[n]];
					memos[v] = n;
				}
			}
		}

		if (pSolution != 0)
		{
			pSolution->clear();

			int v = TOTAL_VALUE;
			while( v > 0)
			{
				int k = memos[v];

				pSolution->push_back(k);

				v = v - VALUE[k];
			}
		}
	
		return min_num[TOTAL_VALUE];
	}

	static int TotalWaysOfUsingCoins_Backtrack_2D(int TOTAL_VALUE, int* VALUE, int N)
	{
		if (TOTAL_VALUE < 0 || N <= 0)
		{
			return 0;
		}

		if (TOTAL_VALUE == 0)
		{
			return 1;
		}

		//
		// ways[n][v] is the total ways we can sum up coin [0...n-1] to make up value v.
		//
		// solution for n can be viewed as the sum of two cases:
		// 
		// 1. excluding coin n-1
		// 2. including coin n-1
		//
		// ways[n][v] = ways[n-1][v] + ways[n][v - VALUE[n-1]]
		//
		vector< vector<int> > ways(N, vector<int>(TOTAL_VALUE + 1, 0));

		// if required total value is 0 then there is 1 solution, that is, do not include any coin
		for (int n = 0; n < N; n++) {
			ways[n][0] = 1;		
		}

		for (int v = 1; v < (TOTAL_VALUE + 1); v++)
		{
			for (int n = 0; n < N; n++)
			{
				int including = (v - VALUE[n] >= 0) ? ways[n][v - VALUE[n]] : 0;
				int excluding = (n >= 1) ? ways[n - 1][v] : 0;

				ways[n][v] = including + excluding;
			}
		}

		return ways[N-1][TOTAL_VALUE];
	}

	static void Test()
	{
		int VALUE[] = { 25, 10, 5, 1 };
		int TOTAL_VALUE = 67;
		int N = sizeof(VALUE) / sizeof(VALUE[0]);

		vector<int> solution;

		printf("MinimumNumberOfCoins_Backtrack_1D: should be 6, the answer is %d\n", CoinChange::MinimumNumberOfCoins_Backtrack_1D(TOTAL_VALUE, VALUE, N, &solution) ); // = 6


		int VALUE2[] = { 2, 5, 3, 6 };
		int TOTAL_VALUE2 = 10;
		int N2 = sizeof(VALUE2) / sizeof(VALUE2[0]);

		printf("TotalWaysOfUsingCoins_Backtrack_2D: should be 5, the answer is %d\n", CoinChange::TotalWaysOfUsingCoins_Backtrack_2D(TOTAL_VALUE2, VALUE2, N2)); // = 5
	}
};


