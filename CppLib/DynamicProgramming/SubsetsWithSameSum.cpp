#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>

using namespace std;



///////////////////////////////////////////////////////////////////////////////////////////


//
// For many sets of consecutive integers from 1 through N (1 <= N <= 39), 
// partition the set into two sets whose sums are identical. 
// Given N, your program should print the number of ways a set containing the integers from 1 through N can be
// partitioned into two sets whose sums are identical. Print 0 if there are no such ways.
//
// Reference: USACO Training Section 2.2
//
class SubsetsWithSameSum
{
public:

	int N = 31;
	long long count = 0;

	int half_total;

	vector<vector<long long>> dp;

	void ComputeBottomUp()
	{
		int total = N*(N + 1) / 2;
		int half_total = total / 2;

		if (total == half_total * 2) // it only makes sense when the sum can be divided by 2
		{
			//ComputeTopDown(0, 0);

			dp.assign(N + 1, vector<long long>(total + 1, 0));

			//
			// dp[n][sum] stores the ways of adding up to a "sum" from any subset formed by numbers up to "n" [1...n] 
			//
			// dp[n][sum] = dp[n-1][sum] + dp[n-1][sum-n]
			//
			//	1. at least dp[n-1][sum], meaning even without using "n", there are already dp[n-1][sum] ways to add up to "sum"
			//  2. plus dp[n-1][sum-n], meaning with "n" in place, it can make new ways to add up to "sum" together with the "sum-n" cases 
			//

			dp[0][0] = 1; // Empty set adds to zero, the solution is an empty set, so set this element to 1.
						  // any other [0][sum] is already set to zero.

			for (int n = 1; n <= N; ++n)
			{
				//
				// 1. at least dp[n-1][sum]
				//
				for (int sum = 0; sum <= half_total; ++sum) {
					dp[n][sum] = dp[n - 1][sum];
				}

				//
				// 2. plus dp[n-1][sum-n]
				//
				for (int sum = 0; sum <= half_total; ++sum) {
					if ((sum - n) >= 0) {
						dp[n][sum] += dp[n - 1][sum - n];
					}
				}
			}

			count = dp[N][half_total];
		}

		count /= 2;

		//
		// N = 7, count = 4
		// N = 31, count = 8273610
		//
		printf("SubsetsWithSameSum, when N = %d, count = %lld \n", N, count);
	}


	// this recursive version is too slow
	bool ComputeTopDown(int sum, int pos) // true means continue to search
	{
		if (sum == half_total) {
			count++;
			return false;
		}

		if (sum > half_total) {
			return false;
		}

		for (int i = pos + 1; i <= N; i++) {
			if (!ComputeTopDown(sum + i, i)) {
				break;
			}
		}

		return true;
	}

	static void Test() {
		SubsetsWithSameSum o;
		o.ComputeBottomUp();
	}
};

