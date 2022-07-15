#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>

using namespace std;



///////////////////////////////////////////////////////////////////////////////////////////


//
// Detect if a subset from a given set of N non-negative integers sums upto a given value S
//
// Reference: https://www.youtube.com/watch?v=5td2QH-x5ck
//

class SubsetSumDetection
{
public:

	//int S = 5;
	//vector<int> A = { 1, 3, 9, 2 }; // => True

	int S = 11;
	vector<int> A = { 2, 3, 7, 8, 10 }; // => True

	void ComputeBottomUp()
	{
		int N = (int) A.size();

		A.insert(A.begin(), 0); // inject a 0 at the beginning to make the computation easier

		vector<vector<bool>> dp;

		dp.assign(N + 1, vector<bool>(S + 1, false));

		//
		// dp[n][sum] indicates if any subset formed by numbers from A[1] to A[n] can add up to "sum" 
		//
		// dp[n][sum] will be:
		//
		//	1. same as dp[n-1][sum], meaning even without using "n", dp[n-1][sum] might already be able to add up to "sum"
		//  2. in addition, consider dp[n-1][sum-n], meaning with "n" in place, it may reach "sum" depending on the "sum-n" case 
		//

		dp[0][0] = true; // Empty set adds to zero, the solution is an empty set, so set this element to true.
						// any other [0][sum] is already set to false.

		for (int n = 1; n <= N; ++n)
		{
			//
			// 1. at least dp[n-1][sum]
			//
			for (int sum = 0; sum <= S; ++sum) {
				dp[n][sum] = dp[n - 1][sum];
			}

			//
			// 2. plus dp[n-1][sum-n]
			//
			for (int sum = 0; sum <= S; ++sum) {
				if ((sum - n) >= 0) {
					dp[n][sum] = dp[n][sum] | dp[n - 1][sum - n];
				}
			}
		}

		printf("SubsetSumDetection, dp[%d][%d] = %s \n", N, S, dp[N][S] ? "True" : "False");
	}

	static void Test() {
		SubsetSumDetection o;
		o.ComputeBottomUp();
	}
};

