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


class RodCutting
{
public:

	//
	// O( 2^n )
	// n represents the length, so pPrices is expected to have the form [0, p1, p2, ..., p(n-1)]
	//
	static double ComputeTopDown(int n, double* pPrices)
	{
		if (n == 0)
		{
			return 0;
		}

		double maxPrice = 0;

		for (int i = 1; i <= n; i++)
		{
			maxPrice = max(maxPrice, pPrices[i] + ComputeTopDown(n-i, pPrices));
		}

		return maxPrice;
	}

	//
	// O(n^2)
	//
	static double ComputeBottomUp(const vector<double>& _prices, vector<int>* pSolution = 0)
	{
		vector<double> prices;

		if (_prices.size() == 0)
		{
			return 0;
		}

		//
		// the data is expected to have the form [0, p1, p2, ..., p(n-1)]
		//
		if (!FLOAT_EQUAL_ZERO(_prices[0]))
		{
			prices.push_back(0);
		}

		prices.insert(prices.end(), _prices.begin(), _prices.end());

		int n = prices.size()-1;

		vector<double> maxPrices(n+1, 0);
		vector<int> memos(n+1, -1);

		for (int j = 1; j <= n; j++)
		{
			double maxPrice = 0;

			for (int i = 1; i <= j; i++)
			{
				double p = prices[i] + maxPrices[j - i];

				if (maxPrice < p)
				{
					maxPrice = p;
					memos[j] = i;  // store the size of the first piece
				}
			}

			maxPrices[j] = maxPrice;
		}

		if (pSolution != 0)
		{
			int length = 0;

			int k = n;

			while (k > 0)
			{
				length += memos[k];

				pSolution->push_back(memos[k]);

				k = k - memos[k];
			}
		}

		return maxPrices[n];
	}

	static void Test()
	{
		vector<double> prices = { 0, 3, 5, 10, 12, 14 };
		
		vector<int> solution;

		printf("RodCutting { 0, 3, 5, 10, 12, 14 } max price is %lf \n", ComputeBottomUp(prices, &solution));

		printf("RodCutting { 0, 3, 5, 10, 12, 14 } max price is %lf \n", ComputeTopDown(prices.size()-1, &prices.front()));
	}
};


