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


//
// 0-1 Knapsack problem:
// Given N types of items, each type with a capacity (or cost, or weight) and a value, one item for each type, determine which items to include 
// so that the total capacity (or cost, or weight) is less than or equal to MAX_CAPACITY and the total value is as large as possible.
//
// Complete Knapsack problem (unbounded):
// Given N types of items, each type  with a capacity (or cost, or weight), a value, and infinite quantity, determine which items to include  
// so that the total capacity (or cost, or weight) is less than or equal to MAX_CAPACITY and the total value is as large as possible.
//
// Multiple Knapsack problem (unbounded):
// Given N types of items, each type with a capacity (or cost, or weight), a value, and specified quantity, determine which items to include  
// so that the total capacity (or cost, or weight) is less than or equal to MAX_CAPACITY and the total value is as large as possible.
//
// Exact Knapsack problem:
// Given N types of items, each type with a value, one item for each type, determine which items to include  
// so that the total value is equal to TOTAL_VALUE.
//
// Reference: 
//
// http://www.es.ele.tue.nl/education/5MC10/Solutions/knapsack.pdf
//
//
// http://www.hawstein.com/posts/dp-knapsack.html
//

class Knapsack
{
public:

	//
	// Time: O(2^N)
	//
	static int ZeroOnePack_Recursive(int MAX_CAPACITY, int* VALUE, int* CAPACITY, int n_item)
	{
		if (MAX_CAPACITY <= 0 || n_item <= 0)
		{
			return 0;
		}

		if (CAPACITY[n_item - 1] > MAX_CAPACITY) {
			return ZeroOnePack_Recursive(MAX_CAPACITY, VALUE, CAPACITY, n_item - 1);
		}

		int w1 = ZeroOnePack_Recursive(MAX_CAPACITY, VALUE, CAPACITY, n_item - 1);
		int w2 = VALUE[n_item - 1] + ZeroOnePack_Recursive(MAX_CAPACITY - CAPACITY[n_item - 1], VALUE, CAPACITY, n_item - 1);

		return max(w1, w2);
	}

	//
	// Time: O( N * MAX_CAPACITY )
	// Space: O( N * MAX_CAPACITY )
	//
	static int ZeroOnePack_Backtrack_1D(int MAX_CAPACITY, int* VALUE, int* CAPACITY, int N)
	{
		if (MAX_CAPACITY <= 0 || N <= 0)
		{
			return 0;
		}

		//
		// value[c] stores the max weight of all items with combined capacity at most c (c <= MAX_CAPACITY)
		//
		// Initialization: make sure weight[0...MAX_CAPACITY] is set to zero
		//
		vector<int> value(MAX_CAPACITY + 1, 0);

		//
		// Initializ the case that only the first item is picked
		//

		value[0] = 0;

		for (int c = 1; c < (MAX_CAPACITY + 1); c++) {
			value[c] = (CAPACITY[0] <= c) ? VALUE[0] : 0;
		}

		for (int n = 1; n < N; n++)
		{
			for (int c = MAX_CAPACITY; c >= CAPACITY[n]; c--)	// this loop is reversed compared to the 2D version!
			{
				value[c] = max(value[c], value[c - CAPACITY[n]] + VALUE[n]);
			}
		}

		return value[MAX_CAPACITY];
	}

	//
	// Time: O( N * MAX_CAPACITY )
	// Space: O( MAX_CAPACITY )
	//
	static int ZeroOnePack_Backtrack_2D(int MAX_CAPACITY, int* VALUE, int* CAPACITY, int N, vector<int>* pSolution = 0)
	{
		if (MAX_CAPACITY <= 0 || N <= 0) {
			return 0;
		}

		//
		// value[n][c] stores the max value of the first n-th items {1, 2, 3, ..., n} (n <= N) with combined capacity at most c (c <= MAX_CAPACITY)
		//
		// value[n][c] = max( value[n-1][c], value[n-1][c-CAPACITY[n]] + VALUE[n] );
		//
		vector< vector<int> > value(N + 1, vector<int>(MAX_CAPACITY + 1, 0));

		for (int n = 0; n <= N; n++)
		{
			for (int c = 0; c <= MAX_CAPACITY; c++)
			{
				if (n == 0 || c == 0)
				{
					value[n][c] = 0;
					continue;
				}

				if (CAPACITY[n - 1] <= c) {
					value[n][c] = max(VALUE[n - 1] + value[n - 1][c - CAPACITY[n - 1]], value[n - 1][c]);
				}
				else {
					value[n][c] = value[n - 1][c];
				}
			}
		}

		if (pSolution != 0)
		{
			pSolution->clear();

			//
			// if memos[n][c] is true, then n is a pick, we can then repeat 
			// this argument for memos[n-1][c].
			// if memos[n][c] is false, then we can then repeat 
			// this argument for memos[n-1][c] and so on.
			//
			int c = MAX_CAPACITY;
			for (int n = N; n > 0; n--)
			{
				if (value[n][c] > value[n - 1][c])
				{
					pSolution->push_back(n - 1);
					c = c - CAPACITY[n - 1];
				}
			}
		}

		return value[N][MAX_CAPACITY];
	}

	//
	// Time: O( N * SUM( MAX_CAPACITY/C[n] ) )
	// Space: same as Time
	//
	static int CompletePack_Backtrack_2D(int MAX_CAPACITY, int* VALUE, int* CAPACITY, int N, vector<int>* pSolution = 0)
	{
		if (MAX_CAPACITY <= 0 || N <= 0)
		{
			return 0;
		}

		//
		// value[n][c] stores the max value of the first n-th items {0, 1, 2, 3, ..., n-1} (n < N) with combined capacity at most c (c <= MAX_CAPACITY)
		//
		// value[n][c] = max{ value[n-1][c-k*CAPACITY[n]] + k*VALUE[n] | 0 <= k*CAPACITY[n] <= MAX_CAPACITY};
		//
		vector< vector<int> > value(N, vector<int>(MAX_CAPACITY + 1, 0));

		value[0][0] = 0;

		//
		// Initializ the case that only the first item is picked
		//
		for (int c = 1; c < (MAX_CAPACITY + 1); c++) {
			value[0][c] = VALUE[0] * (c / CAPACITY[0]);
		}

		for (int n = 1; n < N; n++)
		{
			for (int c = 1; c < (MAX_CAPACITY+1); c++)
			{
				int v = value[n - 1][c];

				for (int k = 1; k <= c / CAPACITY[n]; k++)
				{
					int v2 = value[n - 1][c - k * CAPACITY[n]] + k * VALUE[n];

					v = max(v, v2);
				}

				value[n][c] = v;
			}
		}

		return value[N-1][MAX_CAPACITY];
	}

	//
	// Time: O( N * SUM( LOG(MAX_CAPACITY/C[n]) ) )
	// Space: same as Time
	//
	static int CompletePack_Backtrack_ConvertToZeroOne_2D(int MAX_CAPACITY, int* _VALUE, int* _CAPACITY, int N, vector<int>* pSolution = 0)
	{
		if (MAX_CAPACITY <= 0 || N <= 0)
		{
			return 0;
		}

		vector<int> VALUE, CAPACITY;

		vector< pair<int, int> > memo; // memorize how the new VALUE and CAPACITY map to orginal _VALUE and _CAPACITY

		for (int n = 0; n < N; n++)
		{
			int factor = 1;

			while (true)
			{
				if (factor * _CAPACITY[n] > MAX_CAPACITY) {
					break;
				}

				VALUE.push_back(factor * _VALUE[n]);
				CAPACITY.push_back(factor * _CAPACITY[n]);

				memo.push_back(pair<int, int>(n, factor));

				factor *= 2;
			}
		}

		vector<int> solution;

		int results = ZeroOnePack_Backtrack_2D(MAX_CAPACITY, &VALUE.front(), &CAPACITY.front(), (int)VALUE.size(), &solution);

		if (pSolution != 0)
		{
			for (int i = 0; i < (int)solution.size(); i++)
			{
				pair<int, int>& p = memo[solution[i]];

				for(int j=0; j < p.second; j++ ) {
					pSolution->push_back(p.first);
				}
			}
		}

		return results;
	}

	//
	// Time: O( N * MAX_CAPACITY )
	// Space: O( MAX_CAPACITY )
	//
	static int CompleteBacktrack_1D(int MAX_CAPACITY, int* VALUE, int* CAPACITY, int N)
	{
		if (MAX_CAPACITY <= 0 || N <= 0)
		{
			return 0;
		}

		//
		// value[c] stores the max value of all items with combined capacity at most c (c <= MAX_CAPACITY)
		//
		vector<int> value(MAX_CAPACITY + 1, 0);

		value[0] = 0;

		//
		// Initializ the case that only the first item is picked
		//
		for (int c = 1; c < (MAX_CAPACITY + 1); c++) {
			value[c] = VALUE[0] * (c / CAPACITY[0]);
		}

		for (int n = 1; n < N; n++)
		{
			for (int c = CAPACITY[n]; c < (MAX_CAPACITY + 1); c++)	// this loop is reversed compared to 0-1 knapsack!
			{
				value[c] = max(value[c], value[c - CAPACITY[n]] + VALUE[n]);
			}
		}

		return value[MAX_CAPACITY];
	}

	//
	// Time: O( MAX_CAPACITY * SUM( QUANTITY[i] ) )
	// Space: same as Time
	//
	static int MultiplePack_Backtrack_ConvertToZeroOne_2D(int MAX_CAPACITY, int* _VALUE, int* _CAPACITY, int* QUANTITY, int N, vector<int>* pSolution = 0)
	{
		if (MAX_CAPACITY <= 0 || N <= 0)
		{
			return 0;
		}

		vector<int> VALUE, CAPACITY;

		vector< pair<int, int> > memo;	// memorize how the new VALUE and CAPACITY map to orginal _VALUE and _CAPACITY

		for (int n = 0; n < N; n++)
		{
			vector<int> factors;

			int factor = 1;
			int q = QUANTITY[n];

			while (factor < q)		// factor are 2^0, 2^1, ..., 2^(k-1) and QUANTITY[n] - 2^k + 1 so that the sum of these factors is QUANTITY[n]
			{
				if (factor * _CAPACITY[n] > MAX_CAPACITY) {
					break;
				}

				factors.push_back(factor);

				q -= factor;
				factor *= 2;
			}

			if (q > 0) {
				factors.push_back(q);
			}

			for (auto it = factors.begin(); it != factors.end(); it++)
			{
				VALUE.push_back( (*it) * _VALUE[n]);
				CAPACITY.push_back( (*it) * _CAPACITY[n]);

				memo.push_back(pair<int, int>(n, (*it)));
			}
		}

		vector<int> solution;

		int results = ZeroOnePack_Backtrack_2D(MAX_CAPACITY, &VALUE.front(), &CAPACITY.front(), (int)VALUE.size(), &solution);

		if (pSolution != 0)
		{
			for (int i = 0; i < (int)solution.size(); i++)
			{
				pair<int, int>& p = memo[solution[i]];

				for (int j = 0; j < p.second; j++) {
					pSolution->push_back(p.first);
				}
			}
		}

		return results;
	}

	//
	// Time: O( N * TOTAL_VALUE )
	// Space: O( N * TOTAL_VALUE )
	//
	static bool ExactPack_Backtrack_2D(int TOTAL_VALUE, int* VALUE, int N, vector<int>* pSolution = 0)
	{
		if (TOTAL_VALUE <= 0 || N <= 0)
		{
			return 0;
		}

		//
		// flags[n][v] indicates if there is a way to choose the first n items ( 0 <= n <= N) so that
		// the combined value equals to v (v <= TOTAL_VALUE)
		//
		// flags[n][v] = 
		//
		// 1. flags[n-1][v], if flags[n-1][v] is true
		// 2. flags[n-1][v-VALUE[n]], if flags[n-1][v] is false
		//
		vector< vector<bool> > flags(N+1, vector<bool>(TOTAL_VALUE + 1, false));

		//
		// memos[n][v] indicates if the solution includes n-th item
		//
		vector< vector<bool> > memos(N+1, vector<bool>(TOTAL_VALUE + 1, false));

		for (int v = 0; v < (TOTAL_VALUE + 1); v++)
		{
			flags[0][v] = false;				// "n == 0" means no item is picked, so the solution is false 
			memos[0][v] = flags[0][v];

			flags[1][v] = (VALUE[0] == v);		// for "n == 1", there is a solution iff the first item has the right value 
			memos[1][v] = flags[1][v];
		}

		flags[0][0] = true;

		for (int n = 2; n <= N; n++)
		{
			int k = n - 1;

			for (int v = 0; v < (TOTAL_VALUE + 1); v++)
			{
				if (flags[n - 1][v])	// Case 1
				{
					flags[n][v] = true;
					memos[n][v] = false;
				}
				else if ((v - VALUE[k]) >= 0 && flags[n - 1][v - VALUE[k]])		// Case 2
				{
					flags[n][v] = true;
					memos[n][v] = true;
				}
			}
		}

		if (pSolution != 0)
		{
			pSolution->clear();

			//
			// if memos[n][v] is true, then n is a pick, we can then repeat 
			// this argument for memos[n-1][v].
			// if memos[n][v] is false, then we can then repeat 
			// this argument for memos[n-1][v] and so on.
			//
			int v = TOTAL_VALUE;
			for (int n = N; n > 0; n--)
			{
				if (memos[n][v])
				{
					int k = n - 1;

					pSolution->push_back(k);

					v = v - VALUE[k];
				}
			}
		}

		return flags[N][TOTAL_VALUE];
	}

public:

	static void Test_ZeroOnePack()
	{
		vector<int> solution;

		int VALUE1[] = { 9, 6, 1, 4, 1 };
		int CAPACITY1[] = { 4, 3, 5, 2, 5 };
		int MAX_CAPACITY1 = 10;
		int N1 = sizeof(CAPACITY1) / sizeof(CAPACITY1[0]);


		int res1 = ZeroOnePack_Backtrack_1D(MAX_CAPACITY1, VALUE1, CAPACITY1, N1); // = 19
		printf("ZeroOnePack_Backtrack_1D: 19 -- %d\n", res1);
		res1 = ZeroOnePack_Backtrack_2D(MAX_CAPACITY1, VALUE1, CAPACITY1, N1, &solution); // = 19
		printf("ZeroOnePack_Backtrack_2D: 19 -- %d\n", res1);
		res1 = ZeroOnePack_Recursive(MAX_CAPACITY1, VALUE1, CAPACITY1, N1); // = 19
		printf("ZeroOnePack_Recursive: 19 -- %d\n", res1);


		int VALUE2[] = { 6, 12, 24, 10, 20, 12 };
		int CAPACITY2[] = { 1, 2,  4,  2,  4,  3 };
		int MAX_CAPACITY2 = 5;
		int N2 = sizeof(CAPACITY2) / sizeof(CAPACITY2[0]);

		int res2 = ZeroOnePack_Backtrack_1D(MAX_CAPACITY2, VALUE2, CAPACITY2, N2); // = 30
		printf("ZeroOnePack_Backtrack_1D: 30 -- %d\n", res2);
		res2 = ZeroOnePack_Backtrack_2D(MAX_CAPACITY2, VALUE2, CAPACITY2, N2, &solution); // = 30
		printf("ZeroOnePack_Backtrack_2D: 30 -- %d\n", res2);
		res2 = ZeroOnePack_Recursive(MAX_CAPACITY2, VALUE2, CAPACITY2, N2); // = 30
		printf("ZeroOnePack_Recursive: 30 -- %d\n", res2);


		int VALUE3[] = { 20, 6, 20, 4 };
		int CAPACITY3[] = { 4,  3, 4,  2 };
		int MAX_CAPACITY3 = 9;
		int N3 = sizeof(CAPACITY3) / sizeof(CAPACITY3[0]);

		int res3 = ZeroOnePack_Backtrack_1D(MAX_CAPACITY3, VALUE3, CAPACITY3, N3); // = 40
		printf("ZeroOnePack_Backtrack_1D: 40 -- %d\n", res3);
		res3 = ZeroOnePack_Backtrack_2D(MAX_CAPACITY3, VALUE3, CAPACITY3, N3, &solution); // = 40
		printf("ZeroOnePack_Backtrack_2D: 40 -- %d\n", res3);
		res3 = ZeroOnePack_Recursive(MAX_CAPACITY3, VALUE3, CAPACITY3, N3); // = 40
		printf("ZeroOnePack_Recursive: 40 -- %d\n", res3);


		int VALUE4[] = { 6, 3, 5, 4, 6 };
		int CAPACITY4[] = { 2, 2, 6, 5, 4 };
		int MAX_CAPACITY4 = 10;
		int N4 = sizeof(CAPACITY4) / sizeof(CAPACITY4[0]);

		int res4 = ZeroOnePack_Backtrack_1D(MAX_CAPACITY4, VALUE4, CAPACITY4, N4); // = 15
		printf("ZeroOnePack_Backtrack_1D: 15 -- %d\n", res4);
		res4 = ZeroOnePack_Backtrack_2D(MAX_CAPACITY4, VALUE4, CAPACITY4, N4, &solution); // = 15
		printf("ZeroOnePack_Backtrack_2D: 15 -- %d\n", res4);
		res4 = ZeroOnePack_Recursive(MAX_CAPACITY4, VALUE4, CAPACITY4, N4); // = 15
		printf("ZeroOnePack_Recursive: 15 -- %d\n", res4);


		int VALUE5[] = { 60, 100, 120 };
		int CAPACITY5[] = { 10, 20,  30 };
		int MAX_CAPACITY5 = 50;
		int N5 = sizeof(CAPACITY5) / sizeof(CAPACITY5[0]);

		int res5 = ZeroOnePack_Backtrack_1D(MAX_CAPACITY5, VALUE5, CAPACITY5, N5); // = 220
		printf("ZeroOnePack_Backtrack_1D: 220 -- %d\n", res5);
		res5 = ZeroOnePack_Backtrack_2D(MAX_CAPACITY5, VALUE5, CAPACITY5, N5, &solution); // = 220
		printf("ZeroOnePack_Backtrack_2D: 220 -- %d\n", res5);
		res5 = ZeroOnePack_Recursive(MAX_CAPACITY5, VALUE5, CAPACITY5, N5); // = 220 
		printf("ZeroOnePack_Recursive: 220 -- %d\n", res5);
	}

	static void Test_CompletePack()
	{
		int VALUE[] = { 60, 100, 120 };
		int CAPACITY[] = { 10, 20, 30 };
		int MAX_CAPACITY = 50;
		int N = sizeof(CAPACITY) / sizeof(CAPACITY[0]);

		vector<int> solution;

		printf("CompletePack_Backtrack_2D: should be 300, the answer is %d\n", Knapsack::CompletePack_Backtrack_2D(MAX_CAPACITY, VALUE, CAPACITY, N)); // = 300

		printf("CompletePack_Backtrack_ConvertToZeroOne_2D: should be 300, the answer is %d\n", Knapsack::CompletePack_Backtrack_ConvertToZeroOne_2D(MAX_CAPACITY, VALUE, CAPACITY, N, &solution)); // = 300

		printf("CompleteBacktrack_1D: should be 300, the answer is %d\n", Knapsack::CompleteBacktrack_1D(MAX_CAPACITY, VALUE, CAPACITY, N)); // = 300

		int VALUE2[] = { 30, 14, 16, 9 };
		int CAPACITY2[] = { 6, 3, 4, 2 };
		MAX_CAPACITY = 10;
		N = sizeof(CAPACITY2) / sizeof(CAPACITY2[0]);

		printf("CompleteBacktrack_1D: should be 48, the answer is %d\n", Knapsack::CompleteBacktrack_1D(MAX_CAPACITY, VALUE2, CAPACITY2, N)); // = 48
	}

	static void Test_MultiplePack()
	{
		int VALUE[] =    { 60, 100, 120 };
		int CAPACITY[] = { 1,  2,   3 };
		int QUANTITY[] = { 5,  0,   0 };
		int MAX_CAPACITY = 5;
		int N = sizeof(CAPACITY) / sizeof(CAPACITY[0]);

		vector<int> solution;

		printf("MultiplePack_Backtrack_ConvertToZeroOne_2D: should be 300, the answer is %d\n", Knapsack::MultiplePack_Backtrack_ConvertToZeroOne_2D(MAX_CAPACITY, VALUE, CAPACITY, QUANTITY, N, &solution)); // = 300


		int QUANTITY2[] = { 2, 3, 0 };

		solution.clear();

		printf("MultiplePack_Backtrack_ConvertToZeroOne_2D: should be 260, the answer is %d\n", Knapsack::MultiplePack_Backtrack_ConvertToZeroOne_2D(MAX_CAPACITY, VALUE, CAPACITY, QUANTITY2, N, &solution)); // = 300
	}

	static void Test_ExactPack()
	{
		int VALUE[] = { 2, 3, 7 };
		int TOTAL_VALUE = 5;
		int N = sizeof(VALUE) / sizeof(VALUE[0]);

		vector<int> solution;

		printf("ExactPack_Backtrack_2D: should be Yes, the answer is %s\n", Knapsack::ExactPack_Backtrack_2D(TOTAL_VALUE, VALUE, N, &solution) ? "Yes" : "No"); // = Yes

		TOTAL_VALUE = 8;

		printf("ExactPack_Backtrack_2D: should be No, the answer is %s\n", Knapsack::ExactPack_Backtrack_2D(TOTAL_VALUE, VALUE, N, &solution) ? "Yes" : "No"); // = Yes
	}

	static void Test()
	{
		//Test_ZeroOnePack();

		//Test_CompletePack();

		//Test_MultiplePack();

		Test_ExactPack();
	}
};


class MultipleBagKnapsack
{
public:

	int BagCount, ItemCount;

	vector<int> vecCapacity;
	int* capacity;
	int total_capacity;

	vector<int> vecWeight;
	int* weight;
	vector<int> vecWeightSum;
	int* weight_sum;

	vector<int> vecChosenBag;
	int* chosen_bag;
	int max_waste_capacity, current_waste_capacity;

	int Dfs(int item)
	{
		if (current_waste_capacity > max_waste_capacity) {
			return 0;
		}

		if (item < 0) {
			return 1;
		}

		int start;

		if (weight[item] == weight[item + 1]) {
			start = chosen_bag[item + 1];
		}
		else {
			start = 0;
		}

		for (int bag = start; bag < BagCount; bag++)
		{
			if (capacity[bag] >= weight[item])
			{
				capacity[bag] -= weight[item];

				chosen_bag[item] = bag;

				if (capacity[bag] < weight[0]) {
					current_waste_capacity += capacity[bag];
				}

				int remain = Dfs(item - 1);

				if (capacity[bag] < weight[0]) {
					current_waste_capacity -= capacity[bag];
				}

				chosen_bag[item] = -1;

				capacity[bag] += weight[item];

				// We only need one of the solution to put remain rails, once we find it, return true right away. 
				if (remain) {
					return 1;
				}
			}
		}

		return 0;
	}

	int BinarySearch()
	{
		int left, right, mid;

		left = 0;

		for (right = ItemCount - 1; right >= 0 && weight_sum[right] > total_capacity; right--);

		while (left < right)
		{
			mid = (left + right + 1) / 2;

			max_waste_capacity = total_capacity - weight_sum[mid];
			current_waste_capacity = 0;

			if (Dfs(mid)) {
				left = mid;
			}
			else {
				right = mid - 1;
			}
		}

		return left + 1;
	}

	int Compute(int _BagCount, int* pBagCapacity, int _ItemCount, int* pItemWeight)
	{
		BagCount = _BagCount;
		ItemCount = _ItemCount;

		vecCapacity.assign(BagCount, 0);

		total_capacity = 0;
		for (int i = 0; i < BagCount; i++)
		{
			vecCapacity[i] = pBagCapacity[i];
			total_capacity += vecCapacity[i];
		}

		struct greater
		{
			bool operator()(int const &a, int const &b) const { return a > b; }
		};

		sort(vecCapacity.begin(), vecCapacity.end(), greater());
		capacity = &vecCapacity.front();

		vecChosenBag.assign(ItemCount, -1);
		chosen_bag = &vecChosenBag.front();

		vecWeight.assign(ItemCount, 0);
		for (int i = 0; i < ItemCount; i++) {
			vecWeight[i] = pItemWeight[i];
		}

		sort(vecWeight.begin(), vecWeight.end());
		weight = &vecWeight.front();

		vecWeightSum.assign(ItemCount, 0);

		int sum = 0;
		for (int i = 0; i < ItemCount; i++)
		{
			sum += weight[i];

			vecWeightSum[i] = sum;
		}

		weight_sum = &vecWeightSum.front();

		while (ItemCount > 0 && weight[ItemCount - 1] > capacity[0])  ItemCount--;
		while (BagCount > 0 && capacity[BagCount - 1] < weight[0])  BagCount--;

		int count = BinarySearch();

		return count;
	}

public:

	static void Test()
	{
		int BagCapacity[] = { 30, 40, 50, 25 };
		int ItemWeight[] = { 15, 16, 17, 18, 19, 20, 21, 25, 24, 30 };

		MultipleBagKnapsack* pMBK = new MultipleBagKnapsack();

		int max_item_count = pMBK->Compute(sizeof(BagCapacity)/sizeof(BagCapacity[0]), BagCapacity, sizeof(ItemWeight) / sizeof(ItemWeight[0]), ItemWeight);
		bool b = (max_item_count == 7);

		printf("MultipleBagKnapsack, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), max_item_count);

		delete pMBK;
	}
};
