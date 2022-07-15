#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;


///////////////////////////////////////////////////////////////////////////////////////////


class BinarySearchDemo
{
public:

	static int Demo()
	{
		const int N = 100;

		int start = 0;
		int end = N - 1;

		int hit = -1;

		while (start < end)
		{
			int cursor = (start + end) / 2;

			printf("\n BinarySearchDemo, searching on %d ...", cursor);

			int cmp = YourLogicToCheckCondition(cursor);

			if (cmp == 0)
			{
				hit = cursor;

				end = cursor;
			}
			else if (cmp > 0)
			{
				start = cursor;
			}
			else {
				end = cursor;
			}

			if (start == end || (start + 1) == end)
			{
				break;
			}
		}

		printf("\n BinarySearchDemo, hit is %d \n", hit);

		return hit;
	}

	static int YourLogicToCheckCondition(int cursor)
	{
		int target = 53;

		if (cursor < target) {
			return 1;
		}
		if (cursor > target) {
			return -1;
		}
		return 0;
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


//
// Choose K items from N items
//
class IterateCombinationDemo
{
public:

	static void YourLogicToDealWithComb(int* comb, int count)
	{
		//
		// comb now contains chosen indices, replace the following with real handling logic
		//

		if (comb == 0 || count <= 0) {
			return;
		}

		for (int i = 0; i < count; i++) {
			cout << " " << comb[i];
		}
		cout << endl;
	}

	static void Iterate(int K, int N)
	{
		if (K > N) {
			return;
		}

		if (K == 0) {
			YourLogicToDealWithComb(0, 0);
			return;
		}

		vector<unsigned char> vecMask(K, 1);	// K trailing 1's
		vecMask.resize(N, 0);					// N-K leading 0's
		unsigned char* mask = &vecMask.front();

		vector<int> vecComb(K, 0);
		int* comb = &vecComb.front();

		do
		{
			int count = 0;
			for (int i = 0; i < N; i++)		// [0..N-1] integers
			{
				if (mask[i]) {
					comb[count++] = i;
				}
			}

			YourLogicToDealWithComb(comb, K);

		} while (prev_permutation(vecMask.begin(), vecMask.end()));
	}

	static void Iterate(int N)
	{
		for (int k = 0; k <= N; k++) {
			Iterate(k, N);
		}
	}

	static void IterateBruteForce(int N)
	{
		vector<int> vecMask(N, 0);
		int* mask = &vecMask.front();

		vector<int> vecComb(N, 0);
		int* comb = &vecComb.front();

		for (int i = 0; ; i++)
		{
			int k = 1;
			for (int j = 0; j < N; j++)
			{
				mask[j] += k;

				k = mask[j] / 2;

				mask[j] %= 2;
			}

			if (k) {
				break;
			}

			int count = 0;
			for (int j = 0; j < N; j++)		// [0..N-1] integers
			{
				if (mask[j]) {
					comb[count++] = j;
				}
			}

			YourLogicToDealWithComb(comb, count);
		}
	}

	static void BeginIterateRecursive(int K, int N)
	{
		vector<int> vecComb;
		vecComb.reserve(K);

		IterateRecursive(0, K, N, vecComb);
	}

	static void BeginIterateRecursive(int N)
	{
		YourLogicToDealWithComb(0, 0);

		for (int k = 1; k <= N; k++)
		{
			vector<int> vecComb;
			vecComb.reserve(k);

			IterateRecursive(0, k, N, vecComb);
		}
	}

	static void IterateRecursive(int offset, int K, int N, vector<int>& vecComb)
	{
		if (K == 0) {
			YourLogicToDealWithComb(&vecComb.front(), (int)vecComb.size());
			return;
		}

		for (int i = offset; i <= N - K; i++)
		{
			vecComb.push_back(i);

			IterateRecursive(i + 1, K - 1, N, vecComb);

			vecComb.pop_back();
		}
	}

	static void Demo()
	{
		Iterate(4);

		IterateBruteForce(4);

		BeginIterateRecursive(4);
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


static void SampleTester()
{
	BinarySearchDemo::Demo();

	IterateCombinationDemo::Demo();
}


///////////////////////////////////////////////////////////////////////////////////////////

