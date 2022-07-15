#include <cstdio>
#include <vector>
#include <algorithm>

using namespace std;




///////////////////////////////////////////////////////////////////////////////////////////


//
// Reference: http://algorithms.tutorialhorizon.com/dynamic-programming-longest-increasing-subsequence/
// Reference: http://www.geeksforgeeks.org/longest-monotonically-increasing-subsequence-size-n-log-n/
//
class LongestSubsequence
{
public:

	//
	// O(n^2)
	//
	static int LIS(int* A, int size, vector<int>* pSolution = 0)
	{
		//
		// Length of lis with A[i] as its last element.
		//
		// lis[i] = 1 + max of { lis[j] where 1<j<i, if A[i] > A[j] }
		//				 = 1 if no such j exists.
		// 
		vector<int> lis(size, 1);
		
		for (int i = 0; i < size; i++) // i is the next element to be considered
		{
			int maxLen = -1;

			for (int j = 0; j < i; j++) // check against all A[j] before A[i] 
			{
				if (A[i] > A[j]) {	// A[i] should be considered since its value increases
					maxLen = max(maxLen, lis[j] + 1); // lis[i] should build on the longest lis[j] before i
				}
			}

			if (maxLen > 0) {
				lis[i] = maxLen;
			}
		}

		int largest = max_element(lis.begin(), lis.end()) - lis.begin();

		if (pSolution != 0)
		{
			pSolution->push_back(largest);

			int q = lis[largest] - 1;

			for (int k = largest - 1; k >= 0; k--)
			{
				if (lis[k] == q)
				{
					pSolution->push_back(k);
					q--;
				}
			}

			reverse(pSolution->begin(), pSolution->end());
		}

		return lis[largest];
	}

	static int LIS_Recursive(int* A, int size, int startPosition, int subSeqLengthSoFar, int largestNumberSoFar)
	{
		int longestSubSeqLength = subSeqLengthSoFar;

		for (int i = startPosition; i < size; i++)
		{
			if (A[i] > largestNumberSoFar)
			{
				int len = LIS_Recursive(A, size, i, subSeqLengthSoFar + 1, A[i]);

				if (len > longestSubSeqLength) {
					longestSubSeqLength = len;
				}
			}
		}

		return longestSubSeqLength;
	}

	//
	// O(n^2)
	//
	// This enhanced version also counts the number of each lds[i] with repeated sequences excluded.
	// Reference: http://blog.csdn.net/dingyaguang117/article/details/5836918
	//
	static int LDS(vector<int>& A, long long& num)
	{
		A.push_back(0);  //for simplicity

		int size = (int)A.size();
		
		//
		// Length of LDS with A[i] as its last element.
		//
		// lds[i] = 1 + max of { lds[j] where 1<j<i, if A[i] < A[j] }
		//				 = 1 if no such j exists.
		// 
		vector<int> lds(size, 1);

		vector<long long> count(size, 1);

		vector<int> next(size, 0);

		for (int i = 0; i < size; i++) {
			for (int j = i + 1; j < size; j++)
			{
				if (A[i] == A[j])
				{
					next[i] = j;
					break;
				}
			}
		}

		for (int i = 0; i < size; i++) 
		{
			for (int j = 0; j < i; j++) 
			{
				if ( A[i] < A[j]  && (next[j] == 0 || next[j] > i) )
				{
					if(lds[j] + 1 > lds[i])
					{
						lds[i] = lds[j] + 1;
						count[i] = count[j];
					}
					else if(lds[j] + 1 == lds[i])
					{
						count[i] += count[j];
					}
				}
			}
		}

		A.pop_back();

		num = count[size - 1];

		return lds[size-1]-1;
	}

	static void Test()
	{
		Test_01();
		Test_02();
	}

	static void Test_01()
	{
		int A[] = { 1, 12, 7, 0, 23, 11, 52, 31, 61, 69, 70, 2 };

		int len = LIS_Recursive(A, sizeof(A) / sizeof(A[0]), 0, 0, -1);
		bool b = (len == 7);

		printf("LongestSubsequence::LIS_Recursive, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), len);
		

		vector<int> solution; // = {0, 2, 5, 7, 8, 9, 10}  => { 1, 7, 11, 31, 61, 69, 70 }
		len = LIS(A, sizeof(A) / sizeof(A[0]), &solution);
		b = (len == 7);

		printf("LongestSubsequence::LIS, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), len);
	}

	static void Test_02()
	{
		vector<int> A = { 68, 69, 54, 64, 68, 64, 70, 67, 78, 62, 98, 87 };

		long long count;
		int len = LDS(A, count);
		bool b = (len == 4) && (count == 2);;

		printf("LongestSubsequence::LDS, %s, the length is %d, and the count is %lld \n", (b ? "Correct" : "Wrong"), len, count);


		vector<int> A2 = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

		len = LDS(A2, count);
		b = (len == 1) && (count == 10);

		printf("LongestSubsequence::LDS, %s, the length is %d, and the count is %lld \n", (b ? "Correct" : "Wrong"), len, count);
	}
};

