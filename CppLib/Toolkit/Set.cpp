
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>

using namespace std;


///////////////////////////////////////////////////////////////////////////////////////////


class UnionFind
{
public:

	vector<int> Parents;
	int* pParents;
	vector<int> ElementCounts;
	int* pElementCounts;

	int TotalSetCount;

public:

	UnionFind(int setCount = 0)
	{
		Initialize(setCount);
	}

	void Initialize(int setCount)
	{
		pParents = 0;
		pElementCounts = 0;
		TotalSetCount = setCount;

		if (TotalSetCount > 0)
		{
			ElementCounts.assign(TotalSetCount, 1);
			pElementCounts = &ElementCounts.front();

			Parents.assign(TotalSetCount, 0);
			pParents = &Parents.front();

			for (int i = 0; i < TotalSetCount; i++) {
				pParents[i] = i;
			}
		}
	}

	int FindSet(int e) const
	{
		while (e != pParents[e])
		{
			e = pParents[e];
		}

		return e;
	}

	int FindSetSize(int e) const
	{
		return pElementCounts[FindSet(e)];
	}

	bool AreInSameSet(int e1, int e2) const
	{
		return FindSet(e1) == FindSet(e2);
	}

	void Union(int e1, int e2)
	{
		int s1 = FindSet(e1);
		int s2 = FindSet(e2);

		if (s1 == s2)
		{
			return;
		}

		//
		// Weighted Quick Union. Make smaller root point to larger one.
		//
		if (pElementCounts[s1] < pElementCounts[s2])
		{
			pParents[s1] = s2;
			pElementCounts[s2] += pElementCounts[s1];
		}
		else
		{
			pParents[s2] = s1;
			pElementCounts[s1] += pElementCounts[s2];
		}

		TotalSetCount--;
	}

public:

	static void Test()
	{
		Test("input_union_find_tiny.txt", 2);
		Test("input_union_find_medium.txt", 3);
		//Test("input_union_find_large.txt", 6);
	}

	static void Test(const char* dataFileName, int expectedSetCount)
	{
		ifstream fin(dataFileName);

		int N;
		fin >> N;

		UnionFind UF(N);

		int e1, e2;
		while (!fin.eof())
		{
			e1 = e2 = -1;

			fin >> e1 >> e2;

			if (e1 < 0 || e2 < 0) {
				break;
			}

			if (UF.AreInSameSet(e1, e2)) {
				continue;
			}

			UF.Union(e1, e2);
		}

		fin.close();

		printf("UnionFind: %s, Total Set Count is %d \n", (UF.TotalSetCount == expectedSetCount) ? "Correct" : "Wrong", UF.TotalSetCount);
	}
};


///////////////////////////////////////////////////////////////////////////////////////////



static void SetTester()
{
	UnionFind::Test();
}


///////////////////////////////////////////////////////////////////////////////////////////
