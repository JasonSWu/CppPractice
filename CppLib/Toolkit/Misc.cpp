#include <algorithm>
#include <cstdio>
#include <vector>
#include <queue>
#include <map>
#include <limits>
#include <stack>

using namespace std;




class Misc
{
public:

	static int MaxRectangularAreaOfHistogramBars(vector<int>& vecHeight)
	{
		int count = (int)vecHeight.size();
		int* height = &vecHeight.front();

		int ma = 0;

		stack<int> st;

		for (int x = 0; x < count; ++x)
		{
			if (!st.empty() && height[x] < height[st.top()])
			{
				while (!st.empty() && height[st.top()] > height[x])
				{
					int h = height[st.top()];
					st.pop();

					int area = h * (st.empty() ? x : x - st.top() - 1);

					ma = max(ma, area);
				}
			}

			st.push(x);
		}

		while (!st.empty())
		{
			int h = height[st.top()];
			st.pop();

			int area = h * (st.empty() ? count : count - st.top() - 1);

			ma = max(ma, area);
		}

		return ma;
	}

public:

	static void TestMaxRectangularAreaOfHistogramBars()
	{
		vector<int> vecHeight = { 2, 1, 5, 6, 2, 3 };

		int ma = MaxRectangularAreaOfHistogramBars(vecHeight);

		printf("MaxRectangularAreaOfHistogramBars: %s, Max area is %d \n", (ma == 10) ? "Correct" : "Wrong", 10);
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


static void MiscTester()
{
	Misc::TestMaxRectangularAreaOfHistogramBars();
}


///////////////////////////////////////////////////////////////////////////////////////////


