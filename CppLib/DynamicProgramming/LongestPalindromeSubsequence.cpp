#include <cstdio>
#include <vector>
#include <algorithm>

using namespace std;




///////////////////////////////////////////////////////////////////////////////////////////


//
// Reference: http://www.geeksforgeeks.org/?p=19155
//
class LongestPalindromeSubsequence
{
public:

	static void LPS(const char* str, int& pos, int& len)
	{
		int length = strlen(str);

		//
		// IsPal[i][j] indicates if substring str[i..j] is a palindrome
		//
		vector<vector<unsigned char>> IsPal(length, vector<unsigned char>(length, 0));

		for (int i = 0; i < length; i++) {
			IsPal[i][i] = true;
		}

		pos = 0;
		len = 1;

		//
		// check for sub-string of length 2
		//
		for (int i = 0; i < length - 1; ++i)
		{
			if (str[i] == str[i + 1])
			{
				IsPal[i][i + 1] = true;
				pos = i;
				len = 2;
			}
		}

		//
		// check for sub-string of length greater than 2
		//
		// Note that the lower diagonal values of table are useless and not filled in the process. 
		// The values are filled in a manner similar to Matrix Chain Multiplication DP solution (See
		// http://www.geeksforgeeks.org/archives/15553). k is length of substring
		//
		for (int k = 3; k <= length; k++)
		{
			for (int i = 0; i<length - k + 1; i++)
			{
				int j = i + k - 1;

				if (IsPal[i + 1][j - 1] && str[i] == str[j])
				{
					IsPal[i][j] = true;

					if (k > len)
					{
						pos = i;
						len = k;
					}
				}
			}
		}
	}


	static void Test()
	{
		string str;
		int pos, len;

		str = "AXYA";
		LPS(str.c_str(), pos, len);

		bool b = (len == 1);
		printf("LongestPalindromeSubsequence::LPS, %s, the answer is %d %s \n", (b ? "Correct" : "Wrong"), len, str.substr(pos, len).c_str());


		str = "forgeeksskeegfor";
		LPS(str.c_str(), pos, len);

		b = (len == 10);
		printf("LongestPalindromeSubsequence::LPS, %s, the answer is %d %s \n", (b ? "Correct" : "Wrong"), len, str.substr(pos, len).c_str());
	}

};

