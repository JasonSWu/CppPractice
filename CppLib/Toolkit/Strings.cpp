#include <cstdio>
#include <string>
#include <vector>
#include <set>
#include <algorithm>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////


//
// Find patterns (sub-strings) in a string.
// KMP algorithm complexity is O(n) in worst case with n being the text length
//
class KnuthMorrisPrattAlgorithm
{
public:

	int PatternLength;
	int TextLength;

	//
	// LPS[i] is the length of the longest proper prefix of Pattern[0..i] which is also a suffix of Pattern[0..i]
	//
	vector<int> LPS;

	vector<int> Hits;

public:

	KnuthMorrisPrattAlgorithm()
	{
		Clear();
	}

	void Clear()
	{
		this->PatternLength = 0;
		this->TextLength = 0;

		this->LPS.clear();
		this->Hits.clear();
	}

	//
	// Return value indicates the index of the first match
	//
	int Search(const char* pPattern, const char* pText, bool bFindFirstMatchOnly = false)
	{
		Clear();

		if (pPattern == 0 || pText == 0)
		{
			return -1;
		}

		int M = strlen(pPattern);
		int N = strlen(pText);

		if (M > N)
		{
			return -1;
		}

		this->PatternLength = M;
		this->TextLength = N;

		ComputeLPS(pPattern);

		int iIndexOfFirstMatch = -1;

		int j = 0;  // index for pattern
		int i = 0;  // index for text
		int* pLPS = &this->LPS.front();

		while (i < this->TextLength)
		{
			if (pPattern[j] == pText[i])
			{
				j++;
				i++;
			}

			if (j == this->PatternLength)
			{
				//printf("Found pattern at index %d \n", i - j);

				this->Hits.push_back(i - j);

				if (iIndexOfFirstMatch < 0)
				{
					iIndexOfFirstMatch = i - j;
				}

				j = pLPS[j - 1];

				if (bFindFirstMatchOnly)
				{
					return iIndexOfFirstMatch;
				}
			}
			else if (i < this->TextLength && pPattern[j] != pText[i])  // mismatch after j matches
			{
				// Do not match LPS[0..LPS[j-1]] characters, they will match anyway
				if (j != 0)
				{
					j = pLPS[j - 1];
				}
				else
				{
					i++;
				}
			}
		}

		return iIndexOfFirstMatch;
	}

	void ComputeLPS(const char* pPattern)
	{
		if (pPattern == 0 || pPattern[0] == 0)
		{
			return;
		}

		if ( this->PatternLength == 0)
		{
			this->PatternLength = strlen(pPattern);
		}

		this->LPS.assign(this->PatternLength, 0);
		this->LPS[0] = 0; // lps[0] is always 0

		int len = 0;  // length of the previous longest prefix suffix
		int i = 1;
		int* pLPS = &this->LPS.front();

		// the loop calculates LPS[i] for i = 1 to M-1
		while (i < this->PatternLength)
		{
			if (pPattern[i] == pPattern[len])
			{
				len++;
				pLPS[i] = len;
				i++;
			}
			else // (pPattern[i] != pPattern[len])
			{
				if (len != 0)
				{
					// This is tricky. Consider the example 
					// AAACAAAA and i = 7.
					len = pLPS[len - 1];

					// Also, note that we do not increment i here
				}
				else // if (len == 0)
				{
					pLPS[i] = 0;
					i++;
				}
			}
		}
	}

public:

	static void Test()
	{
		KnuthMorrisPrattAlgorithm KMP;

		KMP.Search("AABA", "AABAACAADAABAAABAA"); // should return position 0, 9, 13
	}
};


///////////////////////////////////////////////////////////////////////////////////////////

//
// Find patterns (sub-strings) in a string.
// Rabin-Karp algorithm complexity is O(n+m) with n being the text length and m being the pattern length,
// but the worse case is O(mn).
//
class RabinKarpAlgorithm
{
public:

	vector<int> Hits;

public:

	RabinKarpAlgorithm()
	{
		Clear();
	}

	void Clear()
	{
		this->Hits.clear();
	}

	//
	// q is a prime number
	// 
	void Search(const char* pPattern, const char* pText, int mod = 101)
	{
		Clear();

		const int base = 256;

		int M = strlen(pPattern);
		int N = strlen(pText);

		int i, j;
		int hash_pattern = 0; // hash value for pattern
		int hash_source = 0; // hash value for txt
		int hipow = 1;

		// The value of hipow would be "pow(base, M-1)%mod"
		for (i = 0; i < M - 1; i++) {
			hipow = (hipow*base) % mod;
		}

		// Calculate the hash value of pattern and first window of text
		for (i = 0; i < M; i++)
		{
			hash_pattern = (base*hash_pattern + pPattern[i]) % mod;
			hash_source = (base*hash_source + pText[i]) % mod;
		}

		// Slide the pattern over text one by one
		for (i = 0; i <= N - M; i++)
		{

			// Check the hash values of current window of text
			// and pattern. If the hash values match then only
			// check for characters on by one
			if (hash_pattern == hash_source)
			{
				// Check for characters one by one 
				for (j = 0; j < M; j++)
				{
					if (pText[i + j] != pPattern[j])
					{
						break;
					}
				}

				// if hash_pattern == hash_source and pPattern[0...M-1] = pText[i, i+1, ...i+M-1]
				if (j == M)
				{
					//printf("Pattern found at index %base \n", i);

					this->Hits.push_back(i);
				}
			}

			// Calculate hash value for next window of text: Remove
			// leading digit, add trailing digit
			if (i < N - M)
			{
				hash_source = (base*(hash_source - pText[i] * hipow) + pText[i + M]) % mod;

				// We might get negative value of hash_source, converting it
				// to positive
				if (hash_source < 0)
				{
					hash_source = (hash_source + mod);
				}
			}
		}
	}


public:

	static void Test()
	{
		RabinKarpAlgorithm BK;

		BK.Search("AABA", "AABAACAADAABAAABAA"); // should return position 0, 9, 13
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class SuffixArray
{
public:

	char* T;			// the input string, up to 100K characters
	int N;              // the length of input string
	vector<char> Text;

	vector<char> Text2;
	int N2;             // the length of input string

	const char* P;		// the pattern string (for string matching)
	int M;              // the length of pattern string

	vector<int> Ranks, TempRanks;    // rank array and temporary rank array

	vector<int> Suffix, TempSuffix;  // suffix array and temporary suffix array

	vector<int> Radix;  // for counting/radix sort

	vector<int> Phi;    // for computing longest common prefix
	vector<int> PLCP;

	vector<int> LCP;	// LCP[i] stores the LCP between previous suffix T+Suffix[i-1] and current suffix T+Suffix[i]

public:

	SuffixArray(const char* pText = 0, char endingChar = '$')
	{
		this->N2 = 0;

		this->P = 0;
		this->M = 0;

		SetText(pText, endingChar);
	}

	void SetText(const char* pText, char endingChar = '$')
	{
		this->N = (pText == 0) ? 0 : strlen(pText);
		this->T = 0;

		if (this->N == 0)
		{
			return;
		}

		//const int MAX_N = (100 * 1024);  
		//this->Text.assign(MAX_N, 0);

		this->Text.assign(this->N + 2, 0);
		this->T = &this->Text.front();

		strcpy(this->T, pText);

		//if (endingChar != 0)
		{
			this->T[this->N++] = endingChar;
		}
		
		//
		// if '\n' is read, uncomment the next line
		//
		// this->T[this->N-1] = endingChar; this->T[this->N] = 0;
		//

		this->Ranks.assign(this->N, 0);
		this->TempRanks.assign(this->N, 0);

		this->Suffix.assign(this->N, 0);
		this->TempSuffix.assign(this->N, 0);

		this->Radix.assign(max(300, this->N), 0);  // up to 255 ASCII chars or length of this->N

		this->Phi.assign(this->N, 0);
		this->PLCP.assign(this->N, 0);

		this->LCP.assign(this->N, 0);

		BuildSuffixArray();

		ComputeLongestCommonPrefix();
	}

	void SetText2(const char* pText2, char endingChar = '#')
	{
		if (pText2 == 0 || pText2[0] == 0)
		{
			return;
		}

		this->N2 = strlen(pText2);

		this->Text2.assign(this->N2 + 1, 0);

		char* T2 = &this->Text2.front();
		strcpy(T2, pText2);
		//
		// if '\n' is read, uncomment the next line
		//
		// T2[this->N2-1] = 0; this->N2--;
		//

		vector<char> str( this->N + this->N2 + 1);
		char* p = &str.front();

		strcpy(p, this->T);  
		strcat(p, T2);

		SetText(p, endingChar);
	}

	// O(n^2 log n)
	void BuildSuffixArrayNaive()                // cannot go beyond 1000 characters
	{
		for (int i = 0; i < this->N; i++)		// initial SA: {0, 1, 2, ..., this->N-1}
		{
			this->Suffix[i] = i;
		}

		struct StringLessThanCompare
		{
			const char* T;

			StringLessThanCompare(const char* t)
			{
				T = t;
			}

			bool operator()(int a, int b)
			{
				return strcmp(T + a, T + b) < 0;
			}
		};

		StringLessThanCompare SLTC(this->T);

		sort(this->Suffix.begin(), this->Suffix.end(), SLTC); // sort: O(n log n) * compare: O(n) = O(n^2 log n)
	}

	void CountingSort(int k)  // O(this->N)
	{                                          
		int i, sum, maxi = (int)this->Radix.size();   

		fill(this->Radix.begin(), this->Radix.end(), 0);

		for (i = 0; i < this->N; i++)			// count the frequency of each integer rank
		{
			this->Radix[i + k < this->N ? this->Ranks[i + k] : 0]++;
		}
		
		for (i = sum = 0; i < maxi; i++) 
		{
			int t = this->Radix[i]; 
			this->Radix[i] = sum; 
			sum += t;
		}

		for (i = 0; i < this->N; i++)          // shuffle the suffix array if necessary
		{
			this->TempSuffix[this->Radix[this->Suffix[i] + k < this->N ? this->Ranks[this->Suffix[i] + k] : 0]++] = this->Suffix[i];
		}
		
		for (i = 0; i < this->N; i++)         // update the suffix array SA
		{
			this->Suffix[i] = this->TempSuffix[i];
		}
	}

	// O(n log n)
	void BuildSuffixArray()          // this version can go up to 100000 characters
	{	
		int i, k, r;

		for (i = 0; i < this->N; i++) // initial rankings
		{
			this->Ranks[i] = this->T[i];                 
		}

		for (i = 0; i < this->N; i++) // initial SA: {0, 1, 2, ..., this->N-1}
		{
			this->Suffix[i] = i;
		}

		for (k = 1; k < this->N; k <<= 1)        // repeat sorting process log this->N times
		{
			CountingSort(k);	// actually radix sort: sort based on the second item
			CountingSort(0);    // then (stable) sort based on the first item

			this->TempRanks[this->Suffix[0]] = r = 0;      // re-ranking; start from rank r = 0

			for (i = 1; i < this->N; i++)        // compare adjacent suffixes
			{
				this->TempRanks[this->Suffix[i]] = // if same pair => same rank r; otherwise, increase r
					(this->Ranks[this->Suffix[i]] == this->Ranks[this->Suffix[i - 1]] && this->Ranks[this->Suffix[i] + k] == this->Ranks[this->Suffix[i - 1] + k]) ? r : ++r;
			}
			
			for (i = 0; i < this->N; i++)                   // update the rank array RA
			{
				this->Ranks[i] = this->TempRanks[i];
			}
			
			if (this->Ranks[this->Suffix[this->N - 1]] == this->N - 1)// nice optimization trick
			{
				break;
			}
		}
	}

	void ComputeLongestCommonPrefixNaive() 
	{
		this->LCP[0] = 0;                                 // default value

		for (int i = 1; i < this->N; i++)                 // compute LCP by definition
		{
			int L = 0;                                       // always reset L to 0
			while (this->T[this->Suffix[i] + L] == this->T[this->Suffix[i - 1] + L])      // same L-th char, L++
			{
				L++;
			}
			
			this->LCP[i] = L;
		}
	}

	// O(n)
	void ComputeLongestCommonPrefix()
	{
		int i, L;
		this->Phi[this->Suffix[0]] = -1;                         // default value

		for (i = 1; i < this->N; i++)                            // compute Phi in O(this->N)
		{
			this->Phi[this->Suffix[i]] = this->Suffix[i - 1];    // remember which suffix is behind this suffix
		}

		for (i = L = 0; i < this->N; i++)              // compute Permuted LCP in O(this->N)
		{
			if (this->Phi[i] == -1)
			{
				this->PLCP[i] = 0;
				continue;
			}   // special case

			while (this->T[i + L] == this->T[this->Phi[i] + L])       // L increased max this->N times
			{
				L++;
			}

			this->PLCP[i] = L;
			L = max(L - 1, 0);                             // L decreased max this->N times
		}

		for (i = 0; i < this->N; i++)                      // compute LCP in O(this->N)
		{
			this->LCP[i] = this->PLCP[this->Suffix[i]];    // put the permuted LCP to the correct position
		}
	}

	//
	// string matching in O(m log n)
	// returns a pair ( lower and upperbound of the Suffix)
	//
	pair<int, int> FindPattern(const char* pPattern, vector<string>* pMatches = 0)        
	{
		if (pMatches != 0)
		{
			pMatches->clear();
		}

		this->P = pPattern;
		this->M = (this->P == 0) ? 0 : strlen(this->P);

		if (this->M == 0)
		{
			return pair<int, int>(-1, -1);    // if not found
		}

		//
		// if '\n' is read, uncomment the following lines
		//
		// this->P[--this->M] = 0; 
		//
		// if (this->M == 0)
		// {
		//	  return pair<int, int>(-1, -1);    // if not found
		// }
		//

		vector<char> str;
		str.assign(this->N + 1, 0);
		char* p = &str.front();

		int lo = 0, hi = this->N - 1, mid = lo;              // valid matching = [0..n]

		while (lo < hi)                                      // find lower bound
		{	
			mid = (lo + hi) / 2;                             // this is round down
			
			int res = strncmp(this->T + this->Suffix[mid], this->P, this->M);  // try to find P in suffix 'mid'
			
			if (res >= 0) hi = mid;        // prune upper half (notice the >= sign)
			else          lo = mid + 1;    // prune lower half including mid
		}                                  // observe `=' in "res >= 0" above

		if (strncmp(this->T + this->Suffix[lo], this->P, this->M) != 0)
		{
			return pair<int, int>(-1, -1);    // if not found
		}
		
		pair<int, int> ans; 
		ans.first = lo;
		
		lo = 0; hi = this->N - 1; mid = lo;
		
		while (lo < hi)            // if lower bound is found, find upper bound
		{
			mid = (lo + hi) / 2;

			int res = strncmp(this->T + this->Suffix[mid], this->P, this->M);

			if (res > 0) hi = mid;                // prune upper half
			else         lo = mid + 1;            // prune lower half including mid
		}                           // (notice the selected branch when res == 0)
		
		if (strncmp(this->T + this->Suffix[hi], this->P, this->M) != 0)     // special case
		{
			hi--;
		}

		ans.second = hi;
		
		if (pMatches != 0 && ans.first >= 0 && ans.second >= 0)
		{
			for (int i = ans.first; i <= ans.second; i++)
			{
				strcpy(p, this->T + this->Suffix[i]);
				p[strlen(p) - 1] = 0;

				pMatches->push_back(p);
			}
		}

		return ans;

	} 

	//
	// returns the LRS index and its length 
	// O(n)
	//
	pair<int, int> FindLongestRepeatingSuffix() const    
	{
		int i, idx = 0, maxLCP = -1;
		
		for (i = 1; i < this->N; i++)      // O(n), start from i = 1
		{
			if (this->LCP[i] > maxLCP)
			{
				maxLCP = this->LCP[i], idx = i;
			}
		}

		return pair<int, int>(idx, maxLCP);
	}

	//
	// returns the LRS string
	//
	void GetLRS(string& lrs, int* pIndex = 0, int* pLength = 0) const
	{
		pair<int, int> ans = FindLongestRepeatingSuffix();

		if (ans.first < 0 || ans.second <= 0)
		{
			lrs = "";
			return;
		}

		vector<char> vec;
		vec.assign(ans.second + 1, 0);

		char* p = &vec.front();

		strncpy(p, this->T + this->Suffix[ans.first], ans.second);

		lrs = p;

		if (pIndex != 0)
		{
			*pIndex = ans.first;
		}

		if (pLength != 0)
		{
			*pLength = ans.second;
		}
	}

	//
	// Belong to first or second text ( for computing Longest Common Substring )
	//
	int Owner(int idx) const
	{ 
		return (idx < this->N - this->N2 - 1) ? 1 : 2; 
	}

	//
	// returns a pair (the index and the length)
	//
	pair<int, int> ComputeLongestCommonSubstring()                  
	{	
		int i, idx = 0, maxLCP = -1;
		
		for (i = 1; i < this->N; i++)     // O(n), start from i = 1
		{
			if (Owner(this->Suffix[i]) != Owner(this->Suffix[i - 1]) && this->LCP[i] > maxLCP)
			{
				maxLCP = this->LCP[i], idx = i;
			}
		}
		
		return pair<int, int>(idx, maxLCP);
	}

	void GetLCS(string& lcs, int* pIndex = 0, int* pLength = 0)
	{
		pair<int, int> ans = ComputeLongestCommonSubstring();

		if (ans.first < 0 || ans.second <= 0)
		{
			lcs = "";
			return;
		}

		vector<char> vec;
		vec.assign(ans.second + 1, 0);

		char* p = &vec.front();

		strncpy(p, this->T + this->Suffix[ans.first], ans.second);

		lcs = p;

		if (pIndex != 0)
		{
			*pIndex = ans.first;
		}

		if (pLength != 0)
		{
			*pLength = ans.second;
		}
	}

	static void GetLongestPalindrome(const char* pText, string& lp, char endingChar = '#')
	{
		lp = "";

		if (pText == 0 || pText[0] == 0)
		{
			return;
		}

		int originalLen = strlen(pText);

		char ending[2];
		ending[0] = endingChar;
		ending[1] = 0;

		vector<char> text;
		text.assign(originalLen*2 + 1 + 1, 0);

		string rev(pText);

		reverse(rev.begin(), rev.end());

		char* p = &text.front();

		strcpy(p, pText);
		strcat(p, ending);
		strcat(p, rev.c_str());

		SuffixArray SA(p, 0);

		int longestLength = 0;
		int pos = 0;

		for (int i = 1; i<SA.N; i++)
		{
			if ( SA.LCP[i] > longestLength)
			{
				if ((SA.Suffix[i - 1]<originalLen && SA.Suffix[i]>originalLen) || (SA.Suffix[i]<originalLen && SA.Suffix[i - 1]>originalLen))
				{
					longestLength = SA.LCP[i];

					pos = SA.Suffix[i];
				}
			}
		}

		vector<char> palindrome;
		palindrome.assign(longestLength + 1, 0);

		strncpy(&palindrome.front(), &p[pos], longestLength - 1);

		lp = &palindrome.front();
	}

public:

	static void Test()
	{
		SuffixArray SA("GATAGACA");

		printf("\nThe Suffix Array of string T = '%s' is shown below O(n log n) version:\n", SA.T);
		printf("i\tSA[i]\tSuffix\n");
		for (int i = 0; i < SA.N; i++) 
		{
			printf("%2d\t%2d\t%s\n", i, SA.Suffix[i], SA.T + SA.Suffix[i]);
		}
		
		string LRS;
		SA.GetLRS(LRS);
		printf("\nThe Longest Repeating String (LRS) is '%s' with length %d\n\n", LRS.c_str(), LRS.length());

		//printf("\nNow, enter a string P below, we will try to find P in T:\n");
		vector<string> matches;
		pair<int, int> pos = SA.FindPattern("A", &matches);
		if (matches.size() > 0)
		{
			printf("'%s' is found at Suffix[%d..%d]:\n", SA.P, pos.first, pos.second);
			for (int i = 0; i < (int)matches.size(); i++)
			{
				printf("%d  %s\n", SA.Suffix[pos.first+i], matches[i].c_str());
			}
		}
		else printf("%s is not found\n", SA.P);


		SA.SetText2("CATA");

		printf("\nThe LCP information of 'T+P' = '%s':\n", SA.T);
		printf("i\tSA[i]\tLCP[i]\tOwner\tSuffix\n");
		for (int i = 0; i < SA.N; i++)
		{
			printf("%2d\t%2d\t%2d\t%2d\t%s\n", i, SA.Suffix[i], SA.LCP[i], SA.Owner(SA.Suffix[i]), SA.T + SA.Suffix[i]);
		}

		string LCS;
		SA.GetLCS(LCS);

		printf("\nThe Longest Common Substring (LCS) is '%s' with length = %d\n", LCS.c_str(), LCS.length());

		string palindrome;
		GetLongestPalindrome("banana", palindrome);

		printf("\nThe longest palindrome of 'banana' is '%s'\n", palindrome.c_str());
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class Trie
{
public:

	struct TrieComparator {
		bool operator() (Trie* p1, Trie* p2) const
		{
			return p1->c < p2->c;
		}
	};


	struct TrieFindPredicate {
		char c;

		TrieFindPredicate(char ch) {
			c = ch;
		}

		bool operator() (Trie* p) const {
			return p->c == c;
		}
	};

public:

	char c;
	bool isKeyNode;
	set<Trie*, TrieComparator> children;

public:

	Trie() {
		clear();
	}

	~Trie() {
		clear();
	}

	void clear()
	{
		c = 0;
		isKeyNode = false;

		for (auto it = children.begin(); it != children.end(); it++) {
			delete *it;
		}

		children.clear();
	}

	void AddKey(const char* pszKey)
	{
		Trie* pTrie = this;
		auto pKey = pszKey;

		while (*pKey)
		{
			auto it = find_if(pTrie->children.begin(), pTrie->children.end(), TrieFindPredicate(*pKey));

			if (it == pTrie->children.end())
			{
				Trie* pNewTrie = new Trie();
				pNewTrie->c = *pKey;

				it = pTrie->children.insert(pNewTrie).first;
			}

			pTrie = *it;
			pKey++;
		}

		pTrie->isKeyNode = true;
	}

	bool HasKey(const char* pszKey)
	{
		Trie* pTrie = this;
		auto pKey = pszKey;

		while (*pKey)
		{
			auto it = find_if(pTrie->children.begin(), pTrie->children.end(), TrieFindPredicate(*pKey));

			if (it == pTrie->children.end()) {
				return false;
			}

			pTrie = *it;
			pKey++;
		}

		return pTrie->isKeyNode;
	}

	static void Test()
	{
		Trie root;

		root.AddKey("the");
		root.AddKey("there");
		root.AddKey("bye");

		printf("\"the\" is%s a key \n", root.HasKey("the") ? "" : " not");
		printf("\"th\" is%s a key \n", root.HasKey("th") ? "" : " not");
		printf("\"byebye\" is%s a key \n", root.HasKey("byebye") ? "" : " not");
	}
};


///////////////////////////////////////////////////////////////////////////////////////////

class StringUtility
{
public:

	// This function prints the longest palindrome substring (LPS)
	// of str[]. It also returns the length of the longest palindrome
	void LongestPalindromeSubsequence(const char *str, int& pos, int& len)
	{
		pos = 0;
		len = 1;  // The result (length of LPS)

		int length = strlen(str);

		int low, high;

		// One by one consider every character as center point of 
		// even and length palindromes
		for (int i = 1; i < length; ++i)
		{
			// Find the longest even length palindrome with center points
			// as i-1 and i.  
			low = i - 1;
			high = i;
			while (low >= 0 && high < length && str[low] == str[high])
			{
				if (high - low + 1 > len)
				{
					pos = low;
					len = high - low + 1;
				}
				--low;
				++high;
			}

			// Find the longest odd length palindrome with center 
			// point as i
			low = i - 1;
			high = i + 1;
			while (low >= 0 && high < length && str[low] == str[high])
			{
				if (high - low + 1 > len)
				{
					pos = low;
					len = high - low + 1;
				}
				--low;
				++high;
			}
		}
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


static void StringTester()
{
	//KnuthMorrisPrattAlgorithm::Test();

	//RabinKarpAlgorithm::Test();

	//SuffixArray::Test();

	Trie::Test();
}

