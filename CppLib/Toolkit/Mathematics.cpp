#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <limits>
#include <bitset>

using namespace std;


///////////////////////////////////////////////////////////////////////////////////////////


//#define INF (1e9)
//#define EPS (1e-9)
//#define PI (acos(-1.0)) 
//
//#define FLOAT_EQUAL(a, b)				(fabs(a-b) <= EPS)
//#define FLOAT_EQUAL_ZERO(a)				(fabs(a) <= EPS)
//#define FLOAT_GREATER_EQUAL_ZERO(a)		(a > -EPS)
//#define FLOAT_LESS_EQUAL_ZERO(a)		(a < EPS)
//
//#define DEG_to_RAD(d)  (d * PI / 180.0)
//#define RAD_to_DEG(r)  (r * 180.0 / PI)

///////////////////////////////////////////////////////////////////////////////////////////


class LargeNumber
{
public:

	int Sign;
	vector<char> Digits;
	map<int, LargeNumber> PowerMap;

public:

	LargeNumber(int n = 0) {
		*this = n;
	}

	LargeNumber(const LargeNumber& other) {
		*this = other;
	}

	LargeNumber(string str) {
		*this = str;
	}

	LargeNumber& operator = (int x)
	{
		Sign = 1;
		Digits.clear();
		PowerMap.clear();

		if (x == 0) 
		{
			Digits.push_back(0);
			return *this;
		}

		Digits.push_back(1);

		return (*this *= x);
	}

	LargeNumber& operator = (const LargeNumber& other)
	{
		Sign = other.Sign;

		if (this != &other) 
		{
			Digits = other.Digits;
			PowerMap = other.PowerMap;
		}

		return *this;
	}

	LargeNumber& operator = (string str)
	{
		Sign = 1;
		Digits.clear();
		PowerMap.clear();

		int len = str.size();
		if (len == 0) {
			return *this;
		}

		const char* pStart = str.c_str();

		if (*pStart == '-') {
			Sign = -1;
			pStart++;
		}

		while (*pStart)
		{
			if (*pStart == '0' || isblank(*pStart)) {
				pStart++;
			}
			else {
				break;
			}
		}

		const char* pEnd = str.c_str() + len - 1;

		for (const char* p = pEnd; p >= pStart; p--) 
		{
			if (*p >= '0' && *p <= '9') {
				Digits.push_back(*p - '0');
			}
		}

		return *this;
	}

	static int CompareAbsolute(const LargeNumber& ln1, const LargeNumber& ln2)
	{
		int s1 = (int)ln1.Digits.size();
		int s2 = (int)ln2.Digits.size();

		if (s1 != s2) {
			return (s1 > s2) ? 1 : -1;
		}

		for (int i = s1 - 1; i >= 0; i--)
		{
			int d1 = ln1.Digits[i];
			int d2 = ln2.Digits[i];

			if (d1 > d2) {
				return 1;
			}
			else if (d1 < d2) {
				return -1;
			}
		}

		return 0;
	}

	static int Compare(const LargeNumber& ln1, const LargeNumber& ln2)
	{
		if (ln1.Sign > 0 && ln2.Sign < 0) {
			return 1;
		}

		if (ln1.Sign < 0 && ln2.Sign > 0) {
			return -1;
		}

		int cmp = CompareAbsolute(ln1, ln2);

		if (cmp > 0) {
			return (ln1.Sign > 0) ? 1 : -1;
		}
		else if (cmp < 0) {
			return (ln1.Sign > 0) ? -1 : 1;
		}

		return 0;
	}

	LargeNumber operator + (const LargeNumber& other) const {
		return LargeNumber(*this) += other;
	}

	LargeNumber& operator += (const LargeNumber& other)	
	{
		if (Sign != other.Sign) 
		{
			LargeNumber tmp(other);
			tmp.Sign = Sign;
			return *this -= tmp;
		}

		int len_self = (int)Digits.size();
		int len_other = (int)other.Digits.size();
		int len_max = max(len_self, len_other);

		int carry = 0;

		for (int i = 0; i < len_max; i++)
		{
			int n1 = (i >= len_self) ? 0 : Digits[i];
			int n2 = (i >= len_other) ? 0 : other.Digits[i];

			carry += n1 + n2;

			if (i >= len_self) {
				Digits.push_back(carry % 10);
				len_self++;
			}
			else {
				Digits[i] = (carry % 10);
			}

			carry /= 10;
		}

		if (carry > 0) {
			Digits.push_back(carry);
		}

		return *this;
	}

	LargeNumber& operator -= (const LargeNumber& _other)
	{
		LargeNumber other(_other);

		if (Sign > 0 && other.Sign < 0) 
		{
			other.Sign = 1;
			return *this += other;
		}

		if (Sign < 0 && other.Sign > 0)
		{
			other.Sign = -1;
			return *this += other;
		}

		int cmp = CompareAbsolute(*this, other);

		if (Sign > 0 && other.Sign > 0)
		{
			if (cmp < 0)
			{
				other = *this;
				*this = _other;
				Sign = -1;
			}
		}
		else if (Sign < 0 && other.Sign < 0)
		{
			if (cmp < 0)
			{
				other = *this;
				*this = _other;
				Sign = 1;
			}
		}

		vector<char>::iterator it1 = Digits.begin();
		vector<char>::const_iterator it2 = other.Digits.begin();

		int dif = 0;

		while (it1 != Digits.end() || it2 != other.Digits.end())
		{
			if (it1 != Digits.end()) {
				dif += *it1;
				++it1;
			}
			if (it2 != other.Digits.end()) {
				dif -= *it2;
				++it2;
			}
			if (dif < 0) {
				*(it1 - 1) = dif + 10;
				dif = -1;
			}
			else {
				*(it1 - 1) = dif % 10;
				dif /= 10;
			}
		}
		
		//if (dif < 0) {
		//	Sign = -1;
		//}

		if (Digits.size() > 1)
		{
			do
			{
				it1 = Digits.end() - 1;
				if (*it1 == 0) {
					Digits.pop_back();
				}
				else {
					break;
				}
			} while (Digits.size() > 1);
		}

		return *this;
	}

	LargeNumber operator - (const LargeNumber& other) const {
		return LargeNumber(*this) -= other;
	}

	//
	// Reference: http://www.geeksforgeeks.org/factorial-large-number/
	//
	LargeNumber& operator *= (int x)
	{
		if (x == 0)
		{
			Sign = 1;
			Digits.clear();
			return *this;
		}

		Sign *= (x < 0) ? -1 : 1;

		if (x < 0) {
			x = -x;
		}

		int carry = 0;  // Initialize carry

		for (int i = 0; i < (int)Digits.size(); i++)
		{
			int prod = Digits[i] * x + carry;
			Digits[i] = prod % 10;  // Store last digit of 'prod'
			carry = prod / 10;		// Put rest in carry
		}

		// Put carry in res and increase result size
		while (carry)
		{
			Digits.push_back(carry % 10);
			carry = carry / 10;
		}

		return *this;
	}

	LargeNumber operator * (int x) const {
		return LargeNumber(*this) *= x;
	}

	LargeNumber operator * (const LargeNumber& other)
	{
		if (other.Digits.size() == 1) {
			return LargeNumber(*this) *= other.Digits[0];
		}

		LargeNumber sum;

		for (int i = 0; i<(int)other.Digits.size(); i++)
		{
			LargeNumber self(*this);

			self *= other.Digits[i];
			for (int k = 0; k < i; k++) {
				self.Digits.insert(self.Digits.begin(), 0);
			}

			sum += self;
		}

		return sum;
	}

	LargeNumber& operator *= (const LargeNumber& other) {
		return *this = (*this * other);
	}

	LargeNumber& operator /= (const LargeNumber& _other)
	{
		int sign = Sign;

		Sign = 1;

		LargeNumber other(_other);
		other.Sign = 1;

		int len = (int)Digits.size();

		LargeNumber tmp(0);

		for (int i = len - 1; i >= 0; i--)
		{
			tmp *= 10;
			tmp += Digits[i];

			Digits[i] = 0;
			
			while ( CompareAbsolute(tmp, other) >= 0 ) 
			{ 
				tmp -= other;
				Digits[i]++; 
			}
		}

		while ((len = (int)Digits.size()) > 0 && Digits[len - 1] == 0) {
			Digits.pop_back();
		}

		Sign = sign;

		if (Digits.size() == 0) {
			Digits.push_back(0);
		}

		return *this;
	}

	LargeNumber operator / (const LargeNumber& other) const {
		return LargeNumber(*this) /= other;
	}

	LargeNumber operator % (const LargeNumber& other) const 
	{
		return *this - (*this / other) * other;
	}

	LargeNumber operator ^(int power)
	{
		if (power == 1) {
			return *this;
		}

		if (PowerMap.count(power)) {
			return PowerMap[power];
		}

		int closestPower = 1;
		while (closestPower < power) {
			closestPower <<= 1;
		}

		closestPower >>= 1;

		if (power == closestPower) {
			PowerMap[power] = (*this ^ (power / 2)) * (*this ^ (power / 2));
		}
		else {
			PowerMap[power] = (*this ^ closestPower) * (*this ^ (power - closestPower));
		}

		return PowerMap[power];
	}

	string ToString() const
	{
		string str;

		if (Sign < 0) {
			str = "-";
		}

		int start = (int)Digits.size() - 1;
		for (int i = start; i >= 0; i--) 
		{
			if (Digits[i] == 0) {
				start--;
			}
			else {
				break;
			}
		}

		for (int i = start; i >= 0; i--) {
			str += Digits[i] + '0';
		}

		if (str.size() == 0 || (str.size() == 1 && str[0] == '-')) {
			str = "0";
		}

		return str;
	}

	static LargeNumber MultiplyRange(int from, int to)
	{
		LargeNumber lnum(from);

		for (int i = from + 1; i <= to; i++) {
			lnum *= i;
		}

		return lnum;
	}

	static LargeNumber Combination(int n, int k) // choose k out of n
	{
		return MultiplyRange(n - k + 1, n) / MultiplyRange(1, k);
	}

public:

	static void Test()
	{
		Test_1();
		Test_2();
		Test_3();
		Test_4();
		Test_5();
		Test_6();
		Test_7();
	}
	
	static void Test_1()
	{
		LargeNumber lnum(123);
		lnum += 12;
		string str = lnum.ToString();
		bool b = (str == "135");

		printf("LargeNumber::+=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


		lnum = -123;
		lnum += -12;
		str = lnum.ToString();
		b = (str == "-135");

		printf("LargeNumber::+=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


		lnum = 123;
		lnum += -12;
		str = lnum.ToString();
		b = (str == "111");

		printf("LargeNumber::+=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


		lnum = -123;
		lnum += 12;
		str = lnum.ToString();
		b = (str == "-111");

		printf("LargeNumber::+=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());
	}

	static void Test_2()
	{
		LargeNumber lnum(123);
		lnum -= 12;
		string str = lnum.ToString();
		bool b = (str == "111");

		printf("LargeNumber::-=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


		lnum = 123;
		lnum -= 124;
		str = lnum.ToString();
		b = (str == "-1");

		printf("LargeNumber::-=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


		lnum = -123;
		lnum -= 124;
		str = lnum.ToString();
		b = (str == "-247");

		printf("LargeNumber::-=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


		lnum = 123;
		lnum -= -124;
		str = lnum.ToString();
		b = (str == "247");

		printf("LargeNumber::-=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());
	}

	static void Test_3()
	{
		LargeNumber lnum(123);
		lnum *= 11;
		string str = lnum.ToString();
		bool b = (str == "1353");

		printf("LargeNumber::*=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


		lnum = 123;
		lnum *= LargeNumber(123);
		str = lnum.ToString();
	    b = (str == "15129");

		printf("LargeNumber::*=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());

		lnum = 100000;
		lnum *= LargeNumber(100000);
		lnum *= LargeNumber(100000);
		lnum *= LargeNumber(100000);
		str = lnum.ToString();
		b = (str == "100000000000000000000");

		printf("LargeNumber::*=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


		lnum = "\r\n     10,000,000,000  \r\n";
		lnum *= LargeNumber("10000000000");
		str = lnum.ToString();
		b = (str == "100000000000000000000");

		printf("LargeNumber::*=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


		lnum = "-00100";
		str = lnum.ToString();
		b = (str == "-100");

		printf("LargeNumber::*=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());
	}

	static void Test_4()
	{
		LargeNumber lnum = LargeNumber::MultiplyRange(1, 100);

		string str = lnum.ToString();

		bool b = (str == "93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000");

		printf("LargeNumber, 100!, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());
	}

	static void Test_5()
	{
		LargeNumber lnum(42);
		lnum /= 2;
		string str = lnum.ToString();
		bool b = (str == "21");

		printf("LargeNumber::/=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


		lnum = -42;
		lnum /= 2;
		str = lnum.ToString();
		b = (str == "-21");

		printf("LargeNumber::/=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


		lnum = 100;
		lnum /= 200;
		str = lnum.ToString();
		b = (str == "0");

		printf("LargeNumber::/=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


		lnum = "100000000000000000000";
		lnum /= LargeNumber("10000000000");
		str = lnum.ToString();
		b = (str == "10000000000");

		printf("LargeNumber::/=, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());
	}

	static void Test_6()
	{
		LargeNumber lnum(30);
		LargeNumber remainder = lnum % 7;
		string str = remainder.ToString();
		bool b = (str == "2");

		printf("LargeNumber::%%, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());

	
		LargeNumber res = LargeNumber::Combination(100, 80);
		str = res.ToString();
		b = (str == "535983370403809682970");

		printf("LargeNumber, Combination(100, 80), %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());
	}

	static void Test_7()
	{
		LargeNumber lnum(2);
		LargeNumber p = lnum ^ 10;
		string str = p.ToString();
		bool b = (str == "1024");

		printf("LargeNumber::^, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


		lnum = 5;
		p = lnum ^ 20;
		str = p.ToString();
		b = (str == "95367431640625");

		printf("LargeNumber::^, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class PrimeNumber
{
public:

	//
	// Sieve of Eratosthenes algorithm, suitable for finding list of primes in [0..upperBound] whith upperBound <= 10^7
	//
	static void SieveOfEratosthenes(long long upperBound, vector<bool>& flags, vector<long long>& primes )
	{
		flags.assign((int)upperBound + 2, true); // set all flags to 1

		flags[0] = flags[1] = 0;

		primes.clear();

		for (long long i = 2; i <= (upperBound+1); i++) 
		{
			if (flags[(int)i])
			{
				// cross out multiples of i starting from i * i
				for (long long j = i * i; j <= (upperBound + 1); j += i)
				{
					flags[(int)j] = 0;
				}
				
				//printf("Sieve Of Eratosthenes: %lld is a prime number\n", i); 

				primes.push_back(i);
			}
		}
	}

	static bool IsPrimeBySieve(long long n, const vector<bool>& flags, const vector<long long>& primes)
	{
		if (n <= flags.size())
		{
			return flags[(int)n];
		}

		for (int i = 0; i < (int)primes.size(); i++)
		{
			if ((n % primes[i]) == 0 )
			{
				return false;
			}
		}

		return true;
	}

	static bool IsPrimeByBruteForce(long long n)
	{
		if (((!(n & 1)) && n != 2) || (n < 2) || (n % 3 == 0 && n != 3))
		{
			return (false);
		}

		//
		// All prime numbers (except 2 and 3) can be expressed in the form 6k+1 or 6k-1, where k is a positive whole number. 
		// This code uses this fact, and tests all numbers in the form of 6k+1 or 6k-1 less than the square root of the number.
		//
		for (long long k = 1; 36 * k*k - 12 * k < n;++k)
		{
			if ((n % (6 * k + 1) == 0) || (n % (6 * k - 1) == 0))
			{
				return false;
			}
		}

		return true;
	}

	static long long Power(long long a, long long n, long long mod)
	{
		long long power = a, result = 1;

		while (n)
		{
			if (n & 1)
			{
				result = (result*power) % mod;
			}
			
			power = (power*power) % mod;
			
			n >>= 1;
		}

		return result;
	}

	static bool WitnessRabinMiller(long long a, long long n)
	{
		long long t, u, i;
		long long prev, curr;

		u = n / 2;
		t = 1;
		while (!(u & 1))
		{
			u /= 2;
			++t;
		}

		prev = Power(a, u, n);

		for (i = 1;i <= t;++i)
		{
			curr = (prev*prev) % n;
			if ((curr == 1) && (prev != 1) && (prev != n - 1))
			{
				return true;
			}
			prev = curr;
		}

		if (curr != 1)
		{
			return true;
		}
		
		return false;
	}

	//
	// if n < 1,373,653, it is enough to test a = 2 and 3.
	// if n < 9,080,191, it is enough to test a = 31 and 73.
	// if n < 4,759,123,141, it is enough to test a = 2, 7, and 61.
	// if n < 2,152,302,898,747, it is enough to test a = 2, 3, 5, 7, and 11.
	// if n < 3,474,749,660,383, it is enough to test a = 2, 3, 5, 7, 11, and 13.
	// if n < 341,550,071,728,321, it is enough to test a = 2, 3, 5, 7, 11, 13, and 17.
	// 
	static bool IsPrimeByRabinMiller(long long number)
	{
		if (((!(number & 1)) && number != 2) || (number < 2) || (number % 3 == 0 && number != 3))
		{
			return false;
		}

		if (number<1373653)
		{
			for (long long k = 1; 36 * k*k - 12 * k < number;++k)
			{
				if ((number % (6 * k + 1) == 0) || (number % (6 * k - 1) == 0))
				{
					return false;
				}
			}

			return true;
		}

		if (number < 9080191)
		{
			if (WitnessRabinMiller(31, number)) return false;
			if (WitnessRabinMiller(73, number)) return false;
			return true;
		}


		if (WitnessRabinMiller(2, number)) return false;
		if (WitnessRabinMiller(7, number)) return false;
		if (WitnessRabinMiller(61, number)) return false;
		
		return true;
	}

	static map<int, int> GetPrimeFactors(long long N)
	{
		vector<bool> flags;
		vector<long long> primes;

		SieveOfEratosthenes(10000000LL, flags, primes);

		map<int, int> dictionary;

		int primeIdx = 0;
		long long PF = primes[primeIdx]; // primes has been populated by sieve

		while (PF * PF <= N)  // stop at sqrt(N); N can get smaller
		{
			while (N % PF == 0) 
			{ 
				N /= PF; 

				if (dictionary.find((int)PF) == dictionary.end())
				{
					dictionary[(int)PF] = 1;
				}
				else
				{
					dictionary[(int)PF] = dictionary[(int)PF] + 1;
				}
			} 
			
			PF = primes[++primeIdx]; // only consider primes!
		}

		if (N != 1)
		{
			dictionary[(int)N] = 1;
		}

		return dictionary; 
	}

public:

	static void TestPrimeNumberSieve()
	{
		vector<bool> flags;
		vector<long long> primes; 

		SieveOfEratosthenes(10000000LL, flags, primes);

		printf("Sieve Of Eratosthenes: 2147483647 is %sa prime number\n", IsPrimeBySieve(2147483647, flags, primes) ? "" : "NOT "); // 10-digits prime
		printf("Sieve Of Eratosthenes: 136117223861 is %sa prime number\n\n", IsPrimeBySieve(136117223861LL, flags, primes) ? "" : "NOT "); // not a prime, 104729*1299709
	}

	static void TestPrimeNumberBruteForce()
	{
		printf("Brute Force: 2147483647 is %sa prime number\n", IsPrimeByBruteForce(2147483647) ? "" : "NOT "); // 10-digits prime
		printf("Brute Force: 136117223861 is %sa prime number\n\n", IsPrimeByBruteForce(136117223861LL) ? "" : "NOT "); // not a prime, 104729*1299709
	}

	static void TestPrimeNumberRabinMiller()
	{
		printf("Rabin-Miller: 2147483647 is %sa prime number\n", IsPrimeByRabinMiller(2147483647) ? "" : "NOT "); // 10-digits prime
		printf("Rabin-Miller: 136117223861 is %sa prime number\n\n", IsPrimeByRabinMiller(136117223861LL) ? "" : "NOT "); // not a prime, 104729*1299709
	}

	static void TestPrimeFactors()
	{
		map<int, int> dictionary = GetPrimeFactors(142391208960LL);

		for (map<int, int>::iterator i = dictionary.begin(); i != dictionary.end(); i++)
		{
			printf("(%d - %d)\n", i->first, i->second);
		}
	}

	static void Test()
	{
		TestPrimeNumberSieve();

		TestPrimeNumberBruteForce();

		TestPrimeNumberRabinMiller();

		TestPrimeFactors();
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class Computing
{
public:

	inline static int min3(int a, int b, int c)
	{
		if (a <= b && a <= c) {
			return a;
		}
		if (b <= a && b <= c) {
			return b;
		}
		return c;
	}

	inline static int max3(int a, int b, int c)
	{
		if (a >= b && a >= c) {
			return a;
		}
		if (b >= a && b >= c) {
			return b;
		}
		return c;
	}

	static unsigned long long ComputeFactorial(unsigned int n)
	{
		if (n == 0)
		{
			return 0;
		}

		long long res = 1;

		for (unsigned int i = 2; i <= n; i++)
		{
			res *= i;
		}

		return res;
	}

	static unsigned long long ComputeCombination(unsigned int n, unsigned int k)
	{
		if (n < k) return 0;
		if (k == 0) return 1;
		if (n == 0) return 0;
		if (k == 1) return n;
		if (n == k) return 1;

		//
		// the formula for "n choose k" is n! / (k! * (n-k)!) 
		//

		unsigned int small = min(k, (n - k));
		unsigned int large = max(k, (n - k));

		pair<unsigned int, unsigned int> x(large+1, n);
		pair<unsigned int, unsigned int> y(1, small);

		long long numerator = 1;
		long long denominator = 1;

		for (unsigned int i = 0; i < small; i++)
		{
			numerator *= x.first + i;
			denominator *= 1 + i;
		}

		return numerator / denominator;
	}

	static void __AggregateMaps(map<unsigned int, unsigned int>& to, const map<unsigned int, unsigned int>& from)
	{
		for (auto it = from.begin(); it != from.end(); ++it)
		{
			auto it2 = to.find(it->first);

			if (it2 != to.end()) {
				it2->second += it->second;
			}
			else {
				to[it->first] = it->second;
			}
		}
	}

	static void __DuelMaps(map<unsigned int, unsigned int>& m1, map<unsigned int, unsigned int>& m2)
	{
		for (auto it1 = m1.begin(); it1 != m1.end(); ++it1)
		{
			auto it2 = m2.find(it1->first);

			if (it2 == m2.end())
			{
				continue;
			}

			int delta = it1->second - it2->second;

			if (delta == 0)
			{
				auto it = it1;
				it1++;
				m1.erase(it);

				m2.erase(it2);
			}
			else if (delta > 0)
			{
				it1->second = delta;
				m2.erase(it2);
			}
			else
			{
				auto it = it1;
				it1++;
				m1.erase(it);

				it2->second = -delta;
			}

			if (it1 == m1.end())
			{
				break;
			}
		}
	}

	static void __GetPrimeFactors(unsigned int N, map<unsigned int, unsigned int>& dictionary)
	{
		static vector<unsigned int> primes =
		{
			2,  3,  5,  7,
			11, 13, 17, 19,
			23, 29,
			31, 37,
			41, 43, 47,
			53, 59,
			61, 67,
			71, 73, 79,
			83, 89,
			97
		};

		dictionary.clear();

		int primeIdx = 0;
		unsigned int PF = primes[primeIdx]; // primes has been populated by sieve

		while (PF * PF <= N)  // stop at sqrt(N); N can get smaller
		{
			while (N % PF == 0)
			{
				N /= PF;

				if (dictionary.find(PF) == dictionary.end())
				{
					dictionary[PF] = 1;
				}
				else
				{
					dictionary[PF] = dictionary[PF] + 1;
				}
			}

			PF = primes[++primeIdx]; // only consider primes!
		}

		if (N != 1)
		{
			dictionary[N] = 1;
		}
	}

	static unsigned long long ComputeLargeCombination(unsigned int n, unsigned int k)
	{
		if (n < k) return 0;
		if (k == 0) return 1;
		if (n == 0) return 0;
		if (k == 1) return n;
		if (n == k) return 1;

		//
		// the formula for "n choose k" is n! / (k! * (n-k)!) 
		//

		unsigned int small = min(k, (n - k));
		unsigned int large = max(k, (n - k));

		pair<unsigned int, unsigned int> x(large + 1, n);
		pair<unsigned int, unsigned int> y(2, small);

		map<unsigned int, unsigned int> factors_x;
		for (auto i = x.first; i <= x.second; i++)
		{
			map<unsigned int, unsigned int> tmp;
			__GetPrimeFactors(i, tmp);

			__AggregateMaps(factors_x, tmp);
		}

		map<unsigned int, unsigned int> factors_y;
		for (auto i = y.first; i <= y.second; i++)
		{
			map<unsigned int, unsigned int> tmp;
			__GetPrimeFactors(i, tmp);

			__AggregateMaps(factors_y, tmp);
		}

		__DuelMaps(factors_x, factors_y);

		unsigned long long numerator = 1;
		for (auto it = factors_x.begin(); it != factors_x.end(); ++it)
		{
			for (unsigned int i = 0; i < it->second; i++) {
				numerator *= it->first;
			}
		}

		unsigned long long denominator = 1;
		for (auto it = factors_y.begin(); it != factors_y.end(); ++it)
		{
			for (unsigned int i = 0; i < it->second; i++) {
				denominator *= it->first;
			}
		}

		return numerator / denominator;
	}

	static int GCD(int a, int b)
	{
		return b == 0 ? a : GCD(b, a % b);
	}

	static int LCM(int a, int b)
	{
		return a * (b / GCD(a, b));
	}

	static int GetReversedInterger(int num)
	{
		int rev = 0;

		do
		{
			int mod = num % 10;
			rev = (rev * 10) + mod;
			num = num / 10;
		} while (num != 0);

		return rev;
	}

public:

	static void TestComputeFactorial()
	{
		printf("\n %d! = %llu \n", 5, ComputeFactorial(5)); // = 120
	}

	static void TestComputeCombination()
	{
		int k = 3;
		int n = 10000;

		printf("\n Choose %d from %d C(%d, %d) is %llu \n", k, n, k, n, ComputeCombination(n, k)); // = 166616670000
	}

	static void TestComputeLargeCombination()
	{
		int k = 26;
		int n = 50;

		printf("\n Choose %d from %d C(%d, %d) is %llu \n", k, n, k, n, ComputeLargeCombination(n, k)); // = 121548660036300
	}

	static void Test()
	{
		TestComputeFactorial();

		TestComputeCombination();

		TestComputeLargeCombination();
	}

};


///////////////////////////////////////////////////////////////////////////////////////////


class RomanNumber
{
public:

	static string ConvertToRoman(int v)
	{
		const pair<int, string> rule[] = {
			{ 1000, "M" },{ 900, "CM" },
			{ 500, "D" },{ 400, "CD" },
			{ 100, "C" },{ 90, "XC" },
			{ 50, "L" },{ 40, "XL" },
			{ 10, "X" },{ 9, "IX" },
			{ 5, "V" },{ 4, "IV" },
			{ 1, "I" },
			{ 0, "" } 
		};

		string str;
		for (const pair<int, string>* p = rule; p->first > 0; ++p) {
			while (v >= p->first) {
				str += p->second;
				v -= p->first;
			}
		}

		return str;
	}

public:

	static void Test()
	{
		printf("3549 => %s \n", ConvertToRoman(3549).c_str());	// MMMDXLIX
	}
};



///////////////////////////////////////////////////////////////////////////////////////////


class FibonacciNumber
{
public:

	static long long CalculateFibonacciNumber(int n)
	{
		if (n <= 1) {
			return n;
		}

		long long a = 0, b = 1;
		long long sum = a + b;

		for(int i=2; i<n;i++)
		{
			a = b;
			b = sum;
			sum = a + b;
		}

		return sum;
	}

public:

	static void Test()
	{
		printf("Fibonacci(%d) = %lld \n", 20, CalculateFibonacciNumber(20));	
		printf("Fibonacci(%d) = %lld \n", 50, CalculateFibonacciNumber(50));	
		printf("Fibonacci(%d) = %lld \n", 80, CalculateFibonacciNumber(80));	
	}
};



///////////////////////////////////////////////////////////////////////////////////////////


class CatalanNumber
{
public:

	static long long CalculateCatalanNumber(int n)
	{
		if (n <= 0) {
			return 1;
		}

		long long sum = 0;

		for (int i = 0; i<n/2; i++) {
			sum += CalculateCatalanNumber(i) * CalculateCatalanNumber(n - i - 1);
		}

		sum *= 2;

		if (n % 2) {
			long long c = CalculateCatalanNumber(n / 2);
			sum += c * c;
		}

		return sum;
	}

public:

	static void Test()
	{
		printf("Catalan(%d) = %lld \n", 10, CalculateCatalanNumber(10));
		printf("Catalan(%d) = %lld \n", 20, CalculateCatalanNumber(20));
	}
};



///////////////////////////////////////////////////////////////////////////////////////////


long long FastPower(int base, int exp, int Mod = 1)
{
	if (base == 0) {
		return 0;
	}

	if (exp == 0) {
		return 1;
	}

	base %= Mod;
	long long half = FastPower(base, exp / 2, Mod);

	long long square = (half * half) % Mod;

	if ((exp % 2) == 0) {
		return square;
	}

	if (exp > 0) {
		return (base * square) % Mod;
	}

	return square / base;
}


///////////////////////////////////////////////////////////////////////////////////////////


static void MathematicsTester()
{
	LargeNumber::Test();

	//PrimeNumber::Test();

	//Computing::Test();

	//RomanNumber::Test();

	//FibonacciNumber::Test();

	//CatalanNumber::Test();
}


///////////////////////////////////////////////////////////////////////////////////////////

