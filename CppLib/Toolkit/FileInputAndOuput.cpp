#include <stdio.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

class FileInputAndOutput
{
	static void ReadAndWriteLineByLineDemo()
	{
		//
		// Example: read intergers in input.txt then output them as float to output.txt
		//

		//
		// 3               <-- N, how many lines that follow
		// 1 2 3           <-- Line 1
		// 4 5 6           <-- Line 2
		// 7 8 9           <-- Line 3
		//

		ifstream inputFile("input.txt");
		string str;

		ofstream outputFile("output.txt");

		getline(inputFile, str);

		int N = 0;
		sscanf(str.c_str(), "%d", &N);

		outputFile << N;
		outputFile << endl;

		for (int i = 0; i < N; i++)
		{
			getline(inputFile, str);

			int d1, d2, d3;
			sscanf(str.c_str(), "%d %d %d", &d1, &d2, &d3);

			outputFile << (float)d1 << " " << (float)d2 << " " << (float)d3;

			if (i < (N - 1))
			{
				outputFile << endl;
			}
		}

		inputFile.close();
		
		outputFile.flush();
		outputFile.close();
	}

	static void ReadAllAndWriteAllDemo()
	{
		//
		// Example: read intergers in input.txt then output them as float to output.txt
		//

		//
		// 3               <-- N, how many lines that follow
		// 1 2 3           <-- Line 1
		// 4 5 6           <-- Line 2
		// 7 8 9           <-- Line 3
		//

		ifstream inputFile("input.txt");
		string str;

		int N = 0;
		sscanf(str.c_str(), "%d", &N);

		vector< vector<int> > intergers(N);

		for (int i = 0; i < N; i++)
		{
			getline(inputFile, str);

			int d1, d2, d3;
			sscanf(str.c_str(), "%d %d %d", &d1, &d2, &d3);

			intergers[i].push_back(d1);
			intergers[i].push_back(d2);
			intergers[i].push_back(d3);
		}

		inputFile.close();

		ofstream outputFile("output.txt");

		for (int i = 0; i < N; i++)
		{
			int size = (int)intergers[i].size();

			for (auto j = 0; j < size; j++)
			{
				outputFile << (float)intergers[i][j] << " ";
				
				if (j < (size - 1))
				{
					outputFile << " ";
				}
			}

			if (i < (N - 1))
			{
				outputFile << endl;
			}
		}

		outputFile.flush();
		outputFile.close();
	}

	static void MapToStandardOutputDemo()
	{
		FILE* pInputFile = fopen("input.txt", "r");

		FILE* pOutputFile = freopen("output.txt", "w", stdout);

		int N = 0;
		fscanf(pInputFile, "%d", &N);  // fscanf is used to read from a file handle

		printf("%d\n", N); // Note that printf, instead of fprintf, is used here since the output file is mapped to stdout 

		for (int i = 0; i < N; i++)
		{
			int d1, d2, d3;
			fscanf(pInputFile, "%d %d %d", &d1, &d2, &d3);

			printf("%f %f %f", (float)d1, (float)d3, (float)d3);

			if (i < (N - 1))
			{
				printf("\n");
			}
		}

		fclose(pInputFile);

		fflush(pOutputFile);
		fclose(pOutputFile);
	}
};


class MyCppLib
{
public:

	//
	// member variables
	//

public:

	MyCppLib()
	{
	}

	void RedirectFileAsStandardInput(const char *pInputFileName)
	{
#ifndef ONLINE_JUDGE

		freopen(pInputFileName, "r", stdin); // redirect input file to stdin

#endif
	}

	void RedirectFileAsStandardOutput(const char *pOutputFileName)
	{
#ifndef ONLINE_JUDGE

		freopen(pOutputFileName, "w", stdout); // redirect output file to stdout

#endif
	}

	void ReadInputByStreamExample(const char *pInputFileName)
	{
		ifstream inputFile(pInputFileName);

		stringstream strStream;
		strStream << inputFile.rdbuf();

		inputFile.close();
	}

	void WriteOutputByStreamExample(const char *pOutputFileName, const char* pText)
	{
		ofstream outputFile(pOutputFileName);

		outputFile << pText;

		outputFile.close();
	}

	void ReadWriteByStreamExample(const char *pInputFileName, const char *pOutputFileName)
	{
		ifstream inputFile(pInputFileName);

		stringstream strStream;
		strStream << inputFile.rdbuf();

		inputFile.close();

		ofstream outputFile(pOutputFileName);

		char buf[1024];

		while (strStream.getline(buf, 1024))
		{
			outputFile << buf;
		}

		outputFile.close();
	}

	void Run()
	{
		printf("MyCppLib is running...\n");
	}
};


