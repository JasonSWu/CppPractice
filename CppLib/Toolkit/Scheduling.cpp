
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <queue>
#include <set>
#include <list>

using namespace std;


#define LARGEST_INT32	(numeric_limits<int>::max())


///////////////////////////////////////////////////////////////////////////////////////////


class FlowJob
{
public:

	//
	// Reference: USACO Training Section 4.2 Job Processing
	//
	static int MinimumProcessingTime(
		int job_num,
		vector<int>& machine_proc_time,		// machine_proc_time[m] means the processing time for m-th machine to handle a job
		vector<int>& finish_time)			// finish_time[j] means the earliest finish time for j-th job
	{
		finish_time.assign(job_num, 0);

		sort(machine_proc_time.begin(), machine_proc_time.end());

		int machine_num = (int)machine_proc_time.size();

		vector<int> busy_until(machine_num, 0);

		for (int j = 0; j < job_num; j++)
		{
			int assigned_machine = -1;
			int min_time = LARGEST_INT32;

			for (int m = 0; m < machine_num; m++)
			{
				if (machine_proc_time[m] + busy_until[m] < min_time)
				{
					min_time = machine_proc_time[m] + busy_until[m];
					assigned_machine = m;
				}
			}

			busy_until[assigned_machine] = min_time;

			finish_time[j] = min_time;
		}

		return finish_time[job_num-1];
	}


	//
	// In light of Johnson's Algorithm, if a job prcoessed by machine group B must be processed by machine group A first,
	// what is the minimum processing time for B?
	//
	static int FlowedMinimumProcessingTime(const vector<int>& finish_time_A, const vector<int>& finish_time_B)
	{
		int job_num = (int)finish_time_A.size();  // assuming A & B have the same same

		int max_time = 0;
		for (int j = 0; j < job_num; j++)
		{
			if (finish_time_A[j] + finish_time_B[job_num - j - 1] > max_time) {
				max_time = finish_time_A[j] + finish_time_B[job_num - j - 1];
			}
		}

		return max_time;
	}


	class Entry
	{
	public:

		int job, ma, mb;

		Entry(int j, int a, int b) {
			job = j; ma = a; mb = b;
		}
	};

	//
	// Johnson's Algorithm, each job must go through machine A then B
	//
	static void JobSequenceWithTwoMachines(
		const vector<int>& machine_A_proc_time,		// machine_A_proc_time[j] means the processing time for machine A to handle j-th job
		const vector<int>& machine_B_proc_time,
		vector<int>& job_sequence,
		int& makespan)
	{
		job_sequence.clear();
		makespan = 0;

		int job_num = (int)machine_A_proc_time.size();
		if (job_num == 0 || job_num != (int)machine_B_proc_time.size()) {
			return;
		}

		list<Entry> jobs_pending;
		for (int j = 0; j < job_num; j++) {
			jobs_pending.push_back( Entry(j, machine_A_proc_time[j], machine_B_proc_time[j]) );
		}

		vector<int> seq_B;

		while (jobs_pending.size() > 0)
		{
			Entry min_A(-1, LARGEST_INT32, -1);
			Entry min_B(-1, -1, LARGEST_INT32);

			list<Entry>::iterator iter_A;
			list<Entry>::iterator iter_B;

			for (list<Entry>::iterator iter = jobs_pending.begin(); iter != jobs_pending.end(); iter++)
			{
				if (min_A.ma > iter->ma) 
				{
					min_A.job = iter->job;
					min_A.ma = iter->ma;
					iter_A = iter;
				}
				if (min_B.mb > iter->mb)
				{
					min_B.job = iter->job;
					min_B.mb = iter->mb;
					iter_B = iter;
				}
			}

			if (min_A.ma < min_B.mb)
			{
				job_sequence.push_back(min_A.job);
				jobs_pending.erase(iter_A);
			}
			else
			{
				seq_B.push_back(min_B.job);
				jobs_pending.erase(iter_B);
			}
		}

		job_sequence.insert(job_sequence.end(), seq_B.rbegin(), seq_B.rend());


		vector<int> finish_A(job_num);
		for (int i = 0; i < (int)job_sequence.size(); i++)
		{
			int j = job_sequence[i];

			finish_A[i] = machine_A_proc_time[j];

			if (i > 0) {
				finish_A[i] += finish_A[i-1];
			}
		}

		vector<int> finish_B(job_num);
		for (int i = 0; i < (int)job_sequence.size(); i++)
		{
			int j = job_sequence[i];

			if (i == 0) {
				finish_B[i] = finish_A[i] + machine_B_proc_time[j];
			}
			else {
				finish_B[i] = max( finish_B[i-1], finish_A[i]) + machine_B_proc_time[j];
			}
		}

		makespan = finish_B[job_num-1];
	}

public:

	static void Test()
	{
		Test_01();

		Test_02();
	}

	static void Test_01()
	{
		int job_num = 5;

		vector<int> machine_proc_time_A = { 1, 1 };
		vector<int> finish_time_A;

		int tA = FlowJob::MinimumProcessingTime(job_num, machine_proc_time_A, finish_time_A);

		bool b = (tA == 3);

		printf("FlowJob::MinimumProcessingTime, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), tA);


		vector<int> machine_proc_time_B = { 3, 1, 4 };
		vector<int> finish_time_B;

		int tB = FlowJob::MinimumProcessingTime(job_num, machine_proc_time_B, finish_time_B);

		b = (tB == 4);

		printf("FlowJob::MinimumProcessingTime, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), tB);



		tB = FlowJob::FlowedMinimumProcessingTime(finish_time_A, finish_time_B);

		b = (tB == 5);

		printf("FlowJob::FlowedMinimumProcessingTime, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), tB);
	}

	static void Test_02()
	{
		vector<int> machine_A_proc_time = { 6, 10, 4, 7, 6, 5 };
		vector<int> machine_B_proc_time = { 4,  8, 9, 2, 3, 6 };

		vector<int> job_sequence;
		int makespan;
		FlowJob::JobSequenceWithTwoMachines(machine_A_proc_time, machine_B_proc_time, job_sequence, makespan);

		vector<int> answer = { 2, 5, 1, 0, 4, 3 };

		bool b = (job_sequence == answer);

		printf("FlowJob::JobSequenceWithTwoMachines, %s \n", (b ? "Correct" : "Wrong"));
	}
};


///////////////////////////////////////////////////////////////////////////////////////////



static void SchedulingTester()
{
	FlowJob::Test();
}


///////////////////////////////////////////////////////////////////////////////////////////
