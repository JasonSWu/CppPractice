#include <stdio.h>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <set>
#include <stack>
#include <queue>


using namespace std;



class Digraph
{
public:

	int V, E;
	vector< vector<int> > AdjList;
	vector<int> Indegree;
	vector<bool> Exist;

	Digraph()
	{
		V = E = 0;
	}

	void CreateAdjList(int v)
	{
		E = 0;
		V = v;
		AdjList.assign(V, vector<int>());
		Indegree.assign(V, 0);
		Exist.assign(V, false);
	}

	void AddOneWayEdge(int from, int to)
	{
		E++;
		AdjList[from].push_back(to);
		Indegree[to]++;
		Exist[from] = true;
		Exist[to] = true;
	}

	void AddTwoWayEdge(int from, int to)
	{
		E += 2;
		AdjList[from].push_back(to);
		Indegree[to]++;
		AdjList[to].push_back(from);
		Indegree[from]++;
		Exist[from] = true;
		Exist[to] = true;
	}

	void CreateReverseFrom(const Digraph& other)
	{
		CreateAdjList(other.V);

		for (int from = 0; from < (int)other.AdjList.size(); from++)
		{
			const vector<int>& edges = other.AdjList[from];

			for (int i = 0; i < (int)edges.size(); i++)
			{
				AddOneWayEdge(edges[i], from);
			}
		}
	}
};


class DigraphTopologyHandler
{
public:

	bool IsCyclic;
	vector<int> Visited;

	DigraphTopologyHandler()
	{
	}

	void Prepare(Digraph& graph)
	{
		IsCyclic = false;
		Visited.assign(graph.V, false);
	}

	bool CheckCycleUtil(Digraph& graph, int startVertex, vector<bool>& flags)
	{
		if (!Visited[startVertex])
		{
			Visited[startVertex] = true;
			flags[startVertex] = true;

			for (int j = 0; j < (int)graph.AdjList[startVertex].size(); j++)
			{
				int nextVertex = graph.AdjList[startVertex][j];

				if (!Visited[nextVertex] && CheckCycleUtil(graph, nextVertex, flags)) {
					return true;
				}
				else if (flags[nextVertex]) {
					return true;
				}
			}
		}

		flags[startVertex] = false;

		return false;
	}

	bool CheckCycle(Digraph& graph)
	{
		Prepare(graph);

		vector<bool> flags(graph.V, false);

		for (int i = 0; i < graph.V; i++)
		{
			if (CheckCycleUtil(graph, i, flags))
			{
				IsCyclic = true;
				return true;
			}
		}

		return false;
	}

	void DfsRecursive(Digraph& graph, int startVertex, vector<int>* pPreOrder = 0, vector<int>* pPostOrder = 0, stack<int>* pReversePostOrder = 0)
	{
		if (Visited[startVertex]) {
			return;
		}

		if (pPreOrder != 0) {
			pPreOrder->push_back(startVertex);
		}

		Visited[startVertex] = true;

		for (int j = 0; j < (int)graph.AdjList[startVertex].size(); j++)
		{
			int nextVertex = graph.AdjList[startVertex][j];

			if (!Visited[nextVertex]) {
				DfsRecursive(graph, nextVertex, pPreOrder, pPostOrder, pReversePostOrder);
			}
		}

		if (pPostOrder != 0) {
			pPostOrder->push_back(startVertex);
		}

		if (pReversePostOrder != 0) {
			pReversePostOrder->push(startVertex);
		}
	}

	int LeastCommonParent_DfsRecursive(Digraph& graph, int startVertex, int v1, int v2)
	{
		if (startVertex == v1 || startVertex == v2) {
			return startVertex;
		}

		if (Visited[startVertex]) {
			return -1;
		}

		Visited[startVertex] = true;

		int x1 = -1, x2 = -1;

		for (int j = 0; j < (int)graph.AdjList[startVertex].size(); j++)
		{
			int nextVertex = graph.AdjList[startVertex][j];

			int x = LeastCommonParent_DfsRecursive(graph, nextVertex, v1, v2);

			if (x >= 0) {
				if (x1 >= 0) x2 = x; else x1 = x;
			}

			if (x1 >= 0 && x2 >= 0) {
				break;
			}
		}

		if (x1 >= 0 && x2 >= 0) {
			return startVertex;
		}

		if (x1 >= 0 && x2 < 0) {
			return x1;
		}

		if (x1 < 0 && x2 >= 0) {
			return x2;
		}

		return -1;
	}

	void TopologicalSort_DfsRecursive(Digraph& graph, vector<int>& topology)
	{
		topology.clear();

		Prepare(graph);

		stack<int> reversePostOrder;

		for (int i = 0; i < graph.V; i++)
		{
			if (!Visited[i])
			{
				vector<int> topo;
				DfsRecursive(graph, i, 0, 0, &reversePostOrder);
			}
		}

		while (!reversePostOrder.empty())
		{
			topology.push_back(reversePostOrder.top());
			reversePostOrder.pop();
		}
	}

	void DfsNonRecursive(Digraph& graph, int startVertex, vector<int>* pTopology = 0)
	{
		if (pTopology != 0) {
			pTopology->clear();
		}

		if (Visited[startVertex]) {
			return;
		}

		stack<int> s;

		s.push(startVertex);
		Visited[startVertex] = true;

		while (!s.empty())
		{
			int currentVertex = s.top();
			s.pop();

			if (pTopology != 0) {
				pTopology->push_back(currentVertex);
			}

			vector<int>& adjs = graph.AdjList[currentVertex];
			int size = (int)adjs.size();

			for (int j = 0; j < size; j++)
			{
				int nextVertex = adjs[j];

				if (!Visited[nextVertex])
				{
					s.push(nextVertex);
					Visited[nextVertex] = true;
				}
			}
		}
	}

	void TopologicalSort_DfsNonRecursive(Digraph& graph, vector<int>& topology)
	{
		topology.clear();

		Prepare(graph);

		for (int i = 0; i < graph.V; i++)
		{
			if (!Visited[i])
			{
				vector<int> topo;
				DfsNonRecursive(graph, i, &topo);

				topology.insert(topology.begin(), topo.begin(), topo.end());
			}
		}
	}

	void LexigraphicalSort_PriorityQueue(Digraph& graph, vector<int>& topology)
	{
		topology.assign(graph.V, -1);

		Prepare(graph);

		vector<int> indegree = graph.Indegree;

		priority_queue<int> q;

		// initialize the queue with zero-indegree nodes (nodes not being dependent on)
		for (int i = 0; i < graph.V; i++)
		{
			if (indegree[i] == 0)
			{
				// priority_queue returns max but we need the min
				q.push(-i);
			}
		}

		for (int i = 0; i < graph.V; i++)
		{
			if (q.empty())
			{
				// nothing in queue, topological sort is impossible.

				IsCyclic = true;
				return;
			}

			int v = -q.top();
			q.pop();

			topology[i] = v;

			for (int j = 0; j < (int)graph.AdjList[v].size(); j++)
			{
				int next = graph.AdjList[v][j];

				indegree[next]--;

				if (indegree[next] == 0) {
					q.push(-next);
				}
			}
		}
	}

	bool CheckCycle_LexigraphicalSort_PriorityQueue(Digraph& graph)
	{
		vector<int> order;
		LexigraphicalSort_PriorityQueue(graph, order);

		return IsCyclic;
	}

	static void Test1_CheckCycle_1()
	{
		//
		// https://www.geeksforgeeks.org/detect-cycle-undirected-graph/
		//

		Digraph graph;

		graph.CreateAdjList(5);

		graph.AddOneWayEdge(1, 0);
		graph.AddOneWayEdge(0, 2);
		graph.AddOneWayEdge(2, 0);
		graph.AddOneWayEdge(0, 3);
		graph.AddOneWayEdge(3, 4);

		DigraphTopologyHandler handler;

		bool bHasCycle1 = handler.CheckCycle(graph); // has cycle
		bool bHasCycle2 = handler.CheckCycle_LexigraphicalSort_PriorityQueue(graph); // has cycle
	}

	static void Test1_CheckCycle_2()
	{
		//
		// USACO 2018-03 Gold MilkingOrder Data 1
		//

		Digraph graph;

		graph.CreateAdjList(4);

		graph.AddOneWayEdge(0, 1);
		graph.AddOneWayEdge(1, 2);
		graph.AddOneWayEdge(2, 3);
		graph.AddOneWayEdge(3, 1);

		graph.AddOneWayEdge(2, 3);
		graph.AddOneWayEdge(3, 0);

		DigraphTopologyHandler handler;

		bool bHasCycle1 = handler.CheckCycle(graph); // has cycle
		bool bHasCycle2 = handler.CheckCycle_LexigraphicalSort_PriorityQueue(graph); // has cycle
	}

	static void Test1_CheckCycle_3()
	{
		Digraph graph;

		graph.CreateAdjList(10);

		graph.AddOneWayEdge(5, 9);
		graph.AddOneWayEdge(9, 8);
		graph.AddOneWayEdge(9, 2);
		graph.AddOneWayEdge(8, 1);
		graph.AddOneWayEdge(1, 2);

		DigraphTopologyHandler handler;

		bool bHasCycle1 = handler.CheckCycle(graph); // doesn't have cycle
		bool bHasCycle2 = handler.CheckCycle_LexigraphicalSort_PriorityQueue(graph); // doesn't have cycle


		graph.CreateAdjList(10);

		graph.AddOneWayEdge(5, 9);
		graph.AddOneWayEdge(9, 8);
		graph.AddOneWayEdge(8, 1);
		graph.AddOneWayEdge(1, 2);
		graph.AddOneWayEdge(2, 9);

		bHasCycle1 = handler.CheckCycle(graph); // has cycle
		bHasCycle2 = handler.CheckCycle_LexigraphicalSort_PriorityQueue(graph); // has cycle
	}

	static void Test_TopologicalSort_1()
	{
		//
		// https://www.geeksforgeeks.org/topological-sorting/
		//

		Digraph graph;

		graph.CreateAdjList(6);

		graph.AddOneWayEdge(5, 2);
		graph.AddOneWayEdge(5, 0);
		graph.AddOneWayEdge(4, 0);
		graph.AddOneWayEdge(4, 1);
		graph.AddOneWayEdge(2, 3);
		graph.AddOneWayEdge(3, 1);

		vector<int> topology;
		DigraphTopologyHandler handler;

		handler.TopologicalSort_DfsNonRecursive(graph, topology); // topology is 5 4 2 3 1 0

		handler.TopologicalSort_DfsRecursive(graph, topology); // topology is 5 4 2 3 1 0
	}

	static void Test_LexigraphicalSort_1()
	{
		//
		// https://www.geeksforgeeks.org/topological-sorting/
		//

		Digraph graph;

		graph.CreateAdjList(6);

		graph.AddOneWayEdge(5, 2);
		graph.AddOneWayEdge(5, 0);
		graph.AddOneWayEdge(4, 0);
		graph.AddOneWayEdge(4, 1);
		graph.AddOneWayEdge(2, 3);
		graph.AddOneWayEdge(3, 1);

		vector<int> order;
		DigraphTopologyHandler handler;

		handler.LexigraphicalSort_PriorityQueue(graph, order); // order is 4 5 0 2 3 1, no cycle
	}

	static void Test_LexigraphicalSort_2()
	{
		//
		// USACO 2018-03 Gold MilkingOrder Data 1
		//

		Digraph graph;

		graph.CreateAdjList(4);

		graph.AddOneWayEdge(0, 1);
		graph.AddOneWayEdge(1, 2);
		graph.AddOneWayEdge(3, 1);


		vector<int> order;
		DigraphTopologyHandler handler;

		handler.LexigraphicalSort_PriorityQueue(graph, order); // order is 0 3 1 2, doesn't have cycle

		graph.AddOneWayEdge(2, 3);
		graph.AddOneWayEdge(3, 0);

		handler.LexigraphicalSort_PriorityQueue(graph, order); // has cycle
	}

	static void Test_LexigraphicalSort_3()
	{
		//
		// USACO 2018-03 Gold MilkingOrder Data 3
		//
		//  10 4
		// 	3 6 10 3
		// 	4 10 9 2 3
		// 	3 9 7 8
		// 	4 10 7 9 5

		Digraph graph;

		graph.CreateAdjList(10);

		graph.AddOneWayEdge(5, 9);
		graph.AddOneWayEdge(9, 2);

		graph.AddOneWayEdge(9, 8);
		graph.AddOneWayEdge(8, 1);
		graph.AddOneWayEdge(1, 2);

		graph.AddOneWayEdge(8, 6);
		graph.AddOneWayEdge(6, 7);

		vector<int> order;
		DigraphTopologyHandler handler;

		handler.LexigraphicalSort_PriorityQueue(graph, order); // order 0 3 4 5 9 8 1 2 6 7, doesn't have cycle

		graph.AddOneWayEdge(9, 6);
		graph.AddOneWayEdge(6, 8);
		graph.AddOneWayEdge(8, 4);

		handler.LexigraphicalSort_PriorityQueue(graph, order); // has cycle
	}

	static void Test_LeastCommonParent_DfsRecursive_01()
	{
		Digraph graph;

		graph.CreateAdjList(8);

		graph.AddTwoWayEdge(1, 2);
		graph.AddTwoWayEdge(1, 3);
		graph.AddTwoWayEdge(3, 4);
		graph.AddTwoWayEdge(3, 5);
		graph.AddTwoWayEdge(4, 6);
		graph.AddTwoWayEdge(5, 7);

		DigraphTopologyHandler handler;
		handler.Prepare(graph);

		int p = handler.LeastCommonParent_DfsRecursive(graph, 1, 4, 7); // = 3

	}

	static void Test()
	{
		Test1_CheckCycle_1();
		Test1_CheckCycle_2();
		Test1_CheckCycle_3();
		Test_TopologicalSort_1();
		Test_LexigraphicalSort_1();
		Test_LexigraphicalSort_2();
		Test_LexigraphicalSort_3();
		Test_LeastCommonParent_DfsRecursive_01();
	}
};



///////////////////////////////////////////////////////////////////////////////////////////


static void DigraphTester()
{
	DigraphTopologyHandler::Test();
}


///////////////////////////////////////////////////////////////////////////////////////////


