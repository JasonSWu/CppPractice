#include <algorithm>
#include <cstdio>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <list>
#include <stack>
#include <limits>
#include <cstring>

using namespace std;



///////////////////////////////////////////////////////////////////////////////////////////


#define EPS (1e-9)

#define SMALLEST_DOUBLE				(numeric_limits<double>::lowest())
#define LARGEST_DOUBLE				(numeric_limits<double>::max())

#define FLOAT_EQUAL(a, b)			(fabs(a-b)<=EPS)
#define FLOAT_EQUAL_ZERO(a)			(fabs(a) <= EPS)

#ifndef NULL
	#define NULL 0
#endif


///////////////////////////////////////////////////////////////////////////////////////////


struct FlowEdge
{
	int idx;  // for simplicity, index of this edge in some vector, etc.
	int v_from;
	int v_to;
	double capacity;
	double flow;

	FlowEdge(int from = -1, int to = -1, double c = -1, int index = -1) {
		v_from = from;
		v_to = to;
		capacity = c;
		flow = 0;
		idx = index;
	}

	bool operator < (const FlowEdge& other) const {
		return capacity < other.capacity;
	}

	int OtherVertex(int v)
	{
		if (v == v_from) {
			return v_to;
		}
		else if (v == v_to) {
			return v_from;
		}

		return -1;  // should not happen!
	}

	double ResidualCapacityToVertex(int v) const
	{
		if (v_from == v) {		// forward edge
			return flow;
		}
		else if (v_to == v) {	// backward edge
			// for super source/sink, keep the capacity as infinite
			return FLOAT_EQUAL(LARGEST_DOUBLE, capacity) ? LARGEST_DOUBLE : (capacity - flow);
		}

		return -1;  // should not happen!
	}

	void AddResidualFlowToVertex(int v, double delta)
	{
		if (v_from == v) {		// forward edge	
			flow -= delta;
		}
		else if (v_to == v) {   // backward edge
			flow += delta;;
		}
	}

};



///////////////////////////////////////////////////////////////////////////////////////////


class FlowNetwork
{
public:

	int V, E;

	vector< FlowEdge > vecEdges;
	FlowEdge* Edges;
	vector< set< int > > vecAdjEdgeSets;   // vecAdjEdgeSets[v] stores forward and backward edges incident to v

public:

	FlowNetwork()
	{
		Clear();
	}

	FlowNetwork(const FlowNetwork& other)
	{
		*this = other;
	}

	FlowNetwork& operator=(const FlowNetwork& other)
	{
		if (&other != this)
		{
			V = other.V;
			E = other.E;

			vecEdges = other.vecEdges;
			vecAdjEdgeSets = other.vecAdjEdgeSets;
		}

		return *this;
	}

	void Clear()
	{
		V = 0;
		E = 0;

		vecEdges.clear();
		vecAdjEdgeSets.clear();
	}

	void AddEdge(int from, int to, double capacity) 
	{
		vecEdges.push_back(FlowEdge(from, to, capacity, (int)vecEdges.size()));

		E = (int)vecEdges.size();

		V = max(V, max(from, to));
	}

	void AddSuperEdge(int from, int to)
	{
		vecEdges.push_back(FlowEdge(from, to, LARGEST_DOUBLE, (int)vecEdges.size()));

		E = (int)vecEdges.size();

		V = max(V, max(from, to));
	}

	void BuildAdjList()
	{
		Edges = &vecEdges.front();

		vecAdjEdgeSets.assign(V, set<int>());

		int size = (int)vecEdges.size();
		for (int e = 0; e < size; e++) {
			FlowEdge& edge = Edges[e];
			vecAdjEdgeSets[edge.v_from].insert(e);
			vecAdjEdgeSets[edge.v_to].insert(e);
		}
	}

	void ResetFlow()
	{
		Edges = &vecEdges.front();

		int size = (int)vecEdges.size();
		for (int i = 0; i < size; i++) {
			Edges[i].flow = 0;
		}
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class FlowHandler
{
public:

	vector<int> vecPath;	// vecPath[v] stores edge index
	int* Path;
	vector<unsigned char> vecVisited;
	unsigned char* Visited;

	vector< vector<int> > vecConnectedEdgeIndices;

	vector<int> MinCutVertices;
	vector<FlowEdge> MinCutEdges;

	double MaxFlow;
	int S, T;

public:

	FlowHandler()
	{
		MaxFlow = 0;
		S = T = -1;
	}

	//
	// Worst Case: O(E*E*V)
	//
	bool FindShortestPathByBfs(FlowNetwork& network, int from, int to)
	{
		if (vecPath.size() != network.V) {
			vecPath.assign(network.V, 0);
		}

		if (vecVisited.size() != network.V) {
			vecVisited.assign(network.V, false);
		}

		Path = &vecPath.front();
		memset(Path, 0, sizeof(Path[0])* vecPath.size());

		Visited = &vecVisited.front();
		memset(Visited, 0, sizeof(Visited[0])* vecVisited.size());

		if (vecConnectedEdgeIndices.size() != network.V)
		{
			vecConnectedEdgeIndices.assign(network.V, vector<int>());

			set<int>::iterator iter;

			for (int v = 0; v < network.V; v++)
			{
				vector<int>& edgeIndices = vecConnectedEdgeIndices[v];

				set<int>& adjs = network.vecAdjEdgeSets[v];

				for (iter = adjs.begin(); iter != adjs.end(); iter++) {
					edgeIndices.push_back(*iter);
				}
			}
		}

		network.Edges = &network.vecEdges.front();

		queue<int> q;

		Visited[from] = true;
		q.push(from);

		while (!q.empty())
		{
			int v = q.front();
			q.pop();

			vector<int>& vecEdgeIndices = vecConnectedEdgeIndices[v];
			int* edgeIndices = &vecEdgeIndices.front();

			int size = (int)vecEdgeIndices.size();
			for (int e = 0; e < size; e++)
			{
				int edge_idx = edgeIndices[e];

				FlowEdge& edge = network.Edges[edge_idx];

				int w = (edge.v_from == v) ? edge.v_to : edge.v_from;  // performance tweak, same as edge.OtherVertex(v);

				if (!Visited[w])
				{
					double residual = (edge.v_from == w) ? edge.flow : edge.capacity - edge.flow;  // performance tweak, same as edge.ResidualCapacityToVertex(w) > 0

					if (residual > 0)
					{
						Path[w] = edge_idx;

						Visited[w] = true;
						q.push(w);
					}
				}
			}
		}

		return Visited[to];
	}

	//
	// Worst Case: O(E*E*LogC), C is the max capacity
	//
	bool FindFattestPathByDijkstra(FlowNetwork& network, int from, int to)
	{
		if (vecPath.size() != network.V) {
			vecPath.assign(network.V, 0);
		}

		if (vecVisited.size() != network.V) {
			vecVisited.assign(network.V, false);
		}

		Path = &vecPath.front();
		memset(Path, 0, sizeof(Path[0])* vecPath.size());

		Visited = &vecVisited.front();
		memset(Visited, 0, sizeof(Visited[0])* vecVisited.size());

		if (vecConnectedEdgeIndices.size() != network.V)
		{
			vecConnectedEdgeIndices.assign(network.V, vector<int>());

			set<int>::iterator iter;

			for (int v = 0; v < network.V; v++)
			{
				vector<int>& edgeIndices = vecConnectedEdgeIndices[v];

				set<int>& adjs = network.vecAdjEdgeSets[v];

				for (iter = adjs.begin(); iter != adjs.end(); iter++) {
					edgeIndices.push_back(*iter);
				}
			}
		}

		network.Edges = &network.vecEdges.front();

		vector<FlowEdge> traversal(network.V, FlowEdge(-1, -1, 0));
		traversal[from].capacity = LARGEST_DOUBLE;

		struct EdgeCompare {
			bool operator()(const FlowEdge& l, const FlowEdge& r) {
				return l.capacity < r.capacity;
			}
		};

		priority_queue<FlowEdge, vector<FlowEdge>, EdgeCompare> q;

		q.push(FlowEdge(from, from, 0));
		Visited[from] = true;

		while (!q.empty())
		{
			FlowEdge current_edge = q.top();
			q.pop();

			set<int>& adjs = network.vecAdjEdgeSets[current_edge.v_to];

			vector<int>& vecEdgeIndices = vecConnectedEdgeIndices[current_edge.v_to];
			int* edgeIndices = &vecEdgeIndices.front();

			int size = (int)vecEdgeIndices.size();
			for (int e = 0; e < size; e++)
			{
				int edge_idx = edgeIndices[e];

				FlowEdge& edge = network.Edges[edge_idx];

				int w = edge.OtherVertex(current_edge.v_to);

				double residual = edge.ResidualCapacityToVertex(w);

				if (residual > 0 && traversal[w].capacity < min(traversal[current_edge.v_to].capacity, residual))
				{
					traversal[w].capacity = min(traversal[current_edge.v_to].capacity, residual);

					Path[w] = edge_idx;

					Visited[w] = true;

					q.push(FlowEdge(current_edge.v_to, w, residual));
				}
			}
		}

		return Visited[to];
	}

	void FordFulkersonAlgorithm(FlowNetwork& network, int s, int t)
	{
		MaxFlow = 0;
		MinCutVertices.clear();
		MinCutEdges.clear();

		S = s;
		T = t;

		while (FindShortestPathByBfs(network, s, t))
		{
			double bottleNeck = LARGEST_DOUBLE;
			for (int v = t; v != s; v = network.vecEdges[Path[v]].OtherVertex(v))	{
				bottleNeck = min(bottleNeck, network.vecEdges[Path[v]].ResidualCapacityToVertex(v));
			}

			for (int v = t; v != s; v = network.vecEdges[Path[v]].OtherVertex(v)) {
				network.vecEdges[Path[v]].AddResidualFlowToVertex(v, bottleNeck);
			}

			MaxFlow += bottleNeck;
		}

		ResolveMinCutVertices(network);

		ResolveMinCutEdges(network);
	}

	void ResolveMinCutVertices(FlowNetwork& network)
	{
		MinCutVertices.clear();

		int size = (int)network.vecAdjEdgeSets.size();
		for (int v = 0; v < size; v++) {
			if (IsMinCutVertex(v)) {
				MinCutVertices.push_back(v);
			}
		}
	}

	void ResolveMinCutEdges(FlowNetwork& network)
	{
		MinCutEdges.clear();

		int size = (int)network.vecAdjEdgeSets.size();
		for (int from = 0; from < size; from++)
		{
			set<int>& adjs = network.vecAdjEdgeSets[from];

			for (auto iter = adjs.begin(); iter != adjs.end(); iter++)
			{
				FlowEdge& edge = network.vecEdges[*iter];

				if (edge.v_from != from) {
					continue;
				}

				int to = edge.v_to;

				if (IsMinCutVertex(from) && !IsMinCutVertex(to) && edge.capacity > 0) {
					MinCutEdges.push_back(edge);
				}
			}
		}
	}

	bool IsMinCutVertex(int v) const {
		return Visited[v];
	}

	double MinCutCapacity() const {
		return MaxFlow;
	}

public:

	static void Test()
	{
		Test_1();
		Test_2();
		Test_3();
	}

	//
	// Reference: Addison Wesley - Algorithms, 4th Edition
	//
	static void Test_1()
	{
		//	6			<- V
		//	8			<- E
		//	0 1 2.0
		//	0 2 3.0
		//	1 3 3.0
		//	1 4 1.0
		//	2 3 1.0
		//	2 4 1.0
		//	3 5 2.0
		//	4 5 3.0

		FlowNetwork network;

		network.V = 6;
		network.E = 8;

		network.AddEdge(0, 1, 2.0);
		network.AddEdge(0, 2, 3.0);
		network.AddEdge(1, 3, 3.0);
		network.AddEdge(1, 4, 1.0);
		network.AddEdge(2, 3, 1.0);
		network.AddEdge(2, 4, 1.0);
		network.AddEdge(3, 5, 2.0);
		network.AddEdge(4, 5, 3.0);

		network.BuildAdjList();

		FlowHandler handler;

		handler.FordFulkersonAlgorithm(network, 0, 5);

		bool b = FLOAT_EQUAL(handler.MaxFlow, 4.0);

		printf("FlowHandler, %s, the answer is %lf \n", (b == true) ? "Correct" : "Wrong", handler.MaxFlow);
	}

	//
	// Reference: http://www.geeksforgeeks.org/ford-fulkerson-algorithm-for-maximum-flow-problem/
	//
	static void Test_2()
	{
		//	6			<- V
		//	10			<- E
		//	0 1 16.0
		//	0 2 13.0
		//	1 2 10.0
		//	1 3 12.0
		//	2 1 4.0
		//	2 4 14.0
		//	3 2 9.0
		//	3 5 20.0
		//	4 3 7.0
		//	4 5 4.0

		FlowNetwork network;

		network.V = 6;
		network.E = 10;

		network.AddEdge(0, 1, 16.0);
		network.AddEdge(0, 2, 13.0);
		network.AddEdge(1, 2, 10.0);
		network.AddEdge(1, 3, 12.0);
		network.AddEdge(2, 1, 4.0);
		network.AddEdge(2, 4, 14.0);
		network.AddEdge(3, 2, 9.0);
		network.AddEdge(3, 5, 20.0);
		network.AddEdge(4, 3, 7.0);
		network.AddEdge(4, 5, 4.0);

		network.BuildAdjList();

		FlowHandler handler;

		handler.FordFulkersonAlgorithm(network, 0, 5);

		bool b = FLOAT_EQUAL(handler.MaxFlow, 23.0);

		printf("FlowHandler, %s, the answer is %lf \n", (b == true) ? "Correct" : "Wrong", handler.MaxFlow);
	}

	static void Test_3()
	{
		//	8			<- V
		//	10			<- E
		//	0 1 2.0
		//	0 2 3.0
		//	1 3 3.0
		//	1 4 1.0
		//	2 3 1.0
		//	2 4 1.0
		//	3 5 2.0
		//	4 5 3.0

		//  6 0 INFINITE
		//  5 7 INFINITE

		FlowNetwork network;

		network.V = 8;
		network.E = 10;

		network.AddEdge(0, 1, 2.0);
		network.AddEdge(0, 2, 3.0);
		network.AddEdge(1, 3, 3.0);
		network.AddEdge(1, 4, 1.0);
		network.AddEdge(2, 3, 1.0);
		network.AddEdge(2, 4, 1.0);
		network.AddEdge(3, 5, 2.0);
		network.AddEdge(4, 5, 3.0);

		network.AddSuperEdge(6, 0);
		network.AddSuperEdge(5, 7);

		network.BuildAdjList();

		FlowHandler handler;

		handler.FordFulkersonAlgorithm(network, 6, 7);

		bool b = FLOAT_EQUAL(handler.MaxFlow, 4.0);

		printf("FlowHandler, %s, the answer is %lf \n", (b == true) ? "Correct" : "Wrong", handler.MaxFlow);
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


static void FlowNetworkTester()
{
	FlowHandler::Test();
}


///////////////////////////////////////////////////////////////////////////////////////////
