#include <iostream>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <list>
#include <stack>
#include <limits>

using namespace std;



///////////////////////////////////////////////////////////////////////////////////////////


#define EPS (1e-9)

#define SMALLEST_DOUBLE				(numeric_limits<double>::lowest())
#define LARGEST_DOUBLE				(numeric_limits<double>::max())

#define FLOAT_EQUAL(a, b)  (fabs(a-b)<=EPS)

#ifndef NULL
    #define NULL 0
#endif


///////////////////////////////////////////////////////////////////////////////////////////


class GraphUnionFind
{
public:

    vector<int> Parents;
    vector<int> ElementCounts;

    int TotalSetCount;

public:

    GraphUnionFind(int setCount = 0)
    {
        Initialize(setCount);
    }

    void Clear()
    {
        this->Parents.clear();
        this->ElementCounts.clear();

        this->TotalSetCount = 0;
    }

    void Initialize(int setCount)
    {
        Clear();

        this->TotalSetCount = setCount;

        if (this->TotalSetCount > 0)
        {
            this->ElementCounts.assign(this->TotalSetCount, 1);

            this->Parents.assign(this->TotalSetCount, 0);

            for (int i = 0; i < this->TotalSetCount; i++)
            {
                this->Parents[i] = i;
            }
        }
    }

    int FindSet(int e) const
    {
        while (e != this->Parents[e])
        {
            e = this->Parents[e];
        }

        return e;
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
        if (this->ElementCounts[s1] < this->ElementCounts[s2])
        {
            this->Parents[s1] = s2;
            this->ElementCounts[s2] += this->ElementCounts[s1];
        }
        else
        {
            this->Parents[s2] = s1;
            this->ElementCounts[s1] += this->ElementCounts[s2];
        }

        this->TotalSetCount--;
    }
};


///////////////////////////////////////////////////////////////////////////////////////////


//
// Reference: http://www.geeksforgeeks.org/all-topological-sorts-of-a-directed-acyclic-graph/
// 
class DAG
{
public:

    int V;
    vector< vector<int> > AdjList;
    vector<int> Indegree;
    vector<int> Outdegree;

    vector< vector<int> > TopologicalSorts;
    vector<bool> Visited;

    void Init(int v)
    {
        V = v;

        AdjList.assign(V, vector<int>());

        Indegree.assign(V, 0);
        Outdegree.assign(V, 0);
    }

    void AddEdge(int from, int to)
    {
        AdjList[from].push_back(to);

        Outdegree[from]++;
        Indegree[to]++;
    }

    void AllTopologicalSorts()
    {
        TopologicalSorts.clear();

        Visited.assign(V, false);

        vector<int> res;
        AllTopologicalSorts(res);
    }

    void AllTopologicalSorts(vector<int>& res)
    {
        // To indicate whether all topological are found or not
        bool flag = false;

        for (int v = 0; v < V; v++)
        {
            //  If indegree is 0 and not yet visited then only choose that vertex
            if (Indegree[v] == 0 && !Visited[v])
            {
                //  reducing indegree of adjacent vertices
                for (int i = 0; i < (int)AdjList[v].size(); i++) {
                    Indegree[AdjList[v][i]]--;
                }

                //  including in result
                res.push_back(v);

                Visited[v] = true;

                AllTopologicalSorts(res);

                // resetting visited, res and indegree for backtracking
                Visited[v] = false;
                res.erase(res.end() - 1);

                for (int i = 0; i < (int)AdjList[v].size(); i++) {
                    Indegree[AdjList[v][i]]++;
                }

                flag = true;
            }
        }

        //  We reach here if all vertices are visited.
        if (!flag)
        {
            TopologicalSorts.push_back(res);

            //for (int i = 0; i < (int)res.size(); i++) {
            //	cout << res[i] << " ";
            //}
            //cout << endl;
        }
    }
};


///////////////////////////////////////////////////////////////////////////////////////////


typedef enum
{
    GRAPH_UNVISITED,
    GRAPH_VISITED,
    GRAPH_NO_ACCESS,
    GRAPH_EXPLORED,

} GRAPH_STATE;


class GraphMatrixCell
{
public:

    bool connected;
    double w;

    GraphMatrixCell(bool c = false, double ww = 0) {
        connected = c;
        w = ww;
    }

    bool operator < (const GraphMatrixCell& other) const {
        return w < other.w;
    }
};


struct GraphVertex 
{
    int v;
    double w;

    GraphVertex(int vv = 0, double ww = 0) {
        v = vv;
        w = ww;
    }

    bool operator < (const GraphVertex& other) const {
        return w < other.w;
    }

    static bool IsIndexSmaller(GraphVertex x1, GraphVertex x2)
    {
        return x1.v < x2.v;
    }
};


struct GraphEdge
{
    int v1;
    int v2;
    double w;

    GraphEdge(int vv1 = 0, int vv2 = 0, double ww = 0) {
        v1 = vv1;
        v2 = vv2;
        w = ww;
    }

    bool operator < (const GraphEdge& other) const {
        return w < other.w;
    }
};


struct GraphTraversal
{
    int state;
    int source_v;
    double accum_weight;

    GraphTraversal(int s = GRAPH_UNVISITED, int src = -1, double aw = 0) {
        state = s;
        source_v = src;
        accum_weight = aw;
    }

    bool operator < (const GraphTraversal& other) const {
        return accum_weight < other.accum_weight;
    }
};


class GraphPath
{
public:

    int source_v;
    vector<GraphEdge> edges;
    double weight_sum; // sum of the weights

    GraphPath(int source = 0)
    {
        Clear();
    }

    GraphPath(const GraphPath& other)
    {
        *this = other;
    }

    GraphPath& operator = (const GraphPath& other)
    {
        if (&other != this)
        {
            this->source_v = other.source_v;
            this->weight_sum = other.weight_sum;
            this->edges = other.edges;
        }

        return *this;
    }

    void Clear()
    {
        this->source_v = 0;
        this->weight_sum = 0;
        this->edges.clear();
    }

    bool IsEmpty() const
    {
        return this->edges.size() == 0;
    }

    int GetEdgeCount() const
    {
        return (int)this->edges.size();
    }

    void GetVertices(vector<int>& vertices) const
    {
        vertices.clear();

        for (int i = 0; i < (int)this->edges.size(); i++)
        {
            vertices.push_back(this->edges[i].v1);
        }
    }

    bool operator == (const GraphPath& other) const
    {
        return this->source_v == other.source_v;
    }

    bool operator < (const GraphPath& other) const
    {
        return this->source_v < other.source_v;
    }

};


class GraphEdgeInfo
{
public:

    int index;
    double weight;
    int v[2];				// resolved vertex indices for two ends
    set<int> edges[2];		// connneted edge indices for two ends

    bool operator < (const GraphEdgeInfo& other) const {
        return index < other.index;
    }

    int FindConnectedEndpoint(int edgeIdx)
    {
        for (int ep = 0; ep < 2; ep++) {
            if (edges[ep].find(edgeIdx) != edges[ep].end()) {
                return ep;
            }
        }
        return -1;
    }
};


///////////////////////////////////////////////////////////////////////////////////////////


class Graph
{
public:

    int V, E;

    //
    //  AdjacencyList[v1][j] means the j-th vertex adjacent to vertex v1
    //  v2 = AdjacencyList[v1][k].v, this means vertex v1 is adjacent to vertex v2
    //
    vector< vector< GraphVertex > > AdjacencyList;  

    //
    //  AdjacencyListMap[ (v1,v2) ] gives k where AdjacencyList[v1][k].v is v2
    //  if we frequently need to know the property of edge v1->v2 such as weight, we need this map to speed up the lookup in AdjacencyList[v1]
    //
    map< pair<int, int>, int> AdjacencyListMap;

    vector< vector< GraphMatrixCell > > AdjacencyMatrix;

public:

    Graph()
    {

        Clear();
    }

    Graph(const Graph& other)
    {
        *this = other;
    }

    Graph& operator=(const Graph& other)
    {
        if (&other != this)
        {
            this->V = other.V;
            this->E = other.E;

            this->AdjacencyList = other.AdjacencyList;
            this->AdjacencyListMap = other.AdjacencyListMap;
            this->AdjacencyMatrix = other.AdjacencyMatrix;
        }

        return *this;
    }

	void Clear()
	{
		this->V = 0;
		this->E = 0;

		this->AdjacencyList.clear();
		this->AdjacencyListMap.clear();

		this->AdjacencyMatrix.clear();
	}

	void CreateAdjacencyList(int v)
	{
		V = v;
		AdjacencyList.assign(V, vector<GraphVertex>());
	}

	void CreateReversedGraph(Graph& graph2)
	{
		graph2.V = this->V;
		graph2.E = this->E;

		graph2.AdjacencyList.assign(graph2.V, vector< GraphVertex >());

		for (int i = 0; i < this->V; i++)
		{
			vector< GraphVertex >& list = this->AdjacencyList[i];

			for (int j = 0; j < (int)list.size(); j++) {
				graph2.AdjacencyList[list[j].v].push_back(GraphVertex(i, list[j].w));
			}
		}
	}

    void AddEdgeToAdjacencyList(bool isDirected, int from, int to, double weight = 0)
    {
        this->AdjacencyList[from].push_back(GraphVertex(to, weight));

        if (!isDirected) {
            this->AdjacencyList[to].push_back(GraphVertex(from, weight));
        }

		E++;
    }

    void BuildAdjacencyListFromEdgeCollection(bool isDirected, vector<GraphEdgeInfo>& edge_collection)
    {
        V = 0;
        E = (int)edge_collection.size();
        AdjacencyList.clear();

        for (int i = 0; i < E; i++) {
            edge_collection[i].index = i;
            edge_collection[i].v[0] = edge_collection[i].v[1] = -1;
        }

        for (int current_edge = 0; current_edge < E; current_edge++)
        {
            GraphEdgeInfo& ei = edge_collection[current_edge];

            for (int ep = 0; ep < 2; ep++)
            {
                if (ei.v[ep] < 0)
                {
                    set<int>& connected_edges = ei.edges[ep];

                    for (auto iter = connected_edges.begin(); iter != connected_edges.end(); iter++)
                    {
                        int connected_edge = *iter;

                        int endpoint = edge_collection[connected_edge].FindConnectedEndpoint(current_edge);

                        if (endpoint >= 0 && edge_collection[connected_edge].v[endpoint] >= 0)
                        {
                            ei.v[ep] = edge_collection[connected_edge].v[endpoint];
                            break;
                        }
                    }
                }

                if (ei.v[ep] < 0)
                {
                    AdjacencyList.push_back(vector<GraphVertex>());
                    V++;
                    ei.v[ep] = V - 1;
                }
            }

            AddEdgeToAdjacencyList(isDirected, ei.v[0], ei.v[1], ei.weight);
        }
    }

    void BuildEdgeCollectionFromAdjacencyList(bool isDirected, vector<GraphEdgeInfo>& edge_collection)
    {
        edge_collection.clear();

        for (int v1 = 0; v1 < (int)this->AdjacencyList.size(); v1++)
        {
            for (int k = 0; k < (int)this->AdjacencyList[v1].size(); k++)
            {
                GraphVertex& item = this->AdjacencyList[v1][k];

                int v2 = item.v;

                int idx;
                for (idx = 0; idx < (int)edge_collection.size(); idx++)
                {
                    GraphEdgeInfo& ei = edge_collection[idx];

                    if ( (ei.v[0] == v1 && ei.v[1] == v2) || (!isDirected && ei.v[0] == v2 && ei.v[1] == v1)) {
                        break;
                    }
                }

                if (idx >= (int)edge_collection.size()) 
                {
                    GraphEdgeInfo ei;
                    ei.v[0] = v1; ei.v[1] = v2; ei.weight = item.w; ei.index = idx;

                    edge_collection.push_back(ei);
                }
            }
        }

        for (int i = 0; i < (int)edge_collection.size(); i++)
        {
            for (int j = 0; j < (int)edge_collection.size(); j++)
            {
                if (i == j) {
                    continue;
                }

                GraphEdgeInfo& ei = edge_collection[i];
                GraphEdgeInfo& ej = edge_collection[j];

                if (ei.v[0] == ej.v[0] || ei.v[0] == ej.v[1]) {
                    ei.edges[0].insert(j);
                }

                if (ei.v[1] == ej.v[0] || ei.v[1] == ej.v[1]) {
                    ei.edges[1].insert(j);
                }
            }
        }
    }

    void SetEdgeToAdjacencyMatrix(bool isDirected, int from, int to, double weight = 0)
    {
        this->AdjacencyMatrix[from][to].w = weight;

        if (!isDirected)
        {
            this->AdjacencyMatrix[to][from].w = weight;
        }
    }

    void BuildAdjacencyListMap()
    {
        this->AdjacencyListMap.clear();

        for (int v1 = 0; v1 < this->V; v1++)
        {
            for (int k = 0; k < (int)this->AdjacencyList[v1].size(); k++)
            {
                int v2 = this->AdjacencyList[v1][k].v;

                this->AdjacencyListMap[pair<int, int>(v1, v2)] = k;
            }
        }

        //
        // int k = AdjacencyListMap[ pair<int, int>(v1, v2) ];
        //
        // double weight  = AdjacencyList[v1][k].Weight;
        //
    }

    void BuildAdjacencyMatrixFromAdjacencyList(bool isDirected)
    {
        this->AdjacencyMatrix.clear();

        vector< GraphMatrixCell > vec;
        vec.assign(this->V, GraphMatrixCell());

        this->AdjacencyMatrix.assign(this->V, vec);

        for (int v1 = 0; v1 < this->V; v1++)
        {
            for (int k = 0; k < (int)this->AdjacencyList[v1].size(); k++)
            {
                GraphVertex& item = this->AdjacencyList[v1][k];

                int v2 = item.v;

                this->AdjacencyMatrix[v1][v2].connected = true;
                this->AdjacencyMatrix[v1][v2].w = item.w;

                if (!isDirected)
                {
                    this->AdjacencyMatrix[v2][v1].connected = true;
                    this->AdjacencyMatrix[v2][v1].w = item.w;
                }
            }
        }
    }

    void BuildAdjacencyListFromAdjacencyMatrix()
    {
        this->AdjacencyList.assign(this->V, vector< GraphVertex >() );

        for (int v1 = 0; v1 < this->V; v1++)
        {
            for (int v2 = 0; v2 < this->V; v2++)
            {
                if (this->AdjacencyMatrix[v1][v2].connected)
                {
                    this->AdjacencyList[v1].push_back(GraphVertex(v2, this->AdjacencyMatrix[v1][v2].w));
                }
            }
        }
    }
};


///////////////////////////////////////////////////////////////////////////////////////////


class GraphHelper
{
public:

	FILE* pOutputFile;

public:

	GraphHelper() {
		pOutputFile = 0;
	}

	static void InputWeightedGraph(Graph& graph, bool isDirected, const char* pInputFileName)
	{
		graph.Clear();

		FILE* pInputFile = 0;

		if (pInputFileName != 0) {
			pInputFile = freopen(pInputFileName, "r", stdin);
		}

		int E;
		scanf("%d %d", &graph.V, &E);

		graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());  

		int v1, v2;
		double w;
		for (int i = 0; i < E; i++)
		{
			scanf("%d %d %lf", &v1, &v2, &w);

			graph.AddEdgeToAdjacencyList(isDirected, v1, v2, w);
		}

		graph.E = E;

		if (pInputFile != 0)
		{
			fclose(pInputFile);
			pInputFile = 0;
		}
	}

	static void InputGraph(Graph& graph, bool isDirected, const char* pInputFileName)
	{
		graph.Clear();

		FILE* pInputFile = 0;

		if (pInputFileName != 0) {
			pInputFile = freopen(pInputFileName, "r", stdin);
		}

		int E;
		scanf("%d %d", &graph.V, &E);

		graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

		int v1, v2;
		for (int i = 0; i < E; i++)
		{
			scanf("%d %d", &v1, &v2);

			graph.AddEdgeToAdjacencyList(isDirected, v1, v2);
		}

		graph.E = E;

		if (pInputFile != 0)
		{
			fclose(pInputFile);
			pInputFile = 0;
		}
	}

	void OpenOutput(const char* pOutputFileName)
	{
		if (pOutputFileName != 0) {
			pOutputFile = freopen(pOutputFileName, "w", stdout);
		}
	}

	void CloseOutput()
	{
		fflush(stdout);

		if (pOutputFile != 0)
		{
			fclose(pOutputFile);
			pOutputFile = 0;
		}
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class GraphTraversalHandler
{
public:

    int V;
    int E;

    bool IsBipartite;
    bool HasCycle;
    bool IsDAG;

    vector<int> States;			// "States[i]" is the state of vertex "i"
    vector<int> Parents;		// "Parents[i]" is the parent of vertex "i", value -1 means no parent
    vector<bool> Colors;		// "Colors[i]" is the color of vertex "i"

    vector<int> Layers;			// "Layers[i]" is the layer of vertex "i" in the sense of BFS

	int CCGCount;               // Total number of connected component groups
	vector<int> CCG;			// "CCG[i]" is the connected component group number of vertex "i" 

	int source_v;
    vector<int> VertexQueue;	// "VertexQueue[i]" is the vertex at i-th iteration order for DFS & BFS, also used for TopologicalSort, etc.

public:

    GraphTraversalHandler()
    {
        Clear();
    }

    void Clear()
    {
        this->V = 0;
        this->E = 0;

        this->CCGCount = 0;

        this->IsBipartite = false;
        this->HasCycle = false;
        this->IsDAG = false;

        this->States.clear();
        this->Parents.clear();
        this->Colors.clear();
        this->CCG.clear();

        this->Layers.clear();

        this->source_v = -1;
        this->VertexQueue.clear();
    }

    void DfsComprehensive(Graph& graph, int currentVertex, int iterCounter = 0)
    {

        if (this->source_v < 0) // indicates the first call of the recursive chain
        {
            iterCounter = 0;

            this->source_v = currentVertex;

            this->V = graph.V;
            this->E = graph.E;

            this->States.assign(this->V, GRAPH_UNVISITED);
            this->Parents.assign(this->V, -1);
            this->Colors.assign(this->V, false);
            this->CCG.assign(this->V, 0);

            this->IsBipartite = true;
        }

        this->States[currentVertex] = GRAPH_VISITED;
        this->CCG[currentVertex] = this->CCGCount;

        for (int j = 0; j < (int)graph.AdjacencyList[currentVertex].size(); j++)
        {
            int childVertex = graph.AdjacencyList[currentVertex][j].v;

            if (this->States[childVertex] == GRAPH_UNVISITED)
            {
                this->Parents[childVertex] = currentVertex;
                this->Colors[childVertex] = !this->Colors[currentVertex];

                DfsComprehensive(graph, childVertex, ++iterCounter);
            }
            else
            {
                this->HasCycle = true;

                if (this->Colors[childVertex] == this->Colors[currentVertex])
                {
                    this->IsBipartite = false;
                }
            }
        }

        this->VertexQueue.push_back(currentVertex);
    }

	void BfsComprehensive(Graph& graph, int source)
	{
		this->source_v = source;

		this->V = graph.V;
		this->E = graph.E;

		this->States.assign(this->V, GRAPH_UNVISITED);
		this->Parents.assign(this->V, -1);
		this->Colors.assign(this->V, false);
		this->CCG.assign(this->V, 0);

		this->Layers.assign(this->V, 0);

		this->IsBipartite = true;

		this->VertexQueue.push_back(source);
		this->States[0] = GRAPH_EXPLORED;

		int IterIndex = 0;

		while (IterIndex < (int)this->VertexQueue.size())
		{
			int currentVertex = this->VertexQueue[IterIndex];

			this->States[currentVertex] = GRAPH_VISITED;
			this->CCG[currentVertex] = this->CCGCount;

			for (int j = 0; j < (int)graph.AdjacencyList[currentVertex].size(); j++)
			{
				int childVertex = graph.AdjacencyList[currentVertex][j].v;

				if (this->States[childVertex] == GRAPH_UNVISITED)
				{
					this->States[childVertex] = GRAPH_EXPLORED;
					this->Layers[childVertex] = this->Layers[currentVertex] + 1;
					this->Parents[childVertex] = currentVertex;
					this->Colors[childVertex] = !this->Colors[currentVertex];

					this->VertexQueue.push_back(childVertex);
				}
				else
				{
					this->HasCycle = true;

					if (this->Colors[childVertex] == this->Colors[currentVertex])
					{
						this->IsBipartite = false;
					}
				}
			}

			IterIndex++;
		}
	}

	void ScanAllByDfsComprehensive(Graph& graph)
	{
		Clear();

		//
		// This should be part of the loop below,
		// but our implementation relys on the first-time DfsComprehensive call to allocate memory,
		// hence we single it out here. 
		//
		DfsComprehensive(graph, 0);

		this->CCGCount = 1;

		for (int i = 1; i < this->V; i++)
		{
			if (this->States[i] == GRAPH_UNVISITED)
			{
				this->CCGCount++;
				DfsComprehensive(graph, i);
			}
		}

		if (!this->HasCycle && this->V == (int)this->VertexQueue.size())
		{
			this->IsDAG = true;
		}
	}

	static void DfsSimpleStack(Graph& graph, vector<int>& vertexState, int startVertex, vector<int>* pProcessOrder = 0)
    {
        if (pProcessOrder != 0) {
			pProcessOrder->clear();
        }

        if (vertexState[startVertex] != GRAPH_UNVISITED) {
            return;
        }

        stack<int> s;

        s.push(startVertex);
        vertexState[startVertex] = GRAPH_VISITED;

        while (!s.empty())
        {
            int currentVertex = s.top();
            s.pop();

            //
            // do something on the vertex
            //

			if (pProcessOrder != 0) {
				pProcessOrder->push_back(currentVertex);
			}

            vector<GraphVertex>& adjs = graph.AdjacencyList[currentVertex];
            int size = (int)adjs.size();

            for (int j = 0; j < size; j++)
            {
                int nextVertex = adjs[j].v;

                if (vertexState[nextVertex] == GRAPH_UNVISITED)
                {
                    s.push(nextVertex);
                    vertexState[nextVertex] = GRAPH_VISITED;
                }
            }
        }
    }

	static void ScanAllByDfsSimpleStack(Graph& graph, vector<int>& vertexState)
	{
		for (int i = 0; i < graph.V; i++) {
			if (vertexState[i] == GRAPH_UNVISITED) {
				DfsSimpleStack(graph, vertexState, i);
			}
		}
	}

	static void TopologicalSortByDfsSimpleStack(Graph& graph, vector<int>& vertexState, vector<int>& topology)
	{
		topology.clear();

		vector<int> procOrder;

		for (int i = 0; i < graph.V; i++) {
			if (vertexState[i] == GRAPH_UNVISITED)
			{
				procOrder.clear();
				DfsSimpleStack(graph, vertexState, i, &procOrder);

				topology.insert(topology.begin(), procOrder.begin(), procOrder.end());
			}
		}
	}

	static void DfsSimpleRecursive(Graph& graph, vector<int>& vertexState, int currentVertex, vector<int>* pProcessOrder = 0)
	{
		if (vertexState[currentVertex] != GRAPH_UNVISITED) {
			return;
		}

		vertexState[currentVertex] = GRAPH_VISITED;

		//
		// do something on the vertex
		//

		vector<GraphVertex>& adjs = graph.AdjacencyList[currentVertex];
		int size = (int)adjs.size();

		for (int j = 0; j < size; j++)
		{
			int nextVertex = adjs[j].v;

			if (vertexState[nextVertex] == GRAPH_UNVISITED) {
				DfsSimpleRecursive(graph, vertexState, nextVertex, pProcessOrder);
			}
		}

		if (pProcessOrder != 0) {
			pProcessOrder->push_back(currentVertex);
		}
	}

	static void ScanAllByDfsSimpleRecursive(Graph& graph, vector<int>& vertexState)
	{
		for (int i = 0; i < graph.V; i++) {
			if (vertexState[i] == GRAPH_UNVISITED) {
				DfsSimpleRecursive(graph, vertexState, i);
			}
		}
	}

	static void TopologicalSortByDfsSimpleRecursive(Graph& graph, vector<int>& vertexState, vector<int>& topology)
	{
		topology.clear();

		for (int i = 0; i < graph.V; i++) {
			if (vertexState[i] == GRAPH_UNVISITED) {
				DfsSimpleRecursive(graph, vertexState, i, &topology);
			}
		}

		reverse(topology.begin(), topology.end());
	}

	static void BfsSimpleQueue(Graph& graph, vector<int>& vertexState, int startVertex, vector<int>* pProcessOrder = 0)
	{
		if (pProcessOrder != 0) {
			pProcessOrder->clear();
		}

		if (vertexState[startVertex] != GRAPH_UNVISITED) {
			return;
		}

		queue<int> q;

		q.push(startVertex);
		vertexState[startVertex] = GRAPH_VISITED;

		while (!q.empty())
		{
			int currentVertex = q.front();
			q.pop();

			//
			// do something on the vertex
			//

			if (pProcessOrder != 0) {
				pProcessOrder->push_back(currentVertex);
			}

			vector<GraphVertex>& adjs = graph.AdjacencyList[currentVertex];
			int size = (int)adjs.size();

			for (int j = 0; j < size; j++)
			{
				int nextVertex = adjs[j].v;

				if (vertexState[nextVertex] == GRAPH_UNVISITED)
				{
					q.push(nextVertex);
					vertexState[nextVertex] = GRAPH_VISITED;
				}
			}
		}
	}
	
	static void GetLongestStepsByDfs(Graph& graph, vector<int>& vertexState, int s, int t, int steps, int& max_steps)
	{
		if (vertexState[s] != GRAPH_UNVISITED) {
			return;
		}

		if (s == t)
		{
			if (steps > max_steps) {
				max_steps = steps;
			}

			return;
		}

		vertexState[s] = GRAPH_VISITED;

		vector<GraphVertex>& adjs = graph.AdjacencyList[s];
		int size = (int)adjs.size();

		for (int j = 0; j < size; j++)
		{
			int nextVertex = adjs[j].v;

			if (vertexState[nextVertex] == GRAPH_UNVISITED) {
				GetLongestStepsByDfs(graph, vertexState, nextVertex, t, steps+1, max_steps);
			}
		}

		vertexState[s] = GRAPH_UNVISITED;
	}

	static pair<int, int> GetLongestStepsByDfs(Graph& graph, vector<int>& vertexState, int s)
	{
		pair<int, int> t(s, 0);

		if (vertexState[s] != GRAPH_UNVISITED) {
			return t;
		}

		vector<int> steps(graph.V, 0);

		stack<int> stk;

		stk.push(s);
		vertexState[s] = GRAPH_VISITED;

		while (!stk.empty())
		{
			int currentVertex = stk.top();
			stk.pop();

			vector<GraphVertex>& adjs = graph.AdjacencyList[currentVertex];
			int size = (int)adjs.size();

			for (int j = 0; j < size; j++)
			{
				int nextVertex = adjs[j].v;

				if (vertexState[nextVertex] == GRAPH_UNVISITED)
				{
					stk.push(nextVertex);
					vertexState[nextVertex] = GRAPH_VISITED;

					steps[nextVertex] = steps[currentVertex] + 1;
				}
			}

			vertexState[currentVertex] = GRAPH_VISITED;
		}

		for (int i = 0; i < graph.V; i++)
		{
			if (steps[i] > t.second)
			{
				t.second = steps[i];
				t.first = i;
			}
		}

		return t;
	}

	static pair<int, int> GetLongestStepsByBfs(Graph& graph, vector<int>& vertexState, int s)
	{
		pair<int, int> t(s, 0);

		if (vertexState[s] != GRAPH_UNVISITED) {
			return t;
		}

		vector<int> steps(graph.V, 0);

		queue<int> q;

		q.push(s);
		vertexState[s] = GRAPH_VISITED;

		while (!q.empty())
		{
			int currentVertex = q.front();
			q.pop();

			vector<GraphVertex>& adjs = graph.AdjacencyList[currentVertex];
			int size = (int)adjs.size();

			for (int j = 0; j < size; j++)
			{
				int nextVertex = adjs[j].v;

				if (vertexState[nextVertex] == GRAPH_UNVISITED)
				{
					q.push(nextVertex);
					vertexState[nextVertex] = GRAPH_VISITED;

					steps[nextVertex] = steps[currentVertex] + 1;
				}
			}
		}

		for (int i = 0; i < graph.V; i++)
		{
			if (steps[i] > t.second)
			{
				t.second = steps[i];
				t.first = i;
			}
		}

		return t;
	}
	
	static void GetLongestStepsByBfs(Graph& graph, vector<int>& vertexState, int& s, int& t, int& max_steps)
	{
		s = t = 0;

		pair<int, int> p1 = GetLongestStepsByBfs(graph, vertexState, 0);

		int size = (int)vertexState.size();
		for (int i = 0; i < size; i++) {
			if (vertexState[i] == GRAPH_VISITED) {
				vertexState[i] = GRAPH_UNVISITED;
			}
		}

		pair<int, int> p2 = GetLongestStepsByBfs(graph, vertexState, p1.first);


		s = p2.first;
		t = p1.first;
		max_steps = p2.second;
	}

    //
    // Figure out the Strongly Connected Components. The edges are considered directed.
    //
	static int KosarajuAlgorithm(Graph& graph, vector<int>& sccGroup, vector<vector<int>>* pSCC = 0 )
	{
		Graph graph2;
		graph.CreateReversedGraph(graph2);

		vector<int> vertexState(graph2.V, GRAPH_UNVISITED);
		vector<int> topology;

		TopologicalSortByDfsSimpleStack(graph2, vertexState, topology);

		vertexState.assign(graph.V, GRAPH_UNVISITED);

		vector<int> procOrder;

		sccGroup.assign(graph.V, -1);

		if (pSCC != 0) {
			pSCC->clear();
		}

		int num = 0;

		for (int i = 0; i < graph.V; i++) 
		{
			if (vertexState[topology[i]] == GRAPH_UNVISITED)
			{
				if (pSCC != 0) {
					pSCC->push_back( vector<int>() );
				}

				vector<int>* pProcOrder = (pSCC == 0) ? &procOrder : &((*pSCC)[num]);

				DfsSimpleStack(graph, vertexState, topology[i], pProcOrder);

				int size = (int)pProcOrder->size();
				for (int k = 0; k < size; k++) {
					sccGroup[(*pProcOrder)[k]] = num;
				}

				num++;
			}
		}

		int num2 = num;
		for (int i = 0; i < graph.V; i++)
		{
			if (sccGroup[i] < 0) {
				sccGroup[i] = num2++;
			}
		}

		return num;
	}

	//
	// Taking SSC as a node, create a new graph
	//
	static int ContractSCC(Graph& graph, Graph& graph2)
	{
		vector<int> sccGroup;
		int num = KosarajuAlgorithm(graph, sccGroup);

		if (num == 1)
		{
			graph2.V = 1;
			graph2.E = 0;
			graph2.AdjacencyList.assign(1, vector<GraphVertex>());
			return num;
		}

		graph2.V = num;
		graph2.E = 0;
		graph2.AdjacencyList.assign(graph.V, vector<GraphVertex>());

		set<pair<int, int>> edges;

		for (int from = 0; from < graph.V; from++)
		{
			vector<GraphVertex>& vecAdj = graph.AdjacencyList[from];
			int size = (int)vecAdj.size();
			GraphVertex* adjs = (size == 0) ? 0 : &vecAdj.front();

			for (int j = 0; j < size; j++)
			{
				int to = adjs[j].v;

				pair<int, int> edge(sccGroup[from], sccGroup[to]);

				if (edge.first != edge.second && edges.find(edge) == edges.end())
				{
					graph2.AdjacencyList[edge.first].push_back(edge.second);
					
					edges.insert(edge);
					
					graph2.E++;
				}
			}
		}

		return num;
	}


public:

    static void Test()
    {
        Test_01();
        Test_02();
        Test_03();
        Test_04();
		Test_05();
		Test_06();
	}

    static void Test_01()
    {
        Graph graph;

        //
        // Ref: Algorithms, 4th Edition
        // from input_graph_cg_tiny.txt
        //
        graph.V = 6;
        graph.E = 8;

        graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

        graph.AddEdgeToAdjacencyList(true, 0, 2);
        graph.AddEdgeToAdjacencyList(true, 0, 1);
        graph.AddEdgeToAdjacencyList(true, 0, 5);
        graph.AddEdgeToAdjacencyList(true, 1, 0);
        graph.AddEdgeToAdjacencyList(true, 1, 2);
        graph.AddEdgeToAdjacencyList(true, 2, 0);
        graph.AddEdgeToAdjacencyList(true, 2, 1);
        graph.AddEdgeToAdjacencyList(true, 2, 3);
        graph.AddEdgeToAdjacencyList(true, 2, 4);
        graph.AddEdgeToAdjacencyList(true, 3, 5);
        graph.AddEdgeToAdjacencyList(true, 3, 4);
        graph.AddEdgeToAdjacencyList(true, 3, 2);
        graph.AddEdgeToAdjacencyList(true, 4, 3);
        graph.AddEdgeToAdjacencyList(true, 4, 2);
        graph.AddEdgeToAdjacencyList(true, 5, 3);
        graph.AddEdgeToAdjacencyList(true, 5, 0);


        GraphTraversalHandler gth;
        gth.DfsComprehensive(graph, 0);

        string str;

        for (int i = 0; i < (int)gth.VertexQueue.size(); i++) {
            str += '0' + gth.VertexQueue[i];
        }

        reverse(str.begin(), str.end());

        bool b = str == "023451";

        printf("GraphTraversalHandler::DfsComprehensive, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


        vector<int> vertexState(graph.V, GRAPH_UNVISITED);
        vector<int> processOrder;

        GraphTraversalHandler::DfsSimpleStack(graph, vertexState, 0, &processOrder);

        str = "";
        for (int i = 0; i < (int)processOrder.size(); i++) {
            str += '0' + processOrder[i];
        }

        b = str == "053412";

        printf("GraphTraversalHandler::DfsSimpleStack, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


		vertexState.assign(graph.V, GRAPH_UNVISITED);
		processOrder.clear();

		GraphTraversalHandler::DfsSimpleRecursive(graph, vertexState, 0, &processOrder);

		str = "";
		for (int i = 0; i < (int)processOrder.size(); i++) {
			str += '0' + processOrder[i];
		}

		b = str == "021354";

		printf("GraphTraversalHandler::DfsSimpleRecursive, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());

    }

    static void Test_02()
    {
		Graph graph;

		//
		// Ref: Algorithms, 4th Edition
		// from input_graph_cg_tiny.txt
		//
		graph.V = 6;
		graph.E = 8;

		graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

		graph.AddEdgeToAdjacencyList(true, 0, 2);
		graph.AddEdgeToAdjacencyList(true, 0, 1);
		graph.AddEdgeToAdjacencyList(true, 0, 5);
		graph.AddEdgeToAdjacencyList(true, 1, 0);
		graph.AddEdgeToAdjacencyList(true, 1, 2);
		graph.AddEdgeToAdjacencyList(true, 2, 0);
		graph.AddEdgeToAdjacencyList(true, 2, 1);
		graph.AddEdgeToAdjacencyList(true, 2, 3);
		graph.AddEdgeToAdjacencyList(true, 2, 4);
		graph.AddEdgeToAdjacencyList(true, 3, 5);
		graph.AddEdgeToAdjacencyList(true, 3, 4);
		graph.AddEdgeToAdjacencyList(true, 3, 2);
		graph.AddEdgeToAdjacencyList(true, 4, 3);
		graph.AddEdgeToAdjacencyList(true, 4, 2);
		graph.AddEdgeToAdjacencyList(true, 5, 3);
		graph.AddEdgeToAdjacencyList(true, 5, 0);


		GraphTraversalHandler gth;
		gth.BfsComprehensive(graph, 0);

		string str;

		for (int i = 0; i < (int)gth.VertexQueue.size(); i++) {
			str += '0' + gth.VertexQueue[i];
		}

		bool b = str == "021534";

		printf("GraphTraversalHandler::BfsComprehensive, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());


		vector<int> vertexState(graph.V, GRAPH_UNVISITED);
		vector<int> processOrder;

		GraphTraversalHandler::BfsSimpleQueue(graph, vertexState, 0, &processOrder);

		str = "";
		for (int i = 0; i < (int)processOrder.size(); i++) {
			str += '0' + processOrder[i];
		}

		b = str == "021534";

		printf("GraphTraversalHandler::BfsSimpleQueue, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());

	}

	static void Test_03()
	{
		Graph graph;

		//
		// Ref: http://www.geeksforgeeks.org/longest-path-undirected-tree/
		//
		graph.V = 10;
		graph.E = 9;

		graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

		graph.AddEdgeToAdjacencyList(false, 0, 1);
		graph.AddEdgeToAdjacencyList(false, 1, 2);
		graph.AddEdgeToAdjacencyList(false, 2, 3);
		graph.AddEdgeToAdjacencyList(false, 2, 9);
		graph.AddEdgeToAdjacencyList(false, 2, 4);
		graph.AddEdgeToAdjacencyList(false, 4, 5);
		graph.AddEdgeToAdjacencyList(false, 1, 6);
		graph.AddEdgeToAdjacencyList(false, 6, 7);
		graph.AddEdgeToAdjacencyList(false, 6, 8);

		vector<int> vertexState(graph.V, GRAPH_UNVISITED);
		int max_steps = 0;

		GraphTraversalHandler::GetLongestStepsByDfs(graph, vertexState, 5, 7, 0, max_steps);

		bool b = max_steps == 5;

		printf("GraphTraversalHandler::GetLongestStepsByDfs, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), max_steps);



		vertexState.assign(graph.V, GRAPH_UNVISITED);

		auto p = GraphTraversalHandler::GetLongestStepsByDfs(graph, vertexState, 5);

		b = p.first == 7 && p.second == 5;

		printf("GraphTraversalHandler::GetLongestStepsByDfs, %s, the answer is (%d, %d) \n", (b ? "Correct" : "Wrong"), p.first, p.second);



		vertexState.assign(graph.V, GRAPH_UNVISITED);

		p = GraphTraversalHandler::GetLongestStepsByBfs(graph, vertexState, 5);

		b = p.first == 7 && p.second == 5;

		printf("GraphTraversalHandler::GetLongestStepsByBfs, %s, the answer is (%d, %d) \n", (b ? "Correct" : "Wrong"), p.first, p.second);



		vertexState.assign(graph.V, GRAPH_UNVISITED);

		int s, t;
		GraphTraversalHandler::GetLongestStepsByBfs(graph, vertexState, s, t, max_steps);

		b = s == 7 && t == 5 && max_steps == 5;

		printf("GraphTraversalHandler::GetLongestStepsByBfs, %s, the answer is (%d, %d, %d) \n", (b ? "Correct" : "Wrong"), s, t, max_steps);
	}

	static void Test_04()
	{
		Graph graph;

		//
		// Ref: output_graph_toposort.txt
		//
		graph.V = 8;
		graph.E = 8;

		graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

		graph.AddEdgeToAdjacencyList(true, 0, 1);
		graph.AddEdgeToAdjacencyList(true, 0, 2);
		graph.AddEdgeToAdjacencyList(true, 1, 3);
		graph.AddEdgeToAdjacencyList(true, 1, 2);
		graph.AddEdgeToAdjacencyList(true, 2, 3);
		graph.AddEdgeToAdjacencyList(true, 2, 5);
		graph.AddEdgeToAdjacencyList(true, 3, 4);
		graph.AddEdgeToAdjacencyList(true, 7, 6);

		vector<int> answer = { 7, 6, 0, 1, 2, 5, 3, 4 };

		vector<int> vertexState(graph.V, GRAPH_UNVISITED);
		vector<int> topology;

		GraphTraversalHandler::TopologicalSortByDfsSimpleRecursive(graph, vertexState, topology);

		string str = "( ";
		for (int i = 0; i < (int)topology.size(); i++) {
			str += topology[i] + '0';
			str += ' ';
		}
		str += ")";

		bool b = topology == answer;

		printf("GraphTraversalHandler::TopologicalSortByDfsSimpleRecursive, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());
	}
	
	static void Test_05()
	{
		Graph graph;

		//
		// Ref: output_graph_toposort.txt
		//
		graph.V = 8;
		graph.E = 8;

		graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

		graph.AddEdgeToAdjacencyList(true, 0, 1);
		graph.AddEdgeToAdjacencyList(true, 0, 2);
		graph.AddEdgeToAdjacencyList(true, 1, 3);
		graph.AddEdgeToAdjacencyList(true, 1, 2);
		graph.AddEdgeToAdjacencyList(true, 2, 3);
		graph.AddEdgeToAdjacencyList(true, 2, 5);
		graph.AddEdgeToAdjacencyList(true, 3, 4);
		graph.AddEdgeToAdjacencyList(true, 7, 6);

		vector<int> vertexState(graph.V, GRAPH_UNVISITED);
		vector<int> topology;
		GraphTraversalHandler::TopologicalSortByDfsSimpleStack(graph, vertexState, topology);

		vector<int> answer = { 7, 6, 0, 2, 5, 3, 4, 1 };

		string str = "( ";
		for (int i = 0; i < (int)topology.size(); i++) {
			str += topology[i] + '0';
			str += ' ';
		}
		str += ")";

		bool b = topology == answer;

		printf("GraphTraversalHandler::TopologicalSortByDfsSimpleStack, %s, the answer is %s \n", (b ? "Correct" : "Wrong"), str.c_str());
	}

	static void Test_06()
	{
		Graph graph;

		//
		// Ref: http://www.geeksforgeeks.org/strongly-connected-components/
		//
		graph.V = 5;
		graph.E = 5;

		graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

		graph.AddEdgeToAdjacencyList(true, 0, 3);
		graph.AddEdgeToAdjacencyList(true, 3, 4);
		graph.AddEdgeToAdjacencyList(true, 1, 0);
		graph.AddEdgeToAdjacencyList(true, 0, 2);
		graph.AddEdgeToAdjacencyList(true, 2, 1);

		vector<int> sccGroup;
		vector<vector<int>> scc;
		int num = KosarajuAlgorithm(graph, sccGroup, &scc);
		bool b = num == 3;

		printf("GraphTraversalHandler::KosarajuAlgorithm, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), num);


		Graph graph2;
		ContractSCC(graph, graph2);
		b = graph2.V == 3;

		printf("GraphTraversalHandler::ContractSCC, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), num);


		//
		// Ref: input_graph_toposort.txt
		//
		graph.V = 8;
		graph.E = 8;

		graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

		graph.AddEdgeToAdjacencyList(true, 0, 1);
		graph.AddEdgeToAdjacencyList(true, 0, 2);
		graph.AddEdgeToAdjacencyList(true, 1, 3);
		graph.AddEdgeToAdjacencyList(true, 1, 2);
		graph.AddEdgeToAdjacencyList(true, 2, 3);
		graph.AddEdgeToAdjacencyList(true, 2, 5);
		graph.AddEdgeToAdjacencyList(true, 3, 4);
		graph.AddEdgeToAdjacencyList(true, 7, 6);

		num = KosarajuAlgorithm(graph, sccGroup, &scc);
		b = num == 8;

		printf("GraphTraversalHandler::KosarajuAlgorithm, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), num);
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class GraphArticulationHandler
{
public:

	vector<int> States;	
	vector<int> Parents;	// -1 means no parent

	vector<int> IterationDistance;		
	vector<int> MinIterationDistance;	
	vector<bool> IsAP;  
	set<int> ArticulationPoints;

	int IterationCounter;
	int RootVertex;
	int RootChildrenCount;


public:

	void ArticulationPointByDfsRecursive(Graph& graph, int currentVertex)
	{
		States[currentVertex] = GRAPH_VISITED;
		IterationDistance[currentVertex] = MinIterationDistance[currentVertex] = IterationCounter++;

		for (int j = 0; j < (int)graph.AdjacencyList[currentVertex].size(); j++)
		{
			int childVertex = graph.AdjacencyList[currentVertex][j].v;

			if (States[childVertex] == GRAPH_UNVISITED)
			{
				this->Parents[childVertex] = currentVertex;

				if (currentVertex == RootVertex) {
					RootChildrenCount++;
				}

				ArticulationPointByDfsRecursive(graph, childVertex);

				MinIterationDistance[currentVertex] = min(MinIterationDistance[currentVertex], MinIterationDistance[childVertex]);

				if (Parents[currentVertex] != -1 && MinIterationDistance[childVertex] >= IterationDistance[currentVertex]) 
				{
					IsAP[currentVertex] = true;
					ArticulationPoints.insert(currentVertex);
				} 

				//if (MinIterationDistance[childVertex] > this->IterationDistance[currentVertex])
				//{
				//	// Edge childVertex-currentVertex is a bridge
				//}
			}
			else
			{
				if (childVertex != Parents[currentVertex])  {	// a back edge and not direct cycle
					MinIterationDistance[currentVertex] = min(MinIterationDistance[currentVertex], IterationDistance[childVertex]);
				}
			}
		}
	}

	void FindArticulationPoints(Graph& graph, int startVertex, set<int>* pAP = 0, set<pair<int, int>>* pBridges = 0)
	{
		States.assign(graph.V, GRAPH_UNVISITED);
		Parents.assign(graph.V, -1);

		IterationDistance.assign(graph.V, 0);
		MinIterationDistance.assign(graph.V, 0);
		IsAP.assign(graph.V, false);
		ArticulationPoints.clear();

		IterationCounter = 0;


		for (int v = 0; v < graph.V; v++)
		{
			if (States[v] == GRAPH_UNVISITED)
			{
				RootVertex = v;
				RootChildrenCount = 0;

				ArticulationPointByDfsRecursive(graph, v);

				IsAP[RootVertex] = RootChildrenCount > 1;
				if (IsAP[RootVertex]) {
					ArticulationPoints.insert(RootVertex);
				}
			}
		}
		
		if (pAP != 0)
		{
			*pAP = ArticulationPoints;
		}

		if (pBridges == 0) {
			return;
		}

		pBridges->clear();

		for (int v1 = 0; v1 < graph.V; v1++) {
			if (IsAP[v1]) {
				for (int j = 0; j < (int)graph.AdjacencyList[v1].size(); j++)
				{
					int v2 = graph.AdjacencyList[v1][j].v;

					//
					// Assume we are dealing with undirected graph
					//
					pair<int, int> e(min(v1, v2), max(v1, v2));

					pBridges->insert(e);
				}
			}
		}
	}

	// INCORRECT
	/*
	static void ArticulationPointsByDfsStack(Graph& graph, int rootVertex, set<int>* pAP = 0, set<pair<int, int>>* pBridges = 0)
	{
		vector<int> states(graph.V, GRAPH_UNVISITED);
		vector<int> parents(graph.V, -1);
		vector<int> dist(graph.V, 0);
		vector<int> minDist(graph.V, 0);
		vector<bool> isAP(graph.V, false);
		int iterationCount = -1;
		int rootChildrenCount = 0;

		stack<int> s;

		s.push(rootVertex);

		states[rootVertex] = GRAPH_EXPLORED;
		dist[rootVertex] = minDist[rootVertex] = ++iterationCount;

		while (!s.empty())
		{
			int currentVertex = s.top();
			s.pop();

			vector<GraphVertex>& adjs = graph.AdjacencyList[currentVertex];
			int size = (int)adjs.size();

			for (int j = 0; j < size; j++)
			{
				int nextVertex = adjs[j].v;

				if (states[nextVertex] == GRAPH_UNVISITED)
				{
					if (currentVertex == rootVertex) {
						rootChildrenCount++;
					}

					parents[nextVertex] = currentVertex;

					s.push(nextVertex);

					states[nextVertex] = GRAPH_EXPLORED;
				}
			}
		}

		s.push(rootVertex);

		states[rootVertex] = GRAPH_VISITED;

		while (!s.empty())
		{
			int currentVertex = s.top();
			s.pop();

			vector<GraphVertex>& adjs = graph.AdjacencyList[currentVertex];
			int size = (int)adjs.size();

			for (int j = 0; j < size; j++)
			{
				int nextVertex = adjs[j].v;

				if (states[nextVertex] == GRAPH_EXPLORED)
				{
					dist[nextVertex] = minDist[nextVertex] = ++iterationCount;
					minDist[currentVertex] = min(minDist[currentVertex], minDist[nextVertex]);

					s.push(nextVertex);

					states[nextVertex] = GRAPH_VISITED;
				}
				else if (states[nextVertex] == GRAPH_VISITED)
				{
					if (nextVertex != parents[currentVertex]) 		// a back edge and not direct cycle
					{
						minDist[currentVertex] = min(minDist[currentVertex], dist[nextVertex]);
					}
				}
			}
		}

		for (int v1 = 0; v1 < graph.V; v1++) {
			for (int j = 0; j < (int)graph.AdjacencyList[v1].size(); j++)
			{
				int v2 = graph.AdjacencyList[v1][j].v;

				if (minDist[v2] >= dist[v1]) {
					isAP[v1] = true;
				}
			}
		}

		isAP[rootVertex] = rootChildrenCount > 1;

		if (pAP != 0)
		{
			pAP->clear();

			for (int i = 0; i < (int)isAP.size(); i++) {
				if (isAP[i]) {
					pAP->insert(i);
				}
			}
		}

		if (pBridges == 0) {
			return;
		}

		pBridges->clear();

		for (int v1 = 0; v1 < graph.V; v1++) {
			if (isAP[v1]) {
				for (int j = 0; j < (int)graph.AdjacencyList[v1].size(); j++)
				{
					int v2 = graph.AdjacencyList[v1][j].v;

					pair<int, int> e(min(v1, v2), max(v1, v2));  // Assume we are dealing with undirected graph

					pBridges->insert(e);
				}
			}
		}
	}
	*/

public:

    static void Test()
    {
        Test_01();
		Test_02();
		Test_03();
	}

	static void Test_01()
	{
		Graph graph;

		//
		// Ref: input_graph_articulation_point.txt
		//
		graph.V = 6;
		graph.E = 5;

		graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

		graph.AddEdgeToAdjacencyList(false, 0, 1);
		graph.AddEdgeToAdjacencyList(false, 1, 2);
		graph.AddEdgeToAdjacencyList(false, 1, 4);
		graph.AddEdgeToAdjacencyList(false, 3, 4);
		graph.AddEdgeToAdjacencyList(false, 4, 5);

		set<int> APs;
		set<pair<int, int>> bridges;

		GraphArticulationHandler gah;
		gah.FindArticulationPoints(graph, 0, &APs, &bridges);

		set<int> ap_answer = { 1, 4 };
		set<pair<int, int>> bridge_answer = { { 0, 1 }, {1, 2}, {1, 4}, {3, 4}, {4, 5} };

		bool b = APs == ap_answer && bridges == bridge_answer;

		string str = "Articulation Points ( ";
		for (auto i = APs.begin(); i != APs.end(); i++) {
			str += *i + '0';
			str += ' ';
		}
		str += ")";

		string str2 = "Bridges ( ";
		for (auto i = bridges.begin(); i != bridges.end(); i++) {
			str2 += '[';  str2 += (i->first + '0'); str2 += " -> "; str2 += (i->second + '0'); str2 += ']';
			str2 += ' ';
		}
		str2 += ")";

		printf("GraphArticulationHandler::FindArticulationPoints, %s, the answer is %s, %s \n", (b ? "Correct" : "Wrong"), str.c_str(), str2.c_str());
	}

	static void Test_02()
	{
		Graph graph;

		//
		// Ref: http://www.geeksforgeeks.org/articulation-points-or-cut-vertices-in-a-graph/
		//
		graph.V = 5;
		graph.E = 5;

		graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

		graph.AddEdgeToAdjacencyList(false, 0, 1);
		graph.AddEdgeToAdjacencyList(false, 0, 2);
		graph.AddEdgeToAdjacencyList(false, 0, 3);
		graph.AddEdgeToAdjacencyList(false, 1, 2);
		graph.AddEdgeToAdjacencyList(false, 3, 4);

		set<int> APs;
		set<pair<int, int>> bridges;

		GraphArticulationHandler gah;
		gah.FindArticulationPoints(graph, 0, &APs, &bridges);

		set<int> ap_answer = { 0, 3 };
		set<pair<int, int>> bridge_answer = { { 0, 1 }, { 0, 2 }, { 0, 3 }, {3, 4} };

		bool b = APs == ap_answer && bridges == bridge_answer;

		string str = "Articulation Points ( ";
		for (auto i = APs.begin(); i != APs.end(); i++) {
			str += *i + '0';
			str += ' ';
		}
		str += ")";

		string str2 = "Bridges ( ";
		for (auto i = bridges.begin(); i != bridges.end(); i++) {
			str2 += '[';  str2 += (i->first + '0'); str2 += " -> "; str2 += (i->second + '0'); str2 += ']';
			str2 += ' ';
		}
		str2 += ")";

		printf("GraphArticulationHandler::FindArticulationPoints, %s, the answer is %s, %s \n", (b ? "Correct" : "Wrong"), str.c_str(), str2.c_str());
	}

	static void Test_03()
	{
		Graph graph;

		//
		// Ref: http://www.geeksforgeeks.org/articulation-points-or-cut-vertices-in-a-graph/
		//
		graph.V = 7;
		graph.E = 8;

		graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

		graph.AddEdgeToAdjacencyList(false, 0, 1);
		graph.AddEdgeToAdjacencyList(false, 0, 2);
		graph.AddEdgeToAdjacencyList(false, 1, 2);
		graph.AddEdgeToAdjacencyList(false, 1, 3);
		graph.AddEdgeToAdjacencyList(false, 1, 4);
		graph.AddEdgeToAdjacencyList(false, 1, 6);
		graph.AddEdgeToAdjacencyList(false, 3, 5);
		graph.AddEdgeToAdjacencyList(false, 4, 5);

		set<int> APs;
		set<pair<int, int>> bridges;

		GraphArticulationHandler gah;
		gah.FindArticulationPoints(graph, 0, &APs, &bridges);

		set<int> ap_answer = { 1 };
		set<pair<int, int>> bridge_answer = { { 0, 1 },{ 1, 2 },{ 1, 3 },{ 1, 4 },{ 1, 6 } };

		bool b = APs == ap_answer && bridges == bridge_answer;

		string str = "Articulation Points ( ";
		for (auto i = APs.begin(); i != APs.end(); i++) {
			str += *i + '0';
			str += ' ';
		}
		str += ")";

		string str2 = "Bridges ( ";
		for (auto i = bridges.begin(); i != bridges.end(); i++) {
			str2 += '[';  str2 += (i->first + '0'); str2 += " -> "; str2 += (i->second + '0'); str2 += ']';
			str2 += ' ';
		}
		str2 += ")";

		printf("GraphArticulationHandler::FindArticulationPoints, %s, the answer is %s, %s \n", (b ? "Correct" : "Wrong"), str.c_str(), str2.c_str());
	}
};



///////////////////////////////////////////////////////////////////////////////////////////


class GraphWeightHandler
{
public:

    int V;
    int E;

    int source_v;

    //
    // traversal[i].state is the state of vertex "i"
    // traversal[i].accum_weight is the accumulated weight from source_v to vertex "i", for shortest path, etc.
    // traversal[i].source is the parent of vertex "i" (in directed graph), value -1 means no parent
    //
    vector<GraphTraversal> traversal;

    //
    // For shortest path search such as BellmanFord or Dijkstra, paths[i] means the found path from source_v to vertex "i".
    // For minimum spanning tree such as Prime or Kruskal, Path[0] is used for the found path
    //
    vector<GraphPath> paths;	

    bool HasNegativeCycle;

public:

    GraphWeightHandler()
    {
        Clear();
    }

    void Clear()
    {
        this->V = 0;
        this->E = 0;

        this->source_v = 0;

        this->traversal.clear();

        this->paths.clear();

        this->HasNegativeCycle = false;
    }

public:

    //
    // Find the minimum spanning tree (MST)
    //
    // By nature it's a greedy algorithm
    // Works for undirected graph only
    // Time complexity is O(V^2)
    //
    double PrimAlgorithm(Graph graph, int source)
    {
        Clear();

        this->V = graph.V;
        this->E = graph.E;

        this->source_v = source;
        this->traversal.assign(this->V, GraphTraversal(GRAPH_UNVISITED));

        this->paths.assign(1, GraphPath());

        struct GreaterThanCompare {
            bool operator()(const GraphEdge& l, const GraphEdge& r) {
                return l.w > r.w;
            }
        };

        //
        // GreaterThanCompare is used so that priority_queue.top() returns the smallest instead
        // 
        priority_queue<GraphEdge, vector<GraphEdge>, GreaterThanCompare> minHeap;

        this->traversal[this->source_v].state = GRAPH_VISITED;

        //
        // add all the edges connecting "start" vertext to the priority queue
        //
        for (int j = 0; j < (int)graph.AdjacencyList[this->source_v].size(); j++)
        {
            GraphVertex& vertex = graph.AdjacencyList[this->source_v][j];

            if (this->traversal[vertex.v].state == GRAPH_UNVISITED)
            {
                minHeap.push(GraphEdge(this->source_v, vertex.v, vertex.w));
            }
        }

        while (!minHeap.empty() && ((int)this->paths[0].edges.size() + 1) < this->V) // check if all vertice included
        {
            //
            // pop up the edge with the least weight
            //
            GraphEdge edge = minHeap.top();
            minHeap.pop();

            //
            // if the two vertices from the edge are already picked, we must discard them otherwise a loop is formed 
            //
            if (this->traversal[edge.v1].state != GRAPH_UNVISITED && this->traversal[edge.v2].state != GRAPH_UNVISITED)	{
                continue;
            }

            // 
            // this is an edge in the MST
            //
            this->paths[0].edges.push_back(edge);
            this->paths[0].weight_sum += edge.w;

            //
            // now we examine all edges connecting to the two vertices of the selected edge
            //
            for (int i = 0; i < 2; i++)
            {
                int u = (i==0) ? edge.v1 : edge.v2;

                if (this->traversal[u].state != GRAPH_UNVISITED)
                {
                    continue;
                }

                this->traversal[u].state = GRAPH_VISITED;

                //
                // add all the edges connecting "u" vertex to the priority queue
                //
                for (int j = 0; j < (int)graph.AdjacencyList[u].size(); j++)
                {
                    GraphVertex& vertex = graph.AdjacencyList[u][j];

                    if (this->traversal[vertex.v].state == GRAPH_UNVISITED)
                    {
                        minHeap.push(GraphEdge(u, vertex.v, vertex.w));
                    }
                }
            }
        }

        return this->paths[0].weight_sum;
    }

    //
    // Find the minimum spanning tree (MST)
    //
    double KruskalAlgorithm(Graph graph)
    {
        Clear();

        this->V = graph.V;
        this->E = graph.E;

        this->paths.assign(1, GraphPath());

        struct GreaterThanCompare {
            bool operator()(const GraphEdge& l, const GraphEdge& r) {
                return l.w > r.w;
            }
        };

        //
        // GreaterThanCompare is used so that priority_queue.top() returns the smallest instead
        // 
        priority_queue<GraphEdge, vector<GraphEdge>, GreaterThanCompare> minHeap;

        //
        // add all the edges to the priority queue
        //
        for (int i = 0; i < this->V; i++)
        {
            for (int j = 0; j < (int)graph.AdjacencyList[i].size(); j++)
            {
                GraphVertex& vertex = graph.AdjacencyList[i][j];

                minHeap.push(GraphEdge(i, vertex.v, vertex.w));
            }
        }

        GraphUnionFind UF(this->V);

        //
        // if we have processed V-1 edges, we are done
        //
        while (!minHeap.empty() && ((int)this->paths[0].edges.size() + 1) < this->V)
        {
            //
            // pop up the edge with the least weight
            //
            GraphEdge edge = minHeap.top();
            minHeap.pop();


            //
            // if the two vertices from the edge are already picked, we must discard them otherwise a loop is formed 
            //
            if (UF.AreInSameSet(edge.v1, edge.v2))
            {
                continue;
            }

            UF.Union(edge.v1, edge.v2);

            // 
            // this is an edge in the MST
            //
            this->paths[0].edges.push_back(edge);
            this->paths[0].weight_sum += edge.w;
        }

        return this->paths[0].weight_sum;
    }

    //
    // Find the single-sourced shortest path (SSSP)
    //
    // Time complexity is O(EV)
    // Can't handle negative cycle
    //
    void BellmanFordAlgorithm(Graph& graph, int source)
    {
        Clear();

        this->V = graph.V;
        this->E = graph.E;

        this->source_v = source;

        this->traversal.assign(this->V, GraphTraversal(GRAPH_UNVISITED, -1, LARGEST_DOUBLE));
        this->traversal[this->source_v].accum_weight = 0;

        for (int i = 0; i < V - 1; i++)  // relax all E edges V-1 times, overall O(VE)
        {
            for (int vertexFrom = 0; vertexFrom < V; vertexFrom++)  // these two loops = O(E)
            {
                for (int k = 0; k < (int)graph.AdjacencyList[vertexFrom].size(); k++)
                {
                    GraphVertex& to = graph.AdjacencyList[vertexFrom][k];

                    if (this->traversal[to.v].accum_weight >(this->traversal[vertexFrom].accum_weight + to.w))
                    {
                        this->traversal[to.v].accum_weight = this->traversal[vertexFrom].accum_weight + to.w;  // relax

                        this->traversal[to.v].source_v = vertexFrom;
                    }
                }
            }
        }

        this->HasNegativeCycle = false;

        for (int vertexFrom = 0; (!this->HasNegativeCycle) && vertexFrom < this->V; vertexFrom++)   // one more pass to check
        {
            for (int k = 0; (!this->HasNegativeCycle) && k < (int)graph.AdjacencyList[vertexFrom].size(); k++)
            {
                GraphVertex& to = graph.AdjacencyList[vertexFrom][k];

                if (this->traversal[to.v].accum_weight > this->traversal[vertexFrom].accum_weight + to.w)  // should be false
                {
                    this->HasNegativeCycle = true;     // but if true, then negative cycle exists!
                }
            }
        }
    }

    //
    // Find the single-sourced shortest path (SSSP)
    //
    // Time complexity is O( E*log(E) + V )
    // Can't handle negative edge
    //
    void DijkstraAlgorithm(Graph& graph, int source)
    {
        Clear();

        this->V = graph.V;
        this->E = graph.E;

        this->source_v = source;

        this->traversal.assign(this->V, GraphTraversal(GRAPH_UNVISITED, -1, LARGEST_DOUBLE));
        this->traversal[this->source_v].accum_weight = 0;

        struct GreaterThanCompare {
            bool operator()(const GraphEdge& l, const GraphEdge& r) {
                return l.w > r.w;
            }
        };

        //
        // GreaterThanCompare is used so that priority_queue.top() returns the smallest instead
        // 
        priority_queue<GraphEdge, vector<GraphEdge>, GreaterThanCompare> minHeap;

        minHeap.push(GraphEdge(this->source_v, this->source_v, 0));

        while (!minHeap.empty())
        {
            //
            // pop up the edge with the least weight
            //
            GraphEdge edge = minHeap.top();
            minHeap.pop();

            if (edge.w > this->traversal[edge.v2].accum_weight)
            {
                continue;
            }

            for (int k = 0; k < (int)graph.AdjacencyList[edge.v2].size(); k++)
            {
                GraphVertex& to = graph.AdjacencyList[edge.v2][k];

                if (this->traversal[to.v].accum_weight >(this->traversal[edge.v2].accum_weight + to.w))
                {
                    this->traversal[to.v].accum_weight = this->traversal[edge.v2].accum_weight + to.w;  // relax

                    minHeap.push(GraphEdge(edge.v2, to.v, to.w));

                    this->traversal[to.v].source_v = edge.v2;
                }
            }
        }
    }

    double DijkstraAlgorithm(Graph& graph, int source, int dest)
    {
        DijkstraAlgorithm(graph, source);

        ResolveSearchPaths();

        if (this->paths[dest].IsEmpty())
        {
            return -1;
        }

        return this->paths[dest].weight_sum;
    }

    //
    // Shortest Path Faster Algorithm (SPFA)
    //
    // Can handle negative edge. Suitable for sparse graph. Otherwise in worst cases perform worse than DijkstraAlgorithm.
    //
    void SPFA(Graph& graph, int source)
    {
        Clear();

        this->V = graph.V;
        this->E = graph.E;

        this->source_v = source;

        this->traversal.assign(this->V, GraphTraversal(GRAPH_UNVISITED, -1, LARGEST_DOUBLE));
        this->traversal[this->source_v].accum_weight = 0;

        queue<GraphEdge> que;

        que.push(GraphEdge(this->source_v, this->source_v, 0));

        while (!que.empty())
        {
            GraphEdge edge = que.front();
            que.pop();

            this->traversal[edge.v2].state = GRAPH_UNVISITED;		// the same node can be pushed into the queue again

            for (int k = 0; k < (int)graph.AdjacencyList[edge.v2].size(); k++)
            {
                GraphVertex& to = graph.AdjacencyList[edge.v2][k];

                if (this->traversal[to.v].accum_weight >(this->traversal[edge.v2].accum_weight + to.w))
                {
                    this->traversal[to.v].accum_weight = this->traversal[edge.v2].accum_weight + to.w;  // relax

                    this->traversal[to.v].source_v = edge.v2;

                    if (this->traversal[to.v].state == GRAPH_UNVISITED)
                    {
                        this->traversal[to.v].state = GRAPH_VISITED;

                        que.push(GraphEdge(edge.v2, to.v, to.w));
                    }
                }
            }
        }
    }

    double SPFA(Graph& graph, int source, int dest)
    {
        SPFA(graph, source);

        ResolveSearchPaths();

        if (this->paths[dest].IsEmpty())
        {
            return -1;
        }

        return this->paths[dest].weight_sum;
    }

    //
    // table[i][j] has the shortest path from vertex i to j
    //
    void AllPairsDijkstraAlgorithm(Graph& graph, vector< vector<GraphPath> >& table)
    {
        table.clear();

        this->V = graph.V;
        this->E = graph.E;

        for (int i = 0; i < this->V; i++)
        {
            DijkstraAlgorithm(graph, i);

            ResolveSearchPaths();

            table.push_back(this->paths);
        }
    }

    //
    // Find the all-pairs shortest pathes 
    //
    // Time complexity is O( V*V*V )
    // Suitable for dense graph
    // Can't handle negative cycle
    //
    static void FloydWarshallAlgorithm(Graph& graph, vector< vector<double> >& table)
    {
        table.clear();

        vector<double> vec;
        vec.assign(graph.V, LARGEST_DOUBLE);

        table.assign(graph.V, vec);

        for (int i = 0; i < graph.V; i++) 
        {
            table[i][i] = 0;
        }

        for (int i = 0; i < graph.V; i++)
        {
            vector< GraphVertex >& list = graph.AdjacencyList[i];

            for (int j = 0; j < (int)list.size(); j++)
            {
                GraphVertex& item = list[j];

                table[i][item.v] = item.w;
            }
        }

        for (int k = 0; k < graph.V; k++) // common error: remember that loop order is k->i->j
        {
            for (int i = 0; i < graph.V; i++)
            {
                for (int j = 0; j < graph.V; j++)
                {
                    table[i][j] = min(table[i][j], table[i][k] + table[k][j]);
                }
            }
        }
    }

    //
    // Find the all-pairs shortest pathes making use of both Bellman-Ford and Dijkstra.
    //
    // Time complexity is O( V*V*logV + VE ), for a complete graph (E=V*V), same as Floyd-Warshall
    // Suitable for sparse graph
    // designed to handle negative edges
    //
    void JohnsonAlgorithm(Graph graph, vector< vector<GraphPath> >& table)
    {
        table.clear();

        //
        // Add a new vertex to the graph, add edges (with zero weight) from the new vertex to all vertices
        //
        graph.AdjacencyList.push_back(vector<GraphVertex>());
        
        vector<GraphVertex>& edges = graph.AdjacencyList[graph.V];

        for (int i = 0; i < graph.V; i++)
        {
            edges.push_back( GraphVertex(i, 0) );
        }
        
        graph.E += graph.V;
        graph.V++;

        //
        // Run Bellman-Ford algorithm the new vertex as source
        //
        BellmanFordAlgorithm(graph, graph.V-1);

        graph.AdjacencyList.pop_back();

        graph.V--;
        graph.E -= graph.V;
        
        if (this->HasNegativeCycle)
        {
            return;
        }

        vector<double> weights(this->traversal.size());

        for (int i = 0; i < (int)weights.size(); i++) {
            weights[i] = this->traversal[i].accum_weight;
        }

        //
        // re-weight the graph by:
        //
        // w(u, v) = w(u, v) + h[u] – h[v], h[] are shortest distances computed by Bellman-Ford 
        //
        for (int from = 0; from < graph.V; from++)
        {
            vector<GraphVertex>& edges = graph.AdjacencyList[from];

            for (int j = 0; j < (int)edges.size(); j++)
            {
                GraphVertex& item = edges[j];
                int to = item.v;

                item.w = item.w + weights[from] - weights[to];
            }
        }

        //
        // Now all edges are non-negative, simply run Dijkstra
        //
        AllPairsDijkstraAlgorithm(graph, table);

        //
        // Restore the real weight sum for the paths
        //
        for (int i = 0; i < graph.V; i++)
        {
            for (int j = 0; j < graph.V; j++)
            {
                GraphPath& path = table[i][j];

                path.weight_sum = 0;

                for (int m = 0; m < path.GetEdgeCount(); m++)
                {
                    path.weight_sum += weights[path.edges[m].v2] - weights[path.edges[m].v1];
                }

            }
        }
    }

    void ShortestPathOnDAGAlgorithm(Graph& graph, int source)
    {
		vector<int> vertexState(graph.V, GRAPH_UNVISITED);
		vector<int> topology;

        GraphTraversalHandler::TopologicalSortByDfsSimpleRecursive(graph, vertexState, topology);

        this->V = graph.V;
        this->E = graph.E;

        this->source_v = source;

        this->traversal.assign(this->V, GraphTraversal(GRAPH_UNVISITED, -1, LARGEST_DOUBLE));
        this->traversal[this->source_v].accum_weight = 0;

        for (int i = 0; i < (int)topology.size(); i++)
        {
            int v1 = topology[i];

            for (int k = 0; k < (int)graph.AdjacencyList[v1].size(); k++)
            {
                GraphVertex& item = graph.AdjacencyList[v1][k];

                int v2 = item.v;

                if (item.w > this->traversal[v2].accum_weight)
                {
                    continue;
                }

                if (this->traversal[v2].accum_weight > (this->traversal[v1].accum_weight + item.w))
                {
                    this->traversal[v2].accum_weight = traversal[v1].accum_weight + item.w;  // relax

                    this->traversal[v2].source_v = v1;
                }
            }
        }
    }

    void LongestPathOnDAGAlgorithm(Graph& graph, int source)
    {
		vector<int> vertexState(graph.V, GRAPH_UNVISITED);
		vector<int> topology;

		GraphTraversalHandler::TopologicalSortByDfsSimpleRecursive(graph, vertexState, topology);

        this->V = graph.V;
        this->E = graph.E;

        this->source_v = source;

        this->traversal.assign(this->V, GraphTraversal(GRAPH_UNVISITED, -1, SMALLEST_DOUBLE));
        this->traversal[this->source_v].accum_weight = 0;

        for (int i = 0; i < (int)topology.size(); i++)
        {
            int v1 = topology[i];

            for (int k = 0; k < (int)graph.AdjacencyList[v1].size(); k++)
            {
                GraphVertex& item = graph.AdjacencyList[v1][k];

                int v2 = item.v;

                if (item.w < this->traversal[v2].accum_weight) {
                    continue;
                }

                if (this->traversal[v2].accum_weight < (this->traversal[v1].accum_weight + item.w))
                {
                    this->traversal[v2].accum_weight = this->traversal[v1].accum_weight + item.w;  // relax

                    this->traversal[v2].source_v = v1;
                }
            }
        }
    }

    //
    // After running BellmanFordAlgorithm or DijkstraAlgorithm, 
    // get the edges of each shortest path. 
    // paths[i] are the edges of the shortest path from source_v to vertex i.
    //
    void ResolveSearchPaths()
    {
        this->paths.assign(this->V, GraphPath());

        for (int i = 0; i < this->V; i++)
        {
            this->paths[i].source_v = this->source_v;

            if (!(this->traversal[i].accum_weight < LARGEST_DOUBLE && this->traversal[i].accum_weight > SMALLEST_DOUBLE))
            {
                continue; // unreachable
            }

            this->paths[i].weight_sum = this->traversal[i].accum_weight;

            if (this->paths[i].source_v == i)
            {
                this->paths[i].edges.push_back(GraphEdge(i, i, 0));
                continue;
            }

            int j = i;
            while (j >= 0)
            {
                int parent = this->traversal[j].source_v;

                if (parent >= 0)
                {
                    //
                    // TODO: find the k that AdjacencyList[parent][k] holds the vertex j
                    //
                    //int k = 0;

                    //GraphVertex& v = graph.AdjacencyList[parent][k];

                    //this->paths[i].edges.push_back(GraphEdge(parent, j, v.Weight));

                    this->paths[i].edges.push_back(GraphEdge(parent, j, 0));
                }

                j = parent;
            }

            reverse(this->paths[i].edges.begin(), this->paths[i].edges.end());
        }
    }

public:

	static void Test()
	{
		TestPrimeAlgorithm();

		TestKruskalAlgorithm();

		TestBellmanFordAlgorithm();

		TestDijkstraAlgorithm();

		TestSPFA();

		TestFloydWarshallAlgorithm();

		TestJohnsonAlgorithm();

		TestShortestPathOnDAGAlgorithm();
	}

    static void TestPrimeAlgorithm()
    {
        Graph graph;
		GraphHelper graphHelper;

		graphHelper.InputWeightedGraph(graph, false, "input_graph_mst_tiny.txt");

        GraphWeightHandler gwh;

        gwh.PrimAlgorithm(graph, 0);

		graphHelper.OpenOutput("output_graph_mst_tiny.txt");

        printf("MST = %lf \n", gwh.paths[0].weight_sum);

        for (int i = 0; i < (int)gwh.paths[0].edges.size(); i++)
        {
            printf("(%d %d %lf) \n", gwh.paths[0].edges[i].v1, gwh.paths[0].edges[i].v2, gwh.paths[0].edges[i].w);
        }

		graphHelper.CloseOutput();
    }

    static void TestKruskalAlgorithm()
    {
        Graph graph;
		GraphHelper graphHelper;

		graphHelper.InputWeightedGraph(graph, true, "input_graph_mst_tiny.txt");

        GraphWeightHandler gwh;

        gwh.KruskalAlgorithm(graph);

		graphHelper.OpenOutput("output_graph_mst_tiny.txt");

        printf("MST = %lf \n", gwh.paths[0].weight_sum);

        for (int i = 0; i < (int)gwh.paths[0].edges.size(); i++)
        {
            printf("(%d %d %lf) \n", gwh.paths[0].edges[i].v1, gwh.paths[0].edges[i].v2, gwh.paths[0].edges[i].w);
        }

		graphHelper.CloseOutput();
    }

    static void TestBellmanFordAlgorithm()
    {
        Graph graph;
		GraphHelper graphHelper;

		graphHelper.InputWeightedGraph(graph, false, "input_graph_bellman_ford.txt");

        GraphWeightHandler gwh;

        gwh.BellmanFordAlgorithm(graph, 0);
        gwh.ResolveSearchPaths();

		graphHelper.OpenOutput("output_graph_bellman_ford.txt");

        printf("Negative Cycle Exist? %s\n", gwh.HasNegativeCycle ? "Yes" : "No");

        if (!gwh.HasNegativeCycle)
        {
            for (int i = 0; i < (int)gwh.paths.size(); i++)
            {
                if (gwh.paths[i].IsEmpty())
                {
                    printf("BellmanFord[%d -> %d]: unreachable \n", gwh.paths[i].source_v, i);
                    continue;
                }

                printf("BellmanFord[%d -> %d]: %lf, (%d", gwh.paths[i].source_v, i, gwh.paths[i].weight_sum, gwh.paths[i].source_v);

                for (int j = 0; j < gwh.paths[i].GetEdgeCount(); j++)
                {
                    printf(" - %d", gwh.paths[i].edges[j].v2);
                }

                printf(")\n");
            }
        }

		graphHelper.CloseOutput();
    }

    static void TestDijkstraAlgorithm()
    {
        //
        // Reference: https://www.youtube.com/watch?v=5GT5hYzjNoo
        //

        Graph graph;
		GraphHelper graphHelper;

		graphHelper.InputWeightedGraph(graph, false, "input_graph_dijkstra.txt");

        GraphWeightHandler gwh;

        double weight_0_7 = gwh.DijkstraAlgorithm(graph, 0, 7);

        //gwh.ResolveSearchPaths();

		graphHelper.OpenOutput("output_graph_dijkstra.txt");

        for (int i = 0; i < (int)gwh.paths.size(); i++)
        {
            if (gwh.paths[i].IsEmpty())
            {
                printf("Dijkstra[%d -> %d]: unreachable \n", gwh.paths[i].source_v, i);
                continue;
            }

            printf("Dijkstra[%d -> %d]: %lf, (%d", gwh.paths[i].source_v, i, gwh.paths[i].weight_sum, gwh.paths[i].source_v);

            for (int j = 0; j < gwh.paths[i].GetEdgeCount(); j++)
            {
                printf(" - %d", gwh.paths[i].edges[j].v2);
            }

            printf(")\n");
        }

        printf("The shortest path from 0(A) to 7(H) should be 11, the answer is %lf \n", weight_0_7);

		graphHelper.CloseOutput();
    }

    static void TestSPFA()
    {
        //
        // Reference: https://www.youtube.com/watch?v=5GT5hYzjNoo
        //

        Graph graph;
		GraphHelper graphHelper;

		graphHelper.InputWeightedGraph(graph, false, "input_graph_SPFA.txt");

        GraphWeightHandler gwh;

        double weight_0_7 = gwh.SPFA(graph, 0, 7);

        //gwh.ResolveSearchPaths();

		graphHelper.OpenOutput("output_graph_SPFA.txt");

        for (int i = 0; i < (int)gwh.paths.size(); i++)
        {
            if (gwh.paths[i].IsEmpty())
            {
                printf("SPFA[%d -> %d]: unreachable \n", gwh.paths[i].source_v, i);
                continue;
            }

            printf("SPFA[%d -> %d]: %lf, (%d", gwh.paths[i].source_v, i, gwh.paths[i].weight_sum, gwh.paths[i].source_v);

            for (int j = 0; j < gwh.paths[i].GetEdgeCount(); j++)
            {
                printf(" - %d", gwh.paths[i].edges[j].v2);
            }

            printf(")\n");
        }

        printf("The shortest path from 0(A) to 7(H) should be 11, the answer is %lf\n", weight_0_7);

		graphHelper.CloseOutput();
    }

    static void TestFloydWarshallAlgorithm()
    {
        Graph graph;
		GraphHelper graphHelper;

		graphHelper.InputWeightedGraph(graph, true, "input_graph_floyd_warshall.txt");

        vector< vector<double> > table;
        GraphWeightHandler::FloydWarshallAlgorithm(graph, table);

		graphHelper.OpenOutput("output_graph_floyd_warshall.txt");

        for (int i = 0; i < graph.V; i++)
        {
            for (int j = 0; j < graph.V; j++)
            {
                printf("FloydWarshall[%d, %d] = ", i, j);

                if (FLOAT_EQUAL(table[i][j], LARGEST_DOUBLE))
                {
                    printf("* \n"); // unreachable
                }
                else
                {
                    printf("%lf\n", table[i][j]);
                }
            }
        }

		graphHelper.CloseOutput();
    }

    static void TestJohnsonAlgorithm()
    {
        Graph graph;
		GraphHelper graphHelper;

		graphHelper.InputWeightedGraph(graph, true, "input_graph_johnson.txt");

        vector< vector<GraphPath> > table;
        GraphWeightHandler gwh;
        
        gwh.JohnsonAlgorithm(graph, table);

		graphHelper.OpenOutput("output_graph_johnson.txt");

        for (int k = 0; k < graph.V; k++)
        {
            for (int j = 0; j < graph.V; j++)
            {
                GraphPath& path = table[k][j];

                if (path.IsEmpty())
                {
                    printf("Johnson[%d -> %d]: unreachable \n", k, j);
                    continue;
                }

                printf("Johnson[%d -> %d]: %lf, (%d", k, j, path.weight_sum, path.source_v);

                for (int m = 0; m < path.GetEdgeCount(); m++)
                {
                    printf(" - %d", path.edges[m].v2);
                }

                printf(")\n");
            }
        }

		graphHelper.CloseOutput();
    }

    static void TestShortestPathOnDAGAlgorithm()
    {
        Graph graph;
		GraphHelper graphHelper;

		graphHelper.InputWeightedGraph(graph, true, "input_graph_EWDAG_tiny.txt");

        GraphWeightHandler gwh;

        gwh.ShortestPathOnDAGAlgorithm(graph, 5);
        gwh.ResolveSearchPaths();

		graphHelper.OpenOutput("output_graph_EWDAG_tiny.txt");

        for (int i = 0; i < (int)gwh.paths.size(); i++)
        {
            if (gwh.paths[i].IsEmpty())
            {
                printf("SPDAG[%d -> %d]: unreachable \n", gwh.paths[i].source_v, i);
                continue;
            }

            printf("SPDAG[%d -> %d]: %lf, (%d", gwh.paths[i].source_v, i, gwh.paths[i].weight_sum, gwh.paths[i].source_v);

            for (int j = 0; j < gwh.paths[i].GetEdgeCount(); j++)
            {
                printf(" - %d", gwh.paths[i].edges[j].v2);
            }

            printf(")\n");
        }

		graphHelper.CloseOutput();
    }

};


///////////////////////////////////////////////////////////////////////////////////////////


class GraphEulerPathHandler
{
public:

    int V;
    int E;

    bool GraphConnected;
    int NumOfOddVertices;
    int OddVertices[2];
    bool HasPath;
    bool HasCircuit;
    vector<int> Path;

    vector<int> NumOfRemainingNeighbors;	// performance optimization only, NumOfRemainingNeighbors[i] stores the number of vertex i's un-removed vertexes

public:

    GraphEulerPathHandler()
    {
        this->V = 0;
        this->E = 0;
    }

    bool IsConnected(Graph& graph)
    {
        queue<int> q;

        vector<bool> visited(graph.V, false);

        q.push(0);

        while (!q.empty())
        {
            int v = q.front();
            q.pop();

            visited[v] = true;

            for (int i = 0; i < (int)graph.AdjacencyList[v].size(); i++)
            {
                int u = graph.AdjacencyList[v][i].v;

                if (!visited[u]) {
                    q.push(u);
                }
            }
        }

        for (int v = 0; v < graph.V; v++)
        {
            if (!visited[v]) {
                return false;
            }
        }

        return true;
    }

    void EulerianPathCheck(Graph& graph)
    {
        this->V = graph.V;
        this->E = graph.E;

        this->NumOfOddVertices = 0;
        this->OddVertices[0] = this->OddVertices[1] = -1;
        this->HasPath = false;
        this->HasCircuit = false;
        this->Path.clear();
        
        this->GraphConnected = IsConnected(graph);
        if (!this->GraphConnected) {
            return;
        }

        this->NumOfRemainingNeighbors.assign(graph.V, 0);

        for (int from = 0; from < (int)graph.AdjacencyList.size(); from++)
        {
            vector<GraphVertex>& neighbors = graph.AdjacencyList[from];

            int num_of_neighbors = (int)neighbors.size();

            if (num_of_neighbors & 1) 
            {
                this->NumOfOddVertices++;

                if (this->OddVertices[0] == -1) {
                    this->OddVertices[0] = from;
                }
                else {
                    this->OddVertices[1] = from;
                }
            }

            if (num_of_neighbors > 0) {
                this->NumOfRemainingNeighbors[from] = num_of_neighbors;
            }
        }

        this->HasCircuit = (this->NumOfOddVertices == 0);

        if (this->HasCircuit) {
            this->HasPath = true;
        }
        else {
            this->HasPath = (this->NumOfOddVertices == 2);
        }
    }

    bool IsBridgeEdge(Graph& graph, int from, int to)
    {
        return this->NumOfRemainingNeighbors[to] == 1;
    }

    void RemoveEdge(Graph& graph, int from, int to)
    {
        vector<GraphVertex>& neighbors = graph.AdjacencyList[from];

        for(int i=0; i< (int)neighbors.size(); i++)
        {
            GraphVertex& neighbor = neighbors[i];

            if (neighbor.v == to && neighbor.w >= 0)
            {
                neighbor.w = -1;
                this->NumOfRemainingNeighbors[from]--;
                break;
            }
        }

        vector<GraphVertex>& neighbors2 = graph.AdjacencyList[to];

        for (int i = 0; i< (int)neighbors2.size(); i++)
        {
            GraphVertex& neighbor2 = neighbors2[i];

            if (neighbor2.v == from && neighbor2.w >= 0)
            {
                neighbor2.w = -1;
                this->NumOfRemainingNeighbors[to]--;
                break;
            }
        }
    }

    void AdjustPathByLexicalOrder()
    {
        int size = (int)this->Path.size();
        for (int i = 0; i < size / 2; i++)
        {
            if (this->Path[i] > this->Path[size - i - 1])
            {
                reverse(this->Path.begin(), this->Path.end());
                break;
            }
        }
    }

    //
    // Find Eulerian path or circuit using the concept of "Don't Burn the Bridge"
    //
    void FleuryAlgorithm(Graph& graph)
    {
        EulerianPathCheck(graph);

        if (!this->HasPath) {
            return;
        }

        int hop = this->HasCircuit ? 0 : this->OddVertices[0];

        this->Path.push_back(hop);

        while (this->NumOfRemainingNeighbors[hop] > 0)
        {
            vector<GraphVertex>& neighbors = graph.AdjacencyList[hop];

            bool only_one_left = this->NumOfRemainingNeighbors[hop] == 1;

            int i;
            for (i = 0; i < (int)neighbors.size(); i++)
            {
                if (only_one_left) 
                {
                    if (neighbors[i].w >= 0) { // note we reuse the weight field for the removal flag
                        break;
                    }
                }
                else 
                {
                    if (neighbors[i].w >= 0 && ! IsBridgeEdge(graph, hop, neighbors[i].v)) {
                        break;
                    }
                }
            }

            int next = neighbors[i].v;

            RemoveEdge(graph, hop, next);

            this->Path.push_back(next);

            hop = next;
        }
    }

    //
    // Find Eulerian path or circuit using a stack
    //
    void EulerianPathAlgorithm(Graph& graph)
    {
        EulerianPathCheck(graph);

        if (!this->HasPath) {
            return;
        }

        stack<int> stack;

        int hop = this->HasCircuit ? 0 : min(this->OddVertices[0], this->OddVertices[1]);

        while (hop >= 0)
        {
            if (this->NumOfRemainingNeighbors[hop] == 0)
            {
                this->Path.push_back(hop);

                if (!stack.empty())
                {
                    hop = stack.top();
                    stack.pop();
                }
                else {
                    hop = -1;
                }
            }
            else
            {
                vector<GraphVertex>& neighbors = graph.AdjacencyList[hop];

                int i;
                for (i = 0; i < (int)neighbors.size(); i++)
                {
                    if (neighbors[i].w >= 0) {	// note we reuse the weight field for the removal flag
                        break;
                    }
                }

                int next = neighbors[i].v;

                RemoveEdge(graph, hop, next);

                stack.push(hop);

                hop = next;
            }
        }
    }

public:

    static void Test()
    {
        Test_EulerianPath_1();
        Test_EulerianPath_2();
        Test_EulerianPath_3();
        Test_EulerianPath_4();
        Test_EulerianPath_5();
        Test_EulerianPath_6();
    }

    //
    // reference: http://www.geeksforgeeks.org/eulerian-path-and-circuit/
    //
    static void Test_EulerianPath_1()
    {
        Graph graph;

        graph.V = 5;
        graph.E = 5;

        graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

        graph.AddEdgeToAdjacencyList(false, 0, 1);
        graph.AddEdgeToAdjacencyList(false, 0, 2);
        graph.AddEdgeToAdjacencyList(false, 0, 3);
        graph.AddEdgeToAdjacencyList(false, 1, 2);
        graph.AddEdgeToAdjacencyList(false, 3, 4);

        GraphEulerPathHandler handler;

        handler.FleuryAlgorithm(graph);
    }

    //
    // reference: http://www.geeksforgeeks.org/eulerian-path-and-circuit/
    //
    static void Test_EulerianPath_2()
    {
        Graph graph;

        graph.V = 5;
        graph.E = 6;

        graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

        graph.AddEdgeToAdjacencyList(false, 0, 1);
        graph.AddEdgeToAdjacencyList(false, 0, 2);
        graph.AddEdgeToAdjacencyList(false, 0, 3);
        graph.AddEdgeToAdjacencyList(false, 0, 4);
        graph.AddEdgeToAdjacencyList(false, 1, 2);
        graph.AddEdgeToAdjacencyList(false, 3, 4);

        GraphEulerPathHandler handler;

        handler.FleuryAlgorithm(graph);
    }

    //
    // reference: http://www.geeksforgeeks.org/eulerian-path-and-circuit/
    //
    static void Test_EulerianPath_3()
    {
        Graph graph;

        graph.V = 5;
        graph.E = 6;

        graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

        graph.AddEdgeToAdjacencyList(false, 0, 1);
        graph.AddEdgeToAdjacencyList(false, 0, 2);
        graph.AddEdgeToAdjacencyList(false, 0, 3);
        graph.AddEdgeToAdjacencyList(false, 1, 2);
        graph.AddEdgeToAdjacencyList(false, 1, 3);
        graph.AddEdgeToAdjacencyList(false, 3, 4);

        GraphEulerPathHandler handler;

        handler.FleuryAlgorithm(graph);
    }

    //
    // reference: Ref 02 EulerianPathAndCircuit.pdf
    //
    static void Test_EulerianPath_4()
    {
        Graph graph;

        graph.V = 6;
        graph.E = 10;

        graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

        int F = 0;
        int E = 1;
        int C = 2;
        int D = 3;
        int A = 4;
        int B = 5;

        char cc[] = {'F', 'E', 'C', 'D', 'A', 'B' };

        graph.AddEdgeToAdjacencyList(false, F, E);
        graph.AddEdgeToAdjacencyList(false, F, C);
        graph.AddEdgeToAdjacencyList(false, F, D);
        graph.AddEdgeToAdjacencyList(false, A, E);
        graph.AddEdgeToAdjacencyList(false, A, C);
        graph.AddEdgeToAdjacencyList(false, A, B);
        graph.AddEdgeToAdjacencyList(false, C, B);
        graph.AddEdgeToAdjacencyList(false, C, D);
        graph.AddEdgeToAdjacencyList(false, B, D);
        graph.AddEdgeToAdjacencyList(false, B, D);

        // pick vertices with smaller indices first
        for (int i = 0; i < (int)graph.AdjacencyList.size(); i++) {
            sort(graph.AdjacencyList[i].begin(), graph.AdjacencyList[i].end(), GraphVertex::IsIndexSmaller);
        }

        GraphEulerPathHandler handler;

        handler.FleuryAlgorithm(graph);

        string str;
        for (int i = 0; i < (int)handler.Path.size(); i++) {
            str += cc[handler.Path[i]];
        }

        printf("FleuryAlgorithm: %s, the answer is %s\n", (str == "FEACFDCBDBA") ? "Correct" : "Wrong", str.c_str());

    }

    //
    // reference: Ref 01 EulerianTour.pdf
    //
    static void Test_EulerianPath_5()
    {
        Graph graph;

        graph.V = 7;
        graph.E = 12;

        graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

        graph.AddEdgeToAdjacencyList(false, 0, 4);
        graph.AddEdgeToAdjacencyList(false, 0, 3);
        graph.AddEdgeToAdjacencyList(false, 1, 4);
        graph.AddEdgeToAdjacencyList(false, 1, 6);
        graph.AddEdgeToAdjacencyList(false, 1, 5);
        graph.AddEdgeToAdjacencyList(false, 1, 3);
        graph.AddEdgeToAdjacencyList(false, 4, 6);
        graph.AddEdgeToAdjacencyList(false, 4, 5);
        graph.AddEdgeToAdjacencyList(false, 5, 6);
        graph.AddEdgeToAdjacencyList(false, 5, 3);
        graph.AddEdgeToAdjacencyList(false, 2, 6);
        graph.AddEdgeToAdjacencyList(false, 2, 3);

        // pick vertices with smaller indices first
        for (int i = 0; i < (int)graph.AdjacencyList.size(); i++) {
            sort(graph.AdjacencyList[i].begin(), graph.AdjacencyList[i].end(), GraphVertex::IsIndexSmaller);
        }

        GraphEulerPathHandler handler;

        handler.EulerianPathAlgorithm(graph);

        string str;
        for (int i = 0; i < (int)handler.Path.size(); i++) 
        {
            if (i > 0) {
                str += ' ';
            }

            str += '1' + handler.Path[i];
        }

        printf("EulerianPathAlgorithm: %s, the answer is %s\n", (str == "1 5 7 6 4 3 7 2 6 5 2 4 1") ? "Correct" : "Wrong", str.c_str());

    }

    //
    // reference: USACO Training, section 3.3, Riding the Fences
    //
    static void Test_EulerianPath_6()
    {
        Graph graph;

        graph.V = 7;
        graph.E = 9;

        graph.AdjacencyList.assign(graph.V, vector< GraphVertex >());

        graph.AddEdgeToAdjacencyList(false, 0, 1);
        graph.AddEdgeToAdjacencyList(false, 1, 2);
        graph.AddEdgeToAdjacencyList(false, 2, 3);
        graph.AddEdgeToAdjacencyList(false, 3, 1);
        graph.AddEdgeToAdjacencyList(false, 3, 4);
        graph.AddEdgeToAdjacencyList(false, 1, 4);
        graph.AddEdgeToAdjacencyList(false, 4, 5);
        graph.AddEdgeToAdjacencyList(false, 4, 6);
        graph.AddEdgeToAdjacencyList(false, 3, 5);

        // pick vertices with smaller indices first
        for (int i = 0; i < (int)graph.AdjacencyList.size(); i++) {
            sort(graph.AdjacencyList[i].begin(), graph.AdjacencyList[i].end(), GraphVertex::IsIndexSmaller);
        }

        GraphEulerPathHandler handler;

        handler.EulerianPathAlgorithm(graph);

        handler.AdjustPathByLexicalOrder();

        string str;
        for (int i = 0; i < (int)handler.Path.size(); i++)
        {
            if (i > 0) {
                str += ' ';
            }

            str += '1' + handler.Path[i];
        }

        printf("EulerianPathAlgorithm: %s, the answer is %s\n", (str == "1 2 3 4 2 5 4 6 5 7") ? "Correct" : "Wrong", str.c_str());

    }
};


///////////////////////////////////////////////////////////////////////////////////////////


class GraphHamiltonPathHandler
{
public:

	vector<unsigned char> vecVisited;
	vector<int> vecPath;

	int PathCount;

	GraphHamiltonPathHandler()
	{
	}
	
	static void ProcessHamiltonPath(vector<int>& vecPath, bool isHamiltonCycle)
	{
		for (int i : vecPath) {
			cout << i << " ";
		}
		cout << (isHamiltonCycle ? "Cycle" : " Not Cycle") << endl;
	}

	void FindAllHamiltonPaths(Graph& graph)
	{
		vecVisited.assign(graph.V, false);

		vecPath.clear();
		vecPath.reserve(graph.V);

		vecVisited[0] = true;
		vecPath.push_back(0);

		PathCount = 0;

		FindAllHamiltonPaths(graph, 0);
	}

	void FindAllHamiltonPaths(Graph& graph, int vertex)
	{
		int pathSize = (int)vecPath.size();

		if (pathSize == graph.V)
		{
			bool isHamiltonCycle = false;

			int endVertex = vecPath[pathSize - 1];

			vector<GraphVertex>& adjs = graph.AdjacencyList[vecPath[0]];
			int adjSize = (int)adjs.size();

			for (int i = 0; i < adjSize; i++)
			{
				if (adjs[i].v == endVertex)
				{
					isHamiltonCycle = true;
					break;
				}
			}

			PathCount++;

			ProcessHamiltonPath(vecPath, isHamiltonCycle);

			return;
		}

		vector<GraphVertex>& adjs = graph.AdjacencyList[vertex];
		int adjSize = (int)adjs.size();

		for (int i = 0; i < adjSize; i++)
		{
			int adjVertex = adjs[i].v;

			if (! vecVisited[adjVertex] )
			{
				vecVisited[adjVertex] = true;
				vecPath.push_back(adjVertex);

				FindAllHamiltonPaths(graph, adjVertex);

				//
				// un-comment the following if only wants to find the first path
				//
				//if (PathCount > 0) {
				//	return;
				//}

				vecVisited[adjVertex] = false;
				vecPath.pop_back();
			}
		}
	}

	static void Test()
	{
		Test_HamiltonPath_1();

		cout << endl;

		Test_HamiltonPath_2();
	}

	static void Test_HamiltonPath_1()
	{
		Graph graph;

		graph.CreateAdjacencyList(4);

		graph.AddEdgeToAdjacencyList(false, 0, 1);
		graph.AddEdgeToAdjacencyList(false, 0, 2);
		graph.AddEdgeToAdjacencyList(false, 0, 3);
		graph.AddEdgeToAdjacencyList(false, 1, 2);
		graph.AddEdgeToAdjacencyList(false, 1, 3);
		graph.AddEdgeToAdjacencyList(false, 2, 3);

		GraphHamiltonPathHandler handler;

		handler.FindAllHamiltonPaths(graph);

		//
		// answer should be:
		//
		//	0 1 2 3 Cycle
		//	0 1 3 2	Cycle
		//	0 2 1 3	Cycle
		//	0 2 3 1	Cycle
		//	0 3 1 2	Cycle
		//	0 3 2 1	Cycle
	}

	static void Test_HamiltonPath_2()
	{
		Graph graph;

		graph.CreateAdjacencyList(5);

		graph.AddEdgeToAdjacencyList(false, 0, 1);
		graph.AddEdgeToAdjacencyList(false, 0, 2);
		graph.AddEdgeToAdjacencyList(false, 1, 3);
		graph.AddEdgeToAdjacencyList(false, 2, 3);
		graph.AddEdgeToAdjacencyList(false, 3, 4);

		GraphHamiltonPathHandler handler;

		handler.FindAllHamiltonPaths(graph);

		//
		// answer should be:
		//
		// No output since there is no Hamilton path	
		//
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


static void GraphTester()
{
	//GraphTraversalHandler::Test();

	//GraphArticulationHandler::Test();
	
	//GraphWeightHandler::Test();

    //GraphEulerPathHandler::Test();

	GraphHamiltonPathHandler::Test();
}


///////////////////////////////////////////////////////////////////////////////////////////
