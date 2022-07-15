

#include <stdio.h>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cstring>
#include <set>

using namespace std;

//
// Reference: USACO Gold 2018-12 Fine Dining
//

namespace Quick_Dijkstra
{
	struct Edge
	{
		int v_from;
		int v_to;
		int weight;

		Edge(int f = 0, int t = 0, int w = 0) {
			this->v_from = f;
			this->v_to = t;
			this->weight = w;
		}

		bool operator < (const Edge& other) const {
			return this->weight < other.weight;
		}
	};


	struct Graph
	{
		vector< vector<Edge> > table;

		Graph(int n)
		{
			table.assign(n, vector<Edge>());
		}

		void AddUndirectedEdge(int from, int to, int weight) {
			table[from].push_back(Edge(from, to, weight));
			table[to].push_back(Edge(to, from, weight));
		}

		void AddDirectedEdge(int from, int to, int weight) {
			table[from].push_back(Edge(from, to, weight));
		}
	};


	const int INF = 0x7FFFFFFF;

	void Dijkstra(Graph& graph, int source, vector<int>& dist, vector<int>& parent)
	{
		parent.assign(graph.table.size() + 1, -1);

		dist.assign(graph.table.size() + 1, INF);

		dist[source] = 0;

		set<pair<int, int>> pq; // this makes it faster

		pq.insert(make_pair(0, source));

		while (!pq.empty())
		{
			int i = pq.begin()->second;
			pq.erase(pq.begin());

			for (int k = 0; k < (int)graph.table[i].size(); k++)
			{
				Edge& e = graph.table[i][k];
				int j = e.v_to;
				int dist_i_j = e.weight;
				int new_dist = dist[i] + dist_i_j;

				if (dist[j] > new_dist)
				{
					dist[j] = new_dist;

					parent[j] = i;

					pq.insert(make_pair(dist[j], j));
				}
				//
				// This is to record the path, in lexicographical order, so that 
				// a cartain path is chosen when there are multiple solutions.
				// Remove this if path selection is not needed.
				//
				else if (dist[j] == new_dist)
				{
					parent[j] = min(i, parent[j]);
				}
			}
		}
	}
};


