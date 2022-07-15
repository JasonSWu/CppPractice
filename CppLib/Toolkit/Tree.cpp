#include <algorithm>
#include <cstdio>
#include <vector>
#include <queue>
#include <stack>
#include <map>

using namespace std;


///////////////////////////////////////////////////////////////////////////////////////////

struct SegmentTreeNode
{
	int min_idx, max_idx, sum;

	SegmentTreeNode(int _min_idx = 0, int _max_idx = 0, int _sum = 0)
	{
		min_idx = _min_idx;
		max_idx = _max_idx;
		sum = _sum;
	}
};

class SegmentTree
{
private:

	int N;
	vector<int> Points;
	vector<SegmentTreeNode> Memo;

public:

	SegmentTree(const vector<int>& arr)
	{
		N = (int)arr.size();
		Points = arr;     
		Memo.assign(4 * N, 0);	// create large enough vector of zeroes

		Build(1, 0, N - 1);     // recursive Build
	}

	SegmentTree(const int* p, int n)
	{
		N = n;
		Points.assign(p, p + n);
		Memo.assign(4 * N, 0);	// create large enough vector of zeroes

		Build(1, 0, N - 1);     // recursive Build
	}

	int GetRangeMin(int i, int j)
	{
		return GetRangeMin(1, 0, N - 1, i, j);   
	}

	int GetRangeMax(int i, int j)
	{
		return GetRangeMax(1, 0, N - 1, i, j);
	}

	int GetRangeSum(int i, int j)
	{
		return GetRangeSum(1, 0, N - 1, i, j);
	}

	SegmentTreeNode UpdatePoint(int idx, int new_value)
	{
		return UpdatePoint(1, 0, N - 1, idx, new_value);
	}

private:

	int Left(int parent) { return parent << 1; }     // same as binary heap operations
	int Right(int parent) { return (parent << 1) + 1; }

	void Build(int parent, int L, int R)  // O(N log N)
	{
		if (L == R) // as L == R, either one is fine
		{
			Memo[parent] = SegmentTreeNode(L, L, Points[L]);
		}
		else // recursively compute the values
		{
			int left = Left(parent);
			int right = Right(parent);

			Build(left, L, (L + R) / 2);
			Build(right, (L + R) / 2 + 1, R);

			int i1 = Memo[left].min_idx;
			int i2 = Memo[right].min_idx;

			Memo[parent].min_idx = (Points[i1] <= Points[i2]) ? i1 : i2;

			i1 = Memo[left].max_idx;
			i2 = Memo[right].max_idx;

			Memo[parent].max_idx = (Points[i1] >= Points[i2]) ? i1 : i2;

			Memo[parent].sum = Memo[left].sum + Memo[right].sum;
		}
	}

	// Range Min Query
	int GetRangeMin(int parent, int L, int R, int i, int j)  // O(log N)
	{
		if (i >  R || j <  L) return -1;					// current segment outside query range
		if (L >= i && R <= j) return Memo[parent].min_idx;  // current segment inside query range

		// compute the min position in the Left and Right part of the interval
		int left = Left(parent);
		int right = Right(parent);
		int p1 = GetRangeMin(left, L, (L + R) / 2, i, j);
		int p2 = GetRangeMin(right, (L + R) / 2 + 1, R, i, j);

		if (p1 == -1) return p2;   // if we try to access segment outside query
		if (p2 == -1) return p1;   // same as above

		return (Points[p1] <= Points[p2]) ? p1 : p2;
	}         

	// Range Max Query
	int GetRangeMax(int parent, int L, int R, int i, int j)  // O(log N)
	{
		if (i >  R || j <  L) return -1;					// current segment outside query range
		if (L >= i && R <= j) return Memo[parent].max_idx;  // current segment inside query range

		// compute the min position in the Left and Right part of the interval
		int left = Left(parent);
		int right = Right(parent);
		int p1 = GetRangeMax(left, L, (L + R) / 2, i, j);
		int p2 = GetRangeMax(right, (L + R) / 2 + 1, R, i, j);

		if (p1 == -1) return p2;   // if we try to access segment outside query
		if (p2 == -1) return p1;   // same as above

		return (Points[p1] >= Points[p2]) ? p1 : p2;
	}

	// Range Sum Query
	int GetRangeSum(int parent, int L, int R, int i, int j)  // O(log N)
	{
		if (i >  R || j <  L) return 0;					// current segment outside query range
		if (L >= i && R <= j) return Memo[parent].sum;  // current segment inside query range

		// compute the min position in the Left and Right part of the interval
		int left = Left(parent);
		int right = Right(parent);
		int s1 = GetRangeSum(left, L, (L + R) / 2, i, j);
		int s2 = GetRangeSum(right, (L + R) / 2 + 1, R, i, j);

		return s1 + s2;
	}

	SegmentTreeNode UpdatePoint(int parent, int L, int R, int idx, int new_value)
	{
		// this update code is still preliminary, i == j
		// must be able to update range in the future!
		int i = idx, j = idx;

		// if the current interval does not intersect 
		// the update interval, return this Memo node value!
		if (i > R || j < L)
		{
			return Memo[parent];
		}

		// if the current interval is included in the update range,
		// update that Memo[node]
		if (L == i && R == j)
		{
			Points[i] = new_value; // update the underlying array
			return Memo[parent] = SegmentTreeNode(L, L, new_value);
		}

		// compute the minimum pition in the 
		// Left and Right part of the interval
		SegmentTreeNode n1 = UpdatePoint(Left(parent), L, (L + R) / 2, idx, new_value);
		SegmentTreeNode n2 = UpdatePoint(Right(parent), (L + R) / 2 + 1, R, idx, new_value);

		Memo[parent].min_idx = (Points[n1.min_idx] <= Points[n2.min_idx]) ? n1.min_idx : n2.min_idx;
		Memo[parent].max_idx = (Points[n1.max_idx] >= Points[n2.max_idx]) ? n1.max_idx : n2.max_idx;

		Memo[parent].sum = n1.sum + n2.sum;

		return Memo[parent];
	}

public:

	static void Test()
	{
		int arr[] = { 18, 17, 13, 19, 15, 11, 20 };         // the original array

		SegmentTree st(arr, 7);

		printf("           0,  1,  2,  3,  4,  5, 6 \n");
		printf("Points: { 18, 17, 13, 19, 15, 11, 20 } \n\n");
		printf("GetRangeMin(1, 3), index = %d, value = %d \n", st.GetRangeMin(1, 3), arr[st.GetRangeMin(1, 3)]);    
		printf("GetRangeMin(4, 6), index = %d, value = %d \n", st.GetRangeMin(4, 6), arr[st.GetRangeMin(4, 6)]);    
		printf("GetRangeMin(3, 4), index = %d, value = %d \n", st.GetRangeMin(3, 4), arr[st.GetRangeMin(3, 4)]);    
		printf("GetRangeMin(0, 0), index = %d, value = %d \n", st.GetRangeMin(0, 0), arr[st.GetRangeMin(0, 0)]);    
		printf("GetRangeMin(0, 1), index = %d, value = %d \n", st.GetRangeMin(0, 1), arr[st.GetRangeMin(0, 1)]);    
		printf("GetRangeMin(0, 6), index = %d, value = %d \n", st.GetRangeMin(0, 6), arr[st.GetRangeMin(0, 6)]);    

		printf("\n");
		printf("GetRangeMax(1, 3), index = %d, value = %d \n", st.GetRangeMax(1, 3), arr[st.GetRangeMax(1, 3)]);             
		printf("GetRangeMax(4, 6), index = %d, value = %d \n", st.GetRangeMax(4, 6), arr[st.GetRangeMax(4, 6)]);             
		printf("GetRangeMax(3, 4), index = %d, value = %d \n", st.GetRangeMax(3, 4), arr[st.GetRangeMax(3, 4)]);             
		printf("GetRangeMax(0, 0), index = %d, value = %d \n", st.GetRangeMax(0, 0), arr[st.GetRangeMax(0, 0)]);             
		printf("GetRangeMax(0, 1), index = %d, value = %d \n", st.GetRangeMax(0, 1), arr[st.GetRangeMax(0, 1)]);             
		printf("GetRangeMax(0, 6), index = %d, value = %d \n", st.GetRangeMax(0, 6), arr[st.GetRangeMax(0, 6)]);             

		printf("\n");
		printf("GetRangeSum(1, 3), sum = %d \n", st.GetRangeSum(1, 3));
		printf("GetRangeSum(4, 6), sum = %d \n", st.GetRangeSum(4, 6));
		printf("GetRangeSum(3, 4), sum = %d \n", st.GetRangeSum(3, 4));
		printf("GetRangeSum(0, 0), sum = %d \n", st.GetRangeSum(0, 0));
		printf("GetRangeSum(0, 1), sum = %d \n", st.GetRangeSum(0, 1));
		printf("GetRangeSum(0, 6), sum = %d \n", st.GetRangeSum(0, 6));


		arr[5] = 100;
		st.UpdatePoint(5, 100);  

		printf("\n");
		printf("           0,  1,  2,  3,  4,  5,  6 \n");
		printf("Points: { 18, 17, 13, 19, 15, 100, 20 } \n\n");

		printf("GetRangeMin(1, 3), index = %d, value = %d \n", st.GetRangeMin(1, 3), arr[st.GetRangeMin(1, 3)]);             
		printf("GetRangeMin(4, 6), index = %d, value = %d \n", st.GetRangeMin(4, 6), arr[st.GetRangeMin(4, 6)]);             
		printf("GetRangeMin(3, 4), index = %d, value = %d \n", st.GetRangeMin(3, 4), arr[st.GetRangeMin(3, 4)]);             
		printf("GetRangeMin(0, 0), index = %d, value = %d \n", st.GetRangeMin(0, 0), arr[st.GetRangeMin(0, 0)]);             
		printf("GetRangeMin(0, 1), index = %d, value = %d \n", st.GetRangeMin(0, 1), arr[st.GetRangeMin(0, 1)]);             
		printf("GetRangeMin(0, 6), index = %d, value = %d \n", st.GetRangeMin(0, 6), arr[st.GetRangeMin(0, 6)]);             

		printf("\n");
		printf("GetRangeMax(1, 3), index = %d, value = %d \n", st.GetRangeMax(1, 3), arr[st.GetRangeMax(1, 3)]);
		printf("GetRangeMax(4, 6), index = %d, value = %d \n", st.GetRangeMax(4, 6), arr[st.GetRangeMax(4, 6)]);
		printf("GetRangeMax(3, 4), index = %d, value = %d \n", st.GetRangeMax(3, 4), arr[st.GetRangeMax(3, 4)]);
		printf("GetRangeMax(0, 0), index = %d, value = %d \n", st.GetRangeMax(0, 0), arr[st.GetRangeMax(0, 0)]);
		printf("GetRangeMax(0, 1), index = %d, value = %d \n", st.GetRangeMax(0, 1), arr[st.GetRangeMax(0, 1)]);
		printf("GetRangeMax(0, 6), index = %d, value = %d \n", st.GetRangeMax(0, 6), arr[st.GetRangeMax(0, 6)]);
	
		printf("\n");
		printf("GetRangeSum(1, 3), sum = %d \n", st.GetRangeSum(1, 3));
		printf("GetRangeSum(4, 6), sum = %d \n", st.GetRangeSum(4, 6));
		printf("GetRangeSum(3, 4), sum = %d \n", st.GetRangeSum(3, 4));
		printf("GetRangeSum(0, 0), sum = %d \n", st.GetRangeSum(0, 0));
		printf("GetRangeSum(0, 1), sum = %d \n", st.GetRangeSum(0, 1));
		printf("GetRangeSum(0, 6), sum = %d \n", st.GetRangeSum(0, 6));
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


// ref: https://www.hackerearth.com/practice/data-structures/advanced-data-structures/fenwick-binary-indexed-trees/tutorial/

#define LEAST_SIGNIFICANT_BIT(x) (x & (-x))

class FenwickTree // also called BIT, Binary Indexed Tree
{
public:

	vector<int> BIT;
	int* pBIT;

public:

	FenwickTree(int n = 0)
	{
		if (n > 0) {
			Init(n);
		}
	}

	void Init(int n)
	{
		// n + 1 zeroes, ignore index 0
		BIT.assign(n + 1, 0);
		pBIT = &BIT.front();
	}

	FenwickTree(const int* p, int n)
	{
		// n + 1 zeroes, ignore index 0
		BIT.assign(n + 1, 0);
		pBIT = &BIT.front();

		for (int i = 0; i < n; i++) {
			UpdateDelta(i, p[i]);
		}
	}

	int GetRangeSum(int idx) // idx is 0-based index of raw input array
	{
		int index = idx + 1;

		int sum = 0;

		for (; index > 0; index -= LEAST_SIGNIFICANT_BIT(index)) {
			sum += pBIT[index];
		}

		return sum;
	}

	int GetRangeSum(int idxFrom, int idxTo)
	{
		if (idxFrom == 0) {
			return GetRangeSum(idxTo);
		}

		return GetRangeSum(idxTo) - GetRangeSum(idxFrom-1);
	}

	void UpdateDelta(int idx, int delta)
	{
		int index = idx + 1; // 1-base for BIT array

		int size = (int)BIT.size();
		for (; index < size; index += LEAST_SIGNIFICANT_BIT(index))  // note: n = BIT.size() - 1
		{
			pBIT[index] += delta;
		}
	}

	static void Test()
	{
		const int n = 11;
		int array[n] = {  2, 4, 5, 5, 6, 6, 6, 7, 7, 8, 9 };
					   // 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10

		FenwickTree ft(array, n);

		for (int i = 0; i < n; i++) {
			printf("GetRangeSum(%d) = %d\n", i, ft.GetRangeSum(i));
		}

		for (int i = 0; i < n-1; i++) {
			for (int j = i+1; j < n; j++) {
				printf("GetRangeSum(%d, %d) = %d\n", i, j, ft.GetRangeSum(i, j));
			}
		}

		ft.UpdateDelta(2, 100);

		for (int i = 0; i < n; i++) {
			printf("GetRangeSum(%d) = %d\n", i, ft.GetRangeSum(i));
		}
	}

};


class FenwickTree2D // also called BIT, Binary Indexed Tree
{
public:

	vector<vector<int>> BIT;

public:

	FenwickTree2D(int X = 0, int Y = 0)
	{
		// 1-based arrayes, ignore all (0, y) and (x, 0)
		BIT.assign(Y + 1, vector<int>(X + 1, 0));
	}

	FenwickTree2D(const int* p, int X, int Y)
	{
		// 1-based arrayes, ignore all (0, y) and (x, 0)
		BIT.assign(Y + 1, vector<int>(X + 1, 0));

		for (int y = 0; y < Y; y++) {
			for (int x = 0; x < X; x++) {
				UpdateDelta(x, y, *p++);
			}
		}
	}

	int GetRangeSum(int idx_x, int idx_y)
	{
		int index_x = idx_x + 1;
		int index_y = idx_y + 1;

		int sum = 0;

		for (; index_x > 0; index_x -= LEAST_SIGNIFICANT_BIT(index_x)) {
			for (; index_y > 0; index_y -= LEAST_SIGNIFICANT_BIT(index_y)) {
				sum += BIT[index_y][index_x];
			}
		}

		return sum;
	}

	int GetRangeSum(int idxFrom_x, int idxFrom_y, int idxTo_x, int idxTo_y)
	{
		if (idxFrom_x == 0 && idxFrom_y == 0) {
			return GetRangeSum(idxTo_x, idxTo_y);
		}

		return	GetRangeSum(idxTo_x, idxTo_y) -
			GetRangeSum(idxTo_x, idxFrom_y - 1) -
			GetRangeSum(idxFrom_x - 1, idxTo_y) +
			GetRangeSum(idxFrom_x - 1, idxFrom_y - 1);
	}

	void UpdateDelta(int idx_x, int idx_y, int delta)
	{
		int index_x = idx_x + 1;
		int index_y = idx_y + 1;

		int y = index_y;
		for (; index_x < (int)BIT[y].size(); index_x += LEAST_SIGNIFICANT_BIT(index_x))
		{
			index_y = y;
			for (; index_y < (int)BIT.size(); index_y += LEAST_SIGNIFICANT_BIT(index_y))
			{
				BIT[index_y][index_x] += delta;
			}
		}
	}

	static void Test_RangeSum()
	{
		int p[] = { 2, 4, 8, 9,
					4, 6, 7, 5,
					5, 3, 8, 1,
					1, 2, 3, 4 };

		/*
			y
			/\
		 3  |       1 2 3 4
		 2  |       5 3 8 1      Sub-matrix   --->     3 8 1
		 1  |       4 6 7 5      {1,1,3,2}             6 7 5
		 0  |       2 4 8 9
			|
		  --|------ 0 1 2 3 ----> x
			|

		*/


		FenwickTree2D ft(p, 4, 4);

		int sum = ft.GetRangeSum(1, 1);			// 16
		sum = ft.GetRangeSum(0, 0, 2, 2);		// should be 37 but it returns 23 ???
		sum = ft.GetRangeSum(3, 3);				// 72
		sum = ft.GetRangeSum(0, 0);				// 2
		sum = ft.GetRangeSum(1, 1, 3, 2);		// 30
		sum = ft.GetRangeSum(2, 3, 3, 3);		// 7
		sum = ft.GetRangeSum(1, 1, 1, 1);		// 6
	}

};


///////////////////////////////////////////////////////////////////////////////////////////


class BinaryTreeNode
{
public:

	BinaryTreeNode *left, *right;

	int data;
	int position;

	BinaryTreeNode(int data = 0)
	{
		left = right = 0;
		position = 0;
		this->data = data;
	}
};


class BinaryTree
{
public:

	BinaryTreeNode* root;

	BinaryTree(BinaryTreeNode* root = 0)
	{
		this->root = root;
	}

	~BinaryTree()
	{
		Delete();
	}

	void Delete()
	{
		if (root == 0) {
			return;
		}

		queue<BinaryTreeNode*> q;

		q.push(root);

		while (!q.empty())
		{
			BinaryTreeNode* current = q.front();
			q.pop();

			if (current->right != 0) {
				q.push(current->right);
			}

			if (current->left != 0) {
				q.push(current->left);
			}

			delete current;
		}
	}

	//
	// Reference: http://www.geeksforgeeks.org/inorder-tree-traversal-without-recursion/
	//
	void InOrderTraversal(vector<BinaryTreeNode*>& traversal)
	{
		traversal.clear();

		if (root == 0) {
			return;
		}

		stack<BinaryTreeNode*> stack;

		BinaryTreeNode* current = root;

		bool done = false;

		while (!done)
		{
			if (current != 0)
			{
				stack.push(current);
				current = current->left;
			}
			else
			{
				done = stack.empty();

				if (!done)
				{
					current = stack.top();
					stack.pop();

					traversal.push_back(current);

					current = current->right;
				}
			}
		}
	}

	void PreOrderTraversal(vector<BinaryTreeNode*>& traversal)
	{
		traversal.clear();

		if (root == 0) {
			return;
		}

		stack<BinaryTreeNode*> stack;

		stack.push(root);

		while (!stack.empty())
		{
			BinaryTreeNode* current = stack.top();
			stack.pop();

			traversal.push_back(current);

			if (current->right != 0) {
				stack.push(current->right);
			}

			if (current->left != 0) {
				stack.push(current->left);
			}
		}
	}

	//
	// Reference: http://www.geeksforgeeks.org/iterative-postorder-traversal-using-stack/
	//
	void PostOrderTraversal(vector<BinaryTreeNode*>& traversal)
	{
		traversal.clear();

		if (root == 0) {
			return;
		}

		stack<BinaryTreeNode*> stack;

		BinaryTreeNode* current = root;

		do
		{
			while (current != 0)
			{
				// Push root's right child and then root to stack.
				if (current->right != 0) {
					stack.push(current->right);
				}

				stack.push(current);

				// Set root as root's left child  
				current = current->left;
			}

			current = stack.top();
			stack.pop();

			// If the popped item has a right child and the right child is not
			// processed yet, then make sure right child is processed before root
			if (current->right != 0 && ! stack.empty() && stack.top() == current->right)
			{
				stack.pop();				// remove right child from stack
				stack.push(current);		// push root back to stack
				current = current->right;	// change root so that the right child is processed next
			}
			else  
			{
				traversal.push_back(current);
				current = 0;
			}

		} while (!stack.empty());
	}

public:

	static void TestInOrderTraversal()
	{
		BinaryTree bt(new BinaryTreeNode('C'));

		bt.root->left = new BinaryTreeNode('B');
		bt.root->right = new BinaryTreeNode('G');
		bt.root->left->left = new BinaryTreeNode('A');
		bt.root->left->right = new BinaryTreeNode('D');
		bt.root->left->right->left = new BinaryTreeNode('E');
		bt.root->left->right->right = new BinaryTreeNode('F');
		bt.root->right->left = new BinaryTreeNode('H');

		vector<BinaryTreeNode*> traversal;

		bt.InOrderTraversal(traversal);

		string str;
		for (int i = 0; i < (int)traversal.size(); i++) {
			str += (char)traversal[i]->data;
		}

		bool b = str == "ABEDFCHG";

		printf("BinaryTree::InOrderTraversal, %s, the answer is %s \n", (b == true) ? "Correct" : "Wrong", str.c_str());
	}

	static void TestPreOrderTraversal()
	{
		BinaryTree bt(new BinaryTreeNode('C'));

		bt.root->left = new BinaryTreeNode('B');
		bt.root->right = new BinaryTreeNode('G');
		bt.root->left->left = new BinaryTreeNode('A');
		bt.root->left->right = new BinaryTreeNode('D');
		bt.root->left->right->left = new BinaryTreeNode('E');
		bt.root->left->right->right = new BinaryTreeNode('F');
		bt.root->right->left = new BinaryTreeNode('H');

		vector<BinaryTreeNode*> traversal;

		bt.PreOrderTraversal(traversal);

		string str;
		for (int i = 0; i < (int)traversal.size(); i++) {
			str += (char)traversal[i]->data;
		}

		bool b = str == "CBADEFGH";

		printf("BinaryTree::PreOrderTraversal, %s, the answer is %s \n", (b == true) ? "Correct" : "Wrong", str.c_str());
	}

	static void TestPostOrderTraversal()
	{
		BinaryTree bt(new BinaryTreeNode('C'));

		bt.root->left = new BinaryTreeNode('B');
		bt.root->right = new BinaryTreeNode('G');
		bt.root->left->left = new BinaryTreeNode('A');
		bt.root->left->right = new BinaryTreeNode('D');
		bt.root->left->right->left = new BinaryTreeNode('E');
		bt.root->left->right->right = new BinaryTreeNode('F');
		bt.root->right->left = new BinaryTreeNode('H');

		vector<BinaryTreeNode*> traversal;

		bt.PostOrderTraversal(traversal);

		string str;
		for (int i = 0; i < (int)traversal.size(); i++) {
			str += (char)traversal[i]->data;
		}

		bool b = str == "AEFDBHGC";

		printf("BinaryTree::PostOrderTraversal, %s, the answer is %s \n", (b == true) ? "Correct" : "Wrong", str.c_str());
	}

	static void Test()
	{
		TestInOrderTraversal();

		TestPreOrderTraversal();

		TestPostOrderTraversal();
	}
};



///////////////////////////////////////////////////////////////////////////////////////////


class ArrayTreeNode
{
public:

	int parent, self, left, right;
	int data;

	ArrayTreeNode(int data = 0)
	{
		parent = self = left = right = -1;
		this->data = data;
	}
};



class ArrayTree
{
public:

	vector<ArrayTreeNode> nodes;

public:

	ArrayTree(int size = 0)
	{
		nodes.reserve(size);
	}

	int AddNode(int data = 0)
	{
		int size = (int)nodes.size();

		ArrayTreeNode node(data);
		node.self = size;

		nodes.push_back(node);

		return size;
	}

	int AddLeftNode(int parent, int data = 0)
	{
		int idx = AddNode(data);

		if (parent >= 0) {
			nodes[parent].left = idx;
		}

		nodes[idx].parent = parent;

		return idx;
	}

	int AddRightNode(int parent, int data = 0)
	{
		int idx = AddNode(data);

		if (parent >= 0) {
			nodes[parent].right = idx;
		}

		nodes[idx].parent = parent;

		return idx;
	}

	void InOrderTraversal(vector<ArrayTreeNode>& traversal)
	{
		traversal.clear();

		if (nodes.size() == 0) {
			return;
		}

		stack<ArrayTreeNode> stack;

		stack.push(nodes[0]);

		ArrayTreeNode current = nodes[0];

		bool done = false;

		while (!done)
		{
			if (current.self >= 0)
			{
				stack.push(current);
				current = nodes[current.left];
			}
			else
			{
				done = stack.empty();

				if (!done)
				{
					current = stack.top();
					stack.pop();

					traversal.push_back(current);

					current = nodes[current.right];
				}
			}
		}
	}

	void PreOrderTraversal(vector<ArrayTreeNode>& traversal)
	{
		traversal.clear();

		if (nodes.size() == 0) {
			return;
		}

		stack<ArrayTreeNode> stack;

		stack.push(nodes[0]);

		while (!stack.empty())
		{
			ArrayTreeNode current = stack.top();
			stack.pop();

			traversal.push_back(current);

			if (current.right >= 0) {
				stack.push(nodes[current.right]);
			}

			if (current.left >= 0) {
				stack.push(nodes[current.left]);
			}
		}
	}

	void PostOrderTraversal(vector<ArrayTreeNode>& traversal)
	{
		traversal.clear();

		if (nodes.size() == 0) {
			return;
		}

		stack<ArrayTreeNode> stack;

		ArrayTreeNode current = nodes[0];

		do
		{
			while (current.self >= 0)
			{
				// Push root's right child and then root to stack.
				if (current.right >= 0) {
					stack.push(nodes[current.right]);
				}

				stack.push(current);

				// Set root as root's left child  
				if (current.left >= 0) {
					current = nodes[current.left];
				}
				else {
					current.self = -1;
				}
			}

			current = stack.top();
			stack.pop();

			// If the popped item has a right child and the right child is not
			// processed yet, then make sure right child is processed before root
			if (current.right >= 0 && !stack.empty() && stack.top().self == current.right)
			{
				stack.pop();				// remove right child from stack
				stack.push(current);		// push root back to stack
				current = nodes[current.right];	// change root so that the right child is processed next
			}
			else
			{
				traversal.push_back(current);
				current.self = -1;
			}

		} while (!stack.empty());
	}

public:

	static void TestInOrderTraversal()
	{
		ArrayTree at;

		int root = at.AddNode('C');

		int B = at.AddLeftNode(root, 'B');
		int G = at.AddRightNode(root, 'G');

		at.AddLeftNode(B, 'A');

		int D = at.AddRightNode(B, 'D');

		at.AddLeftNode(D, 'E');
		at.AddRightNode(D, 'F');

		at.AddLeftNode(G, 'H');

		vector<ArrayTreeNode> traversal;

		at.PreOrderTraversal(traversal);

		string str;
		for (int i = 0; i < (int)traversal.size(); i++) {
			str += (char)traversal[i].data;
		}

		bool b = str == "ABEDFCHG";

		printf("ArrayTree::InOrderTraversal, %s, the answer is %s \n", (b == true) ? "Correct" : "Wrong", str.c_str());
	}

	static void TestPreOrderTraversal()
	{
		ArrayTree at;

		int root = at.AddNode('C');

		int B = at.AddLeftNode(root, 'B');
		int G = at.AddRightNode(root, 'G');
		
		at.AddLeftNode(B, 'A');

		int D = at.AddRightNode(B, 'D');

		at.AddLeftNode(D, 'E');
		at.AddRightNode(D, 'F');

		at.AddLeftNode(G, 'H');

		vector<ArrayTreeNode> traversal;

		at.PreOrderTraversal(traversal);

		string str;
		for (int i = 0; i < (int)traversal.size(); i++) {
			str += (char)traversal[i].data;
		}

		bool b = str == "CBADEFGH";

		printf("ArrayTree::PreOrderTraversal, %s, the answer is %s \n", (b == true) ? "Correct" : "Wrong", str.c_str());
	}

	static void TestPostOrderTraversal()
	{
		ArrayTree at;

		int root = at.AddNode('C');

		int B = at.AddLeftNode(root, 'B');
		int G = at.AddRightNode(root, 'G');

		at.AddLeftNode(B, 'A');

		int D = at.AddRightNode(B, 'D');

		at.AddLeftNode(D, 'E');
		at.AddRightNode(D, 'F');

		at.AddLeftNode(G, 'H');

		vector<ArrayTreeNode> traversal;

		at.PostOrderTraversal(traversal);

		string str;
		for (int i = 0; i < (int)traversal.size(); i++) {
			str += (char)traversal[i].data;
		}

		bool b = str == "AEFDBHGC";

		printf("ArrayTree::PostOrderTraversal, %s, the answer is %s \n", (b == true) ? "Correct" : "Wrong", str.c_str());
	}

	static void Test()
	{
		TestInOrderTraversal();

		TestPreOrderTraversal();

		TestPostOrderTraversal();
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


static void TreeTester()
{
	SegmentTree::Test();

	//FenwickTree::Test();

	//BinaryTree::Test();

	//ArrayTree::Test();
}


///////////////////////////////////////////////////////////////////////////////////////////

