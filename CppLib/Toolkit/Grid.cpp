#include <stdio.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <stack>
#include <climits>

using namespace std;



///////////////////////////////////////////////////////////////////////////////////////////


typedef enum
{
	GRID_UNVISITED,
	GRID_VISITED,
	GRID_NO_ACCESS,

} GRID_STATE;


typedef enum
{
	NORTH,
	EAST,
	SOUTH,
	WEST,
	//NE,
	//SE,
	//SW,
	//NW,

	NUM_OF_DRECTION

} GridMoveDirection;


struct GridCoor
{
	int x;
	int y;

	GridCoor(int xx = 0, int yy = 0)
	{
		x = xx;
		y = yy;
	}

	bool operator ==(const GridCoor& other) const
	{
		return x == other.x && y == other.y;
	}

	bool operator <(const GridCoor& other) const
	{
		if (x == other.x) {
			return y < other.y;
		}

		return x < other.x;
	}
};


struct GridCoorWeight
{
	int x;
	int y;
	int weight;

	GridCoorWeight(int xx = 0, int yy = 0, int w = 0)
	{
		x = xx;
		y = yy;
		weight = w;
	}

	bool operator ==(const GridCoorWeight& other) const
	{
		return x == other.x && y == other.y && weight == other.weight;
	}

	bool operator <(const GridCoorWeight& other) const
	{
		return weight < other.weight;
	}
};


struct GridMove
{
	int x;
	int y;
	int direct;  // use GridMoveDirection_t

	GridMove(int xx = 0, int yy = 0, int dd = 0)
	{
		x = xx;
		y = yy;
		direct = dd;
	}

	bool operator ==(const GridMove& other) const
	{
		return x == other.x && y == other.y && direct == other.direct;
	}

	bool operator <(const GridMove& other) const
	{
		if (x == other.x) 
		{
			if (y == other.y) {
				return direct < other.direct;
			}

			return y < other.y;
		}

		return x < other.x;
	}
};


class GridCellData
{
public:

	int state;
	int weight;
	int fence_mask;

	GridCellData(int s = GRID_STATE::GRID_UNVISITED, int w = 0, int fm = 0)
	{
		state = s;
		weight = w;
		fence_mask = fm;
	}

	bool operator == (const GridCellData& other) const {
		return state == other.state && weight == other.weight && fence_mask == other.fence_mask;
	}
};



//
// The top-left cell is considered as (0, 0)
//
class Grid
{
public:

	int COL, ROW;

	vector< vector<GridCellData> > Cells;

	int OppositeDirection[NUM_OF_DRECTION];

	vector<GridCoor> MoveDelta;
	vector<char> MoveSign;

public:

	Grid()
	{
		OppositeDirection[NORTH] = SOUTH;
		OppositeDirection[SOUTH] = NORTH;
		OppositeDirection[EAST] = WEST;
		OppositeDirection[WEST] = WEST;
		//OppositeDirection[NE] = SW;
		//OppositeDirection[SW] = NE;
		//OppositeDirection[NW] = SE;
		//OppositeDirection[SE] = NW;

		MoveDelta.resize(8);

		MoveDelta[NORTH] = GridCoor(0, -1);
		MoveDelta[EAST] = GridCoor(1, 0);
		MoveDelta[SOUTH] = GridCoor(0, 1);
		MoveDelta[WEST] = GridCoor(-1, 0);
		//MoveDelta[NE] = GridCoor(1, -1);
		//MoveDelta[SE] = GridCoor(1, 1);
		//MoveDelta[SW] = GridCoor(-1, 1);
		//MoveDelta[NW] = GridCoor(-1, -1);

		MoveSign.resize(8);

		MoveSign[NORTH] = 'N';		
		MoveSign[EAST] = 'E';		
		MoveSign[SOUTH] = 'S';		
		MoveSign[WEST]	= 'W';		
		//MoveSign[NE] = 'O';			
		//MoveSign[SE] = 'P';			
		//MoveSign[SW] = 'Q';			
		//MoveSign[NW] = 'R';			

		Clear();
	}

	Grid(const Grid& other)
	{
		*this = other;
	}

	Grid& operator = (const Grid& other)
	{
		if (this != &other)
		{
			COL = other.COL;
			ROW = other.ROW;

			Cells = other.Cells;

			MoveDelta = other.MoveDelta;
			MoveSign = other.MoveSign;
		}

		return *this;
	}

	bool operator == (const Grid& other) const
	{
		if (ROW != other.ROW || COL != other.COL) {
			return false;
		}

		for (int y = 0; y < ROW; y++) {
			for (int x = 0; x < COL; x++) {
				if ( !(this->Cells[x][y] == other.Cells[x][y])) {
					return false;
				}
			}
		}
		return true;
	}

	void Clear()
	{
		COL = ROW = 0;

		Cells.clear();
	}

	void Create(int R, int C)
	{
		ROW = R;
		COL = C;

		Cells.assign(COL, vector<GridCellData>(ROW));
	}

	void Fill(const GridCellData& cell)
	{
		for (int y = 0; y < ROW; y++) {
			for (int x = 0; x < COL; x++) {
				this->Cells[x][y] = cell;
			}
		}
	}

	void FillWeight(int w)
	{
		for (int y = 0; y < ROW; y++) {
			for (int x = 0; x < COL; x++) {
				Cells[x][y].weight = w;
			}
		}
	}

	void FillState(int s)
	{
		for (int y = 0; y < ROW; y++) {
			for (int x = 0; x < COL; x++) {
				Cells[x][y].state = s;
			}
		}
	}

	bool IsAdjacent(const GridCoor& c1, const GridCoor& c2) const
	{
		for (int d = 0; d < NUM_OF_DRECTION; d++)
		{
			GridCoor c(c2.x + MoveDelta[d].x, c2.y + MoveDelta[d].y);

			if (c1 == c)	{
				return true;
			}
		}

		return false;
	}

	bool IsAnyAdjacent(const GridCoor& c1, set<GridCoor>& cs) const
	{
		for (int d = 0; d < NUM_OF_DRECTION; d++)
		{
			GridCoor c(c1.x + MoveDelta[d].x, c1.y + MoveDelta[d].y);

			if (cs.find(c) != cs.end()) {
				return true;
			}
		}

		return false;
	}

	bool IsPartOf(const GridCoor& c1, set<GridCoor>& cs) const
	{
		return cs.find(c1) != cs.end();
	}

	bool IsValid(const GridCoor& c) const
	{
		if (c.x < 0 || c.x >= COL || c.y < 0 || c.y >= ROW) {
			return false;
		}

		return Cells[c.x][c.y].state != GRID_NO_ACCESS;
	}

	bool CanMove(const GridCoor& currentCoor, int direct) const
	{
		if (Cells[currentCoor.x][currentCoor.y].fence_mask & (1 << direct)) {
			return false;
		}

		GridCoor nextCoor(currentCoor.x + MoveDelta[direct].x, currentCoor.y + MoveDelta[direct].y);

		if (nextCoor.x < 0 || nextCoor.x >= COL || nextCoor.y < 0 || nextCoor.y >= ROW) {
			return false;
		}

		if (Cells[nextCoor.x][nextCoor.y].fence_mask & (1 << OppositeDirection[direct])) {
			return false;
		}

		return Cells[nextCoor.x][nextCoor.y].state != GRID_NO_ACCESS;
	}

	GridCoor MoveStep(const GridCoor& c, int direct)
	{
		return GridCoor(c.x + MoveDelta[direct].x, c.y + MoveDelta[direct].y);
	}

	GridCoor MovePath(const GridCoor& _c, const char* pPath, vector<GridCoor>* pTrail = 0)
	{
		GridCoor c = _c;

		if (pTrail != 0) {
			pTrail->clear();
			pTrail->push_back(c);
		}

		for (int i = 0; i < (int)strlen(pPath); i++)
		{
			int direct = -1;
			for (int d = 0; d < NUM_OF_DRECTION; d++)
			{
				if (MoveSign[d] == pPath[i]) {
					direct = d;
					break;
				}
			}

			c = MoveStep(c, direct);

			if (pTrail != 0) {
				pTrail->push_back(c);
			}
		}

		return c;
	}
	
	int BfsLeastWeight(const GridCoor& start)
	{
		FillWeight(INT_MAX);

		Cells[start.x][start.y].weight = 0;

		queue<GridCoor> q;

		q.push(start);

		while (!q.empty())
		{
			auto current = q.front();
			q.pop();

			auto& currentCell = Cells[current.x][current.y];

			for (int d = 0; d < NUM_OF_DRECTION; d++)
			{
				if (!CanMove(current, d)) {
					continue;
				}

				auto nextCoor = MoveStep(current, d);

				auto& nextCell = Cells[nextCoor.x][nextCoor.y];

				if (nextCell.weight >(currentCell.weight + 1))
				{
					nextCell.weight = (currentCell.weight + 1);
					q.push(nextCoor);
				}
			}
		}

		int maxLeastWeight = 0;
		for (int y = 0; y < ROW; y++) {
			for (int x = 0; x < COL; x++) {
				if (Cells[x][y].weight < INT_MAX) {
					maxLeastWeight = max(maxLeastWeight, Cells[x][y].weight);
				}
			}
		}

		return maxLeastWeight;
	}

	void BFS(const GridCoor& startCoor)
	{
		auto& startCell = Cells[startCoor.x][startCoor.y];

		if (startCell.state == GRID_VISITED || startCell.state == GRID_NO_ACCESS) {
			return;
		}

		queue<GridCoor> q;

		q.push(startCoor);
		startCell.state = GRID_VISITED;

		while (!q.empty())
		{
			auto currentCoor = q.front();
			q.pop();

			auto& currentCell = Cells[currentCoor.x][currentCoor.y];

			//
			// do something on the cell
			//

			for (int d = 0; d < NUM_OF_DRECTION; d++)
			{
				if (!CanMove(currentCoor, d)) {
					continue;
				}

				auto nextCoor = MoveStep(currentCoor, d);

				auto& nextCell = Cells[nextCoor.x][nextCoor.y];

				if (nextCell.state == GRID_UNVISITED)
				{
					q.push(nextCoor);
					nextCell.state = GRID_VISITED;
				}
			}
		}
	}

	void DFS(const GridCoor& startCoor)
	{
		auto& startCell = Cells[startCoor.x][startCoor.y];

		if (startCell.state == GRID_VISITED || startCell.state == GRID_NO_ACCESS) {
			return;
		}

		stack<GridCoor> s;

		s.push(startCoor);
		startCell.state = GRID_VISITED;

		while (!s.empty())
		{
			auto currentCoor = s.top();
			s.pop();

			auto& currentCell = Cells[currentCoor.x][currentCoor.y];

			//
			// do something on the cell
			//

			for (int d = 0; d < NUM_OF_DRECTION; d++)
			{
				if (!CanMove(currentCoor, d)) {
					continue;
				}

				auto nextCoor = MoveStep(currentCoor, d);

				auto& nextCell = Cells[nextCoor.x][nextCoor.y];

				if (nextCell.state == GRID_UNVISITED)
				{
					s.push(nextCoor);
					nextCell.state = GRID_VISITED;
				}
			}
		}
	}

	// Assuming a visited cell can't be visitd again
	void GetLongestStepsByDfs(const GridCoor& startCoor, int direction, int steps, int& longest_steps)
	{
		if (steps > longest_steps) {
			longest_steps = steps;
		}

		auto& startCell = Cells[startCoor.x][startCoor.y];

		if (startCell.state == GRID_VISITED || startCell.state == GRID_NO_ACCESS) {
			return;
		}

		startCell.state = GRID_VISITED;

		bool canMove = CanMove(startCoor, direction);
		GridCoor nextCoor = MoveStep(startCoor, direction);

		if (canMove && Cells[nextCoor.x][nextCoor.y].state == GRID_UNVISITED)
		{
			GetLongestStepsByDfs(nextCoor, direction, steps+1, longest_steps);
		}
		else if (!canMove)
		{
			for (int d = 0; d < NUM_OF_DRECTION; d++)
			{
				if ( d == direction || !CanMove(startCoor, d)) {
					continue;
				}

				nextCoor = MoveStep(startCoor, d);

				if (Cells[nextCoor.x][nextCoor.y].state == GRID_UNVISITED) {
					GetLongestStepsByDfs(nextCoor, d, steps + 1, longest_steps);
				}
			}
		}

		startCell.state = GRID_UNVISITED;
	}

	int GetLongestStepsByDfs()
	{
		int longest_steps = 0;

		//
		// Assuming diagnosal move is prohibited
		//
		GetLongestStepsByDfs(GridCoor(0, 0), EAST, 1, longest_steps);
		GetLongestStepsByDfs(GridCoor(0, 0), SOUTH, 1, longest_steps);

		return longest_steps;
	}

	//
	// Must set the start coor state correctly before calling this method
	//
	void FindEulerianPaths(const GridCoor& currentCoor, const GridCoor& destCoor, int depth, int& count)
	{
		if (currentCoor == destCoor)
		{
			if (depth < ROW * COL - 1) {
				return;
			}

			count++;

			return;
		}

		bool moveFlags[NUM_OF_DRECTION];

		for (int d = 0; d < NUM_OF_DRECTION; d++) {
			moveFlags[d] = CanMove(currentCoor, d);
		}

		// if left and right can go, but up and down can't, isolated areas will be formed
		if (moveFlags[EAST] && moveFlags[WEST] && !moveFlags[NORTH] && !moveFlags[SOUTH]) {
			return;
		}

		// if up and down can go, but left and right can't, isolated areas will be formed
		if (!moveFlags[EAST] && !moveFlags[WEST] && moveFlags[NORTH] && moveFlags[SOUTH]) {
			return;
		}

		for (int d = 0; d < NUM_OF_DRECTION; d++)
		{
			if (!moveFlags[d]) {
				continue;
			}

			auto nextCoor = MoveStep(currentCoor, d);
			auto& nextCell = Cells[nextCoor.x][nextCoor.y];

			if (nextCell.state == GRID_UNVISITED)
			{
				if (!(nextCoor == destCoor))
				{
					int movableWays = 0;
					for (int d2 = 0; d2 < NUM_OF_DRECTION; d2++) {
						if (CanMove(nextCoor, d2))
						{
							auto coor = MoveStep(nextCoor, d2);
							auto& cell = Cells[coor.x][coor.y];

							if (cell.state == GRID_UNVISITED) {
								movableWays++;
							}
						}
					}

					//no way around new cell,so it is an orphan cell
					if (movableWays == 0) {
						return;
					}

					if (movableWays == 1)
					{
						nextCell.state = GRID_VISITED;

						FindEulerianPaths(nextCoor, destCoor, depth + 1, count);

						nextCell.state = GRID_UNVISITED;

						return;
					}
				}

				nextCell.state = GRID_VISITED;

				FindEulerianPaths(nextCoor, destCoor, depth + 1, count);

				nextCell.state = GRID_UNVISITED;
			}
		}
	}

	void CopySection(Grid& section, const GridCoor& ul, const GridCoor& br)
	{
		int width = br.x - ul.x + 1;
		int height = br.y - ul.y + 1;

		section.Create(height, width);

		for (int x = ul.x; x <= br.x; x++) {
			for (int y = ul.y; y <= br.y; y++) {
				section.Cells[x-ul.x][y-ul.y] = Cells[x][y];
			}
		}
	}

	bool IsSectionSame(const GridCoor& ul1, const GridCoor& br1, const GridCoor& ul2, const GridCoor& br2)
	{
		int width1 = br1.x - ul1.x + 1;
		int height1 = br1.y - ul1.y + 1;

		int width2 = br2.x - ul2.x + 1;
		int height2 = br2.y - ul2.y + 1;

		if (width1 != width2 || height1 != height2) {
			return false;
		}

		for (int x = ul1.x; x <= br1.x; x++) {
			for (int y = ul1.y; y <= br1.y; y++) 
			{
				GridCellData& d1 = Cells[x][y];
				GridCellData& d2 = Cells[x - br1.x + br2.x][y - br1.y + br2.y];

				if (memcmp(&d1, &d2, sizeof(GridCellData)) != 0) {
					return false;
				}
			}
		}

		return true;
	}

	void FlipSectionUpDown(const GridCoor& ul, const GridCoor& br)
	{
		int width = br.x - ul.x + 1;
		int height = br.y - ul.y + 1;

		for (int x = ul.x; x <= br.x - width / 2 ; x++)
		{
			int x2 = br.x - (x - ul.x);

			for (int y = ul.y; y <= br.y; y++)
			{
				GridCellData tmp = Cells[x][y];

				Cells[x][y] = Cells[x2][y];

				Cells[x2][y] = tmp;
			}
		}
	}

	void FlipSectionLeftRight(const GridCoor& ul, const GridCoor& br)
	{
		int width = br.x - ul.x + 1;
		int height = br.y - ul.y + 1;

		for (int y = ul.y; y <= br.y - height / 2; y++)
		{
			int y2 = br.y - (y - ul.y);

			for (int x = ul.x; x <= br.x; x++)
			{
				GridCellData tmp = Cells[x][y];

				Cells[x][y] = Cells[x][y2];

				Cells[x][y2] = tmp;
			}
		}
	}

	void FlipUpDown()
	{
		for (int y = 0; y < ROW; y++) {
			for (int x = 0; x < COL / 2; x++)
			{
				GridCellData tmp = Cells[x][y];
				Cells[x][y] = Cells[COL - x - 1][y];
				Cells[COL - x - 1][y] = tmp;
			}
		}
	}

	void FlipLeftRight()
	{
		for (int y = 0; y < ROW / 2; y++) {
			for (int x = 0; x < COL; x++)
			{
				GridCellData tmp = Cells[x][y];
				Cells[x][y] = Cells[x][ROW - y - 1];
				Cells[x][ROW - y - 1] = tmp;
			}
		}
	}

	void FlipLeftRight(Grid& other)
	{
		other.Create(ROW, COL);

		for (int y = 0; y < ROW; y++) {
			for (int x = 0; x < COL; x++) {
				other.Cells[COL - x - 1][y] = Cells[x][y];
			}
		}
	}

	void FlipUpDown(Grid& other)
	{
		other.Create(ROW, COL);

		for (int y = 0; y < ROW; y++) {
			for (int x = 0; x < COL; x++) {
				other.Cells[x][ROW - y - 1] = Cells[x][y];
			}
		}
	}

	void RotateRight90(Grid& other)
	{
		other.Create(COL, ROW);

		for (int y = 0; y < ROW; y++) {
			for (int x = 0; x < COL; x++) {
				other.Cells[ROW-y-1][x] = Cells[x][y];
			}
		}
	}

	void RotateRight180(Grid& other)
	{
		other.Create(ROW, COL);

		for (int y = 0; y < ROW; y++) {
			for (int x = 0; x < COL; x++) {
				other.Cells[COL - x - 1][ROW - y - 1] = Cells[x][y];
			}
		}
	}

	void RotateRight270(Grid& other)
	{
		other.Create(COL, ROW);

		for (int y = 0; y < ROW; y++) {
			for (int x = 0; x < COL; x++) {
				other.Cells[y][COL - x - 1] = Cells[x][y];
			}
		}
	}


public:

	static void Test()
	{
		Test_01();
		Test_02();
		Test_03();
	}

	static void Test_01()
	{
		Grid grid;

		ifstream fin("input_grid_flood_fill.txt");
		ofstream fout("output_grid_flood_fill.txt");

		int COL, ROW;

		fin >> COL >> ROW;

		grid.Create(ROW, COL);

		int num;
		for (int y = 0; y < ROW; y++) {
			for (int x = 0; x < COL; x++) 
			{
				fin >> num;
				grid.Cells[x][y].state = (num == 0) ? GRID_UNVISITED : GRID_NO_ACCESS;
			}
		}

		grid.BfsLeastWeight( GridCoor(0,0) );

		for (int y = 0; y < ROW; y++) {
			for (int x = 0; x < COL; x++)
			{
				auto w = grid.Cells[x][y].weight;

				if (w < 1000) {
					fout << w;
				}
				else {
					fout << 'x';
				}

				fout << ' ';

			}
			fout << endl;
		}

	}

	static void Test_02()
	{
		Grid grid;

		int COL = 5, ROW = 6;

		grid.Create(ROW, COL);

		grid.Cells[2][2].weight = 100;

		grid.FlipSectionUpDown(GridCoor(1, 1), GridCoor(4, 4));
		bool b = grid.Cells[3][2].weight == 100;

		printf("Grid::FlipSectionUpDown, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), grid.Cells[3][2].weight);


		grid.FlipSectionLeftRight(GridCoor(1, 1), GridCoor(4, 4));
		b = grid.Cells[3][3].weight == 100;

		printf("Grid::FlipSectionLeftRight, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), grid.Cells[3][3].weight);


		grid.Cells[1][1].weight = 100;
		b = grid.IsSectionSame(GridCoor(0, 0), GridCoor(1, 1), GridCoor(2, 2), GridCoor(3, 3));

		printf("Grid::IsSectionSame, %s\n", (b ? "Correct" : "Wrong"));


		Grid section;
		grid.CopySection(section, GridCoor(2, 2), GridCoor(3, 3));
		b = grid.Cells[1][1].weight == 100;

		printf("Grid::CopySection, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), grid.Cells[1][1].weight);
	}

	static void Test_03()
	{
		Grid grid;

		int COL = 4, ROW = 3;

		grid.Create(ROW, COL);

		grid.Cells[0][0].weight = 'a'; grid.Cells[1][0].weight = 'b'; grid.Cells[2][0].weight = 'c'; grid.Cells[3][0].weight = 'd';
		grid.Cells[0][1].weight = '1'; grid.Cells[1][1].weight = '2'; grid.Cells[2][1].weight = '3'; grid.Cells[3][1].weight = '4';
		grid.Cells[0][2].weight = 'A'; grid.Cells[1][2].weight = 'B'; grid.Cells[2][2].weight = 'C'; grid.Cells[3][2].weight = 'D';

		Grid other;
		
		grid.RotateRight90(other);
		bool b = other.Cells[0][3].weight == 'D';

		printf("Grid::RotateRight90, %s, the answer is %c \n", (b ? "Correct" : "Wrong"), other.Cells[0][3].weight);


		grid.RotateRight180(other);
		b = other.Cells[3][2].weight == 'a';

		printf("Grid::RotateRight180, %s, the answer is %c \n", (b ? "Correct" : "Wrong"), other.Cells[3][2].weight);


		grid.RotateRight270(other);
		b = other.Cells[0][3].weight == 'a';

		printf("Grid::RotateRight270, %s, the answer is %c \n", (b ? "Correct" : "Wrong"), other.Cells[0][3].weight);


		grid.FlipLeftRight(other);
		b = other.Cells[3][2].weight == 'A';

		printf("Grid::FlipLeftRight, %s, the answer is %c \n", (b ? "Correct" : "Wrong"), other.Cells[3][2].weight);


		grid.FlipUpDown(other);
		b = other.Cells[3][2].weight == 'd';

		printf("Grid::FlipUpDown, %s, the answer is %c \n", (b ? "Correct" : "Wrong"), other.Cells[3][2].weight);
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


static void GridTester()
{
	Grid::Test();
}


///////////////////////////////////////////////////////////////////////////////////////////
