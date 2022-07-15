#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>
#include <string>

using namespace std;


///////////////////////////////////////////////////////////////////////////////////////////


#define INF (1e9)
#define EPS (1e-9)
#define PI (acos(-1.0)) 

#define FLOAT_EQUAL(a, b)				(fabs(a-b)<=EPS)
#define FLOAT_EQUAL_ZERO(a)				(fabs(a) <= EPS)
#define FLOAT_GREATER_EQUAL_ZERO(a)		(a > -EPS)
#define FLOAT_LESS_EQUAL_ZERO(a)		(a < EPS)

#define DEG_to_RAD(d)  (d * PI / 180.0)
#define RAD_to_DEG(r)  (r * 180.0 / PI)


///////////////////////////////////////////////////////////////////////////////////////////


class GeoPointInt
{
public:

	int X, Y;

public:

	GeoPointInt(int x = 0, int y = 0)
	{
		this->X= x;
		this->Y= y;
	}

	bool operator<(const GeoPointInt& other) const
	{
		if (this->X != other.X)
		{
			return this->X < other.X;
		}

		return this->Y < other.Y;
	}

	string ToString() const
	{
		string str;

		str += to_string(this->X);
		str += " ";
		str += to_string(this->Y);

		return str;
	}

public:

	static double GetDistance(const GeoPointInt& p1, const GeoPointInt& p2)
	{
		return hypot(p1.X - p2.X, p1.Y - p2.Y);
	}

	static int CrossProduct(const GeoPointInt& p1, const GeoPointInt& p2)
	{
		return p1.X * p2.Y - p1.Y * p2.X;
	}

	static int DotProduct(const GeoPointInt& p1, const GeoPointInt& p2)
	{
		return p1.X*p2.X + p1.Y*p2.Y;
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class GeoPoint
{
public:

	double X, Y;

public:

	GeoPoint(double x = 0, double y = 0)
	{
		this->X = x;
		this->Y = y;
	}

	bool operator==(const GeoPoint& other) const
	{
		return FLOAT_EQUAL(this->X, other.X) && FLOAT_EQUAL(this->Y, other.Y);
	}

	bool operator<(const GeoPoint& other) const
	{
		if (!FLOAT_EQUAL(this->X,other.X))
		{
			return this->X < other.X;
		}
		
		return this->Y < other.Y;
	}

	void RotateCounterClockwise(double degree)
	{
		double rad = DEG_to_RAD(degree);    // multiply degree with PI / 180.0

		double x = this->X * cos(rad) - this->Y * sin(rad);
		double y = this->X * sin(rad) + this->Y * cos(rad);

		this->X = x;
		this->Y = y;
	}

	void RotateClockwise(double degree)
	{
		double rad = DEG_to_RAD(degree);    // multiply degree with PI / 180.0

		double x = this->X * cos(rad) + this->Y * sin(rad);
		double y = -this->X * sin(rad) + this->Y * cos(rad);

		this->X = x;
		this->Y = y;
	}

	string ToString() const
	{
		string str;

		str += to_string(this->X);
		str += " ";
		str += to_string(this->Y);

		return str;
	}

public:

	static double GetDistance(const GeoPoint& p1, const GeoPoint& p2)
	{
		return hypot(p1.X - p2.X, p1.Y - p2.Y);
	}

	static double CrossProduct(const GeoPoint& p1, const GeoPoint& p2)
	{
		return p1.X * p2.Y - p1.Y * p2.X;
	}

	static double DotProduct(const GeoPoint& p1, const GeoPoint& p2)
	{
		return p1.X*p2.X + p1.Y*p2.Y;
	}

	//
	//  1: Counter-Clockwise (left turn)
	//  0: collinear
	// -1: Clockwise (right turn)
	//
	static int Orientation(const GeoPoint& p1, const GeoPoint& pm, const GeoPoint& p2)
	{
		double crossProd = (pm.Y - p1.Y) * (p2.X - pm.X) - (p2.Y - pm.Y) * (pm.X - p1.X);

		if (FLOAT_EQUAL_ZERO(crossProd)) {
			return 0;
		}

		return (crossProd > 0) ? -1 : 1;
	}

public:


	static void Test()
	{
		int orientation = Orientation(GeoPoint(0, 0), GeoPoint(4, 4), GeoPoint(1, 2));
		bool b = (1 == orientation);

		printf("GeoPoint::Orientation, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), orientation);


		orientation = Orientation(GeoPoint(0, 0), GeoPoint(4, 0), GeoPoint(4, 4));
		b = (1 == orientation);

		printf("GeoPoint::Orientation, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), orientation);


		orientation = Orientation(GeoPoint(0, 0), GeoPoint(4, 0), GeoPoint(8, 0));
		b = (0 == orientation);

		printf("GeoPoint::Orientation, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), orientation);


		orientation = Orientation(GeoPoint(0, 4), GeoPoint(4, 4), GeoPoint(4, 0));
		b = (-1 == orientation);

		printf("GeoPoint::Orientation, %s, the answer is %d \n", (b ? "Correct" : "Wrong"), orientation);
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class GeoVector
{
public:

	double XM, YM; // magnitude of X and Y

public:

	GeoVector(double xm = 0, double ym = 0)
	{
		this->XM = xm;
		this->YM = ym;
	}

	static GeoVector FromPoints(const GeoPoint& a, const GeoPoint& b)
	{
		return GeoVector(b.X - a.X, b.Y - a.Y);
	}

	bool operator==(const GeoVector& other) const
	{
		return FLOAT_EQUAL(this->XM, other.XM) && FLOAT_EQUAL(this->YM, other.YM);
	}

	bool operator<(const GeoVector& other) const
	{
		if (!FLOAT_EQUAL(this->XM, other.XM))
		{
			return this->XM < other.XM;
		}

		return this->YM < other.YM;
	}

	static GeoVector Scale(const GeoVector& v, double s)
	{
		return GeoVector(v.XM * s, v.YM * s);
	}

	static GeoPoint Translate(const GeoVector& v, const GeoPoint& p)
	{
		return GeoPoint(p.X + v.XM, p.Y + v.YM);
	}

	static double DotProduct(const GeoVector& v1, const GeoVector& v2)
	{
		return v1.XM * v2.XM + v1.YM * v2.YM;
	}

	static double CrossProduct(const GeoVector& v1, const GeoVector& v2)
	{
		return v1.XM * v2.YM - v1.YM * v2.XM;
	}

	double NormSquare() const
	{
		return this->XM * this->XM + this->YM * this->YM;
	}

	// returns angle aob in rad
	static double GetAngle(const GeoPoint& a, const GeoPoint& o, const GeoPoint& b)
	{ 
		GeoVector oa = FromPoints(o, a);
		
		GeoVector ob = FromPoints(o, b);

		return acos(DotProduct(oa, ob) / sqrt(oa.NormSquare() * ob.NormSquare()));
	}
	
	static bool IsCounterClockwise(const GeoPoint& a, const GeoPoint& b, const GeoPoint& r)
	{
		GeoVector ab = FromPoints(a, b);

		GeoVector ar = FromPoints(a, r);

		double cp = CrossProduct(ab, ar);

		return cp > 0;
	}

	//
	// The turning test.
	// If the cross product is positive/zero/negative,
	// then we know that a -> b -> r is a left turn / collinear / right turn, respectively
	//
	//  1: left turn
	//  0: collinear
	// -1: right turn
	//
	static int GetTurnDirection(const GeoPoint& a, const GeoPoint& b, const GeoPoint& r)
	{
		GeoVector ab = FromPoints(a, b);

		GeoVector ar = FromPoints(a, r);

		double cp = CrossProduct(ab, ar);

		if (FLOAT_EQUAL_ZERO(cp))
		{
			return 0;
		}

		if (cp > 0)
		{
			return 1;
		}

		return -1;
	}

	string ToString() const
	{
		string str;

		str += to_string(this->XM);
		str += " ";
		str += to_string(this->YM);

		return str;
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class GeoLine // ax + by + c = 0
{
public:

	double A, B, C;

public:

	GeoLine(double a = 0, double b = 0, double c = 0)
	{
		this->A = a;
		this->B = b;
		this->C = c;
	}

	bool GetSlopeAndIntercept(double& slope, double& y_intercept) const
	{
		if (IsVertical())
		{
			return false;
		}

		slope = -this->A / this->B;
		y_intercept = -this->C / this->B;

		return true;
	}

	bool GetPerpendicularSlope(double& slope)  const
	{
		if (IsHorizontal())
		{
			return false;
		}

		slope = this->B / this->A;

		return true;
	}

	GeoLine GetPerpendicularLine(const GeoPoint& p) const
	{
		if (IsHorizontal())
		{
			return GeoLine(1, 0, -p.X);
		}

		if (IsVertical())
		{
			return GeoLine(0, 1, -p.Y);
		}

		double slope = -this->A / this->B;

		return GeoLine(1/slope, 1, -( p.X/slope + p.Y) );
	}

	bool IsVertical() const
	{
		return FLOAT_EQUAL_ZERO(this->B);
	}

	bool IsHorizontal() const
	{
		return FLOAT_EQUAL_ZERO(this->A);
	}

	bool operator==(const GeoLine& other) const
	{
		return FLOAT_EQUAL(this->A, other.A) && FLOAT_EQUAL(this->B, other.B) && FLOAT_EQUAL(this->C, other.C);
	}

	bool operator<(const GeoLine& other) const
	{
		if (!FLOAT_EQUAL(this->A, other.A))
		{
			return this->A < other.A;
		}

		if (!FLOAT_EQUAL(this->B, other.B))
		{
			return this->B < other.B;
		}

		return this->C < other.C;
	}

	string ToString() const
	{
		string str;

		str += to_string(this->A);
		str += " ";
		str += to_string(this->B);
		str += " ";
		str += to_string(this->C);

		return str;
	}

public:

	static GeoLine FromPoints(const GeoPoint& p1, const GeoPoint& p2)
	{
		GeoLine line;

		if (FLOAT_EQUAL(p1.X, p2.X))      // vertical line is fine
		{
			line.A = 1.0;
			line.B = 0.0;
			line.C = p1.X;           // default values
		}
		else
		{
			line.A = -(double)(p1.Y - p2.Y) / (p1.X - p2.X);
			line.B = 1.0;              // IMPORTANT: we fix the value of b to 1.0
			line.C = -(double)(line.A * p1.X) - p1.Y;
		}

		return line;
	}

	static bool AreParallel(const GeoLine& l1, const GeoLine& l2)
	{
		return FLOAT_EQUAL(l1.A, l2.A) && FLOAT_EQUAL(l1.B, l2.B);
	}

	static bool AreSame(const GeoLine& l1, const GeoLine& l2)
	{
		return l1 == l2;
	}

	static bool AreIntersected(const GeoLine& l1, const GeoLine& l2, GeoPoint* pIntersectPoint = 0)
	{
		if (AreParallel(l1, l2))
		{
			return false;
		}

		if (pIntersectPoint != 0)
		{
			// solve system of 2 linear algebraic equations with 2 unknowns
			pIntersectPoint->X = (l2.B * l1.C - l1.B * l2.C) / (l2.A * l1.B - l1.A * l2.B);

			// special case: test for vertical line to avoid division by zero
			if (FLOAT_EQUAL(l1.B, 0))
			{
				pIntersectPoint->Y = -(l2.A * pIntersectPoint->X + l2.C);
			}
			else
			{
				pIntersectPoint->Y = -(l1.A * pIntersectPoint->X + l1.C);
			}
		}

		return true;
	}

	static bool AreCollinear(const GeoLine& line, const GeoPoint& point)
	{
		double d = line.A * point.X + line.B * point.Y + line.C;

		return FLOAT_EQUAL_ZERO(d);
	}

	static bool AreCollinear(const GeoPoint& p1, const GeoPoint& p2, const GeoPoint& p3)
	{
		return AreCollinear(FromPoints(p1, p2), p3);
	}

	static double GetShortestDistance(const GeoLine& line, const GeoPoint& sourcePoint, GeoPoint* pIntersectPoint = 0)
	{
		GeoLine pl = line.GetPerpendicularLine(sourcePoint);

		GeoPoint ip;
		AreIntersected(line, pl, &ip);

		if (pIntersectPoint != 0)
		{
			*pIntersectPoint = ip;
		}

		return hypot(sourcePoint.X - ip.X, sourcePoint.Y - ip.Y);
	}

	static bool AreOnSameSide(const GeoLine& line, const GeoPoint& p1, const GeoPoint& p2)
	{
		double s1 = line.A*p1.X + line.B*p1.Y + line.C;
		double s2 = line.A*p2.X + line.B*p2.Y + line.C;

		bool opposite = (s1 < 0 && s2 > 0) || (s1 > 0 && s2 < 0);

		return !opposite;
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class GeoLineSegmentInt
{
public:

	GeoPointInt P1, P2;

public:

	GeoLineSegmentInt(GeoPointInt p1 = GeoPointInt(), GeoPointInt p2 = GeoPointInt())
	{
		this->P1 = p1;
		this->P2 = p2;
	}

	//
	// Normalize the coordinates so that (X1,Y1) is always the left, and (X2, Y2) is the right
	//
	void Normalize()
	{
		if (this->P1.X > this->P2.X)
		{
			GeoPointInt p = this->P1;

			this->P1 = this->P2;

			this->P2 = p;
		}
	}

	bool IsCollinear(const GeoPointInt& p) const
	{
		if (p.X < min(P1.X, P2.X) || p.X > max(P1.X, P2.X)) {
			return false;
		}

		if (p.Y < min(P1.Y, P2.Y) || p.Y > max(P1.Y, P2.Y)) {
			return false;
		}

		return (p.Y - P1.Y) * (P2.X - P1.X) == (P2.Y - P1.Y) * (p.X - P1.X);
	}

	double GetLength() const
	{
		return hypot(P1.X - P2.X, P1.Y - P2.Y);
	}

	bool IsVertical() const
	{
		return this->P1.X == this->P2.X;
	}

	bool IsHorizontal() const
	{
		return this->P1.Y == this->P2.Y;
	}

	bool operator<(const GeoLineSegmentInt& other) const
	{
		if (this->P1.X != other.P1.X)
		{
			return this->P1.X < other.P1.X;
		}

		if (this->P1.Y != other.P1.Y)
		{
			return this->P1.Y < other.P1.Y;
		}

		if (this->P2.X != other.P2.X)
		{
			return this->P2.X < other.P2.X;
		}

		return this->P2.Y < other.P2.Y;
	}

	string ToString() const
	{
		string str;

		str += to_string(this->P1.X);
		str += " ";
		str += to_string(this->P1.Y);
		str += " ";
		str += to_string(this->P2.X);
		str += " ";
		str += to_string(this->P2.Y);

		return str;
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class GeoLineSegment
{
public:

	GeoPoint P1, P2;

public:

	GeoLineSegment(GeoPoint p1 = GeoPoint(), GeoPoint p2 = GeoPoint())
	{
		this->P1 = p1;
		this->P2 = p2;
	}

	//
	// Normalize the coordinates so that (X1,Y1) is always the left, and (X2, Y2) is the right
	//
	void Normalize()
	{
		if (this->P1.X > this->P2.X)
		{
			GeoPoint p = this->P1;

			this->P1 = this->P2;

			this->P2 = p;
		}
	}

	bool IsCollinear(const GeoPoint& p) const
	{
		if (p.X < min(P1.X, P2.X) || p.X > max(P1.X, P2.X)) {
			return false;
		}

		if (p.Y < min(P1.Y, P2.Y) || p.Y > max(P1.Y, P2.Y)) {
			return false;
		}

		return FLOAT_EQUAL((p.Y - P1.Y) * (P2.X - P1.X), (P2.Y - P1.Y) * (p.X - P1.X));
	}

	double GetLength() const
	{
		return hypot(P1.X - P2.X, P1.Y - P2.Y);
	}

	bool IsVertical() const
	{
		return FLOAT_EQUAL(this->P1.X, this->P2.X);
	}

	bool IsHorizontal() const
	{
		return FLOAT_EQUAL(this->P1.Y, this->P2.Y);
	}

	bool operator==(const GeoLineSegment& other) const
	{
		return this->P1 == other.P1 && this->P2 == other.P2;
	}

	bool operator<(const GeoLineSegment& other) const
	{
		if (!(this->P1 == other.P1))
		{
			return this->P1 < other.P1;
		}

		return this->P2 < other.P2;
	}

	string ToString() const
	{
		string str;

		str += to_string(this->P1.X);
		str += " ";
		str += to_string(this->P1.Y);
		str += " ";
		str += to_string(this->P2.X);
		str += " ";
		str += to_string(this->P2.Y);

		return str;
	}

public:

	static double GetShortestDistance(const GeoLineSegment& seg, const GeoPoint& sourcePoint, GeoPoint* pIntersectPoint = 0)
	{
		GeoLine line = GeoLine::FromPoints(seg.P1, seg.P2);

		GeoLine pl = line.GetPerpendicularLine(sourcePoint);

		GeoPoint intersectPoint;
		GeoPoint& ip = (pIntersectPoint == 0) ? intersectPoint : *pIntersectPoint;

		GeoLine::AreIntersected(line, pl, &ip);

		if (seg.IsCollinear(ip))
		{
			return hypot(sourcePoint.X - ip.X, sourcePoint.Y - ip.Y);
		}

		double d1 = (sourcePoint.X - seg.P1.X)*(sourcePoint.X - seg.P1.X) + (sourcePoint.Y - seg.P1.Y)*(sourcePoint.Y - seg.P1.Y);
		double d2 = (sourcePoint.X - seg.P2.X)*(sourcePoint.X - seg.P2.X) + (sourcePoint.Y - seg.P2.Y)*(sourcePoint.Y - seg.P2.Y);

		if (d1 < d2)
		{
			return sqrt(d1);
		}

		return sqrt(d2);
	}

	static bool AreIntersected(const GeoLineSegment& seg, const GeoLine& line2, GeoPoint* pIntersectPoint = 0)
	{
		GeoPoint ip;

		GeoLine line1 = GeoLine::FromPoints(seg.P1, seg.P2);

		if (!GeoLine::AreIntersected(line1, line2, &ip))
		{
			return false;
		}

		if (!seg.IsCollinear(ip))
		{
			return false;
		}

		if (pIntersectPoint != 0)
		{
			*pIntersectPoint = ip;
		}

		return true;
	}

	static bool AreIntersected(const GeoLineSegment& seg1, const GeoLineSegment& seg2, GeoPoint* pIntersectPoint = 0)
	{
		GeoPoint ip;

		GeoLine line1 = GeoLine::FromPoints(seg1.P1, seg1.P2);
		GeoLine line2 = GeoLine::FromPoints(seg2.P1, seg2.P2);

		if (!GeoLine::AreIntersected(line1, line2, &ip))
		{
			return false;
		}

		if (!seg1.IsCollinear(ip) || !seg2.IsCollinear(ip))
		{
			return false;
		}

		if (pIntersectPoint != 0)
		{
			*pIntersectPoint = ip;
		}

		return true;
	}

	static bool AreOnSameSide(const GeoLineSegment& seg, const GeoPoint& p1, const GeoPoint& p2)
	{
		bool opposite = ((seg.P1.Y - seg.P2.Y)*(p1.X-seg.P1.X)+(seg.P2.X-seg.P1.X)*(p1.Y-seg.P1.Y))*((seg.P1.Y-seg.P2.Y)*(p2.X-seg.P1.X)+(seg.P2.X-seg.P1.X)*(p2.Y-seg.P1.Y)) < 0;
	
		return !opposite;
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class GeoRectangleInt
{
public:
	
	int X1, Y1; // top-left
	int X2, Y2; // bottom-right

public:
	
	GeoRectangleInt(int x1 = 0, int y1 = 0, int x2 = 0, int y2 = 0)
	{
		X1 = x1;
		Y1 = y1;
		X2 = x2;
		Y2 = y2;
	}

	//
	// Normalize the coordinates so that (X1,Y1) is always the top-left, and (X2, Y2) is the bottom-right
	//
	void Normalize()
	{
		if (this->X1 > this->X2)
		{
			int x = this->X1;
			int y = this->Y1;

			this->X1 = this->X2;
			this->Y1 = this->Y2;

			this->X2 = x;
			this->Y2 = y;
		}

		if (this->Y1 > this->Y2)
		{
			int height = abs(this->Y1 - this->Y2);

			this->Y1 -= height;
			this->Y2 += height;
		}
	}

	int GetWidth() const
	{
		return abs(X1 - X2);
	}

	int GetHeight() const
	{
		return abs(Y1 - Y2);
	}

	long long GetArea() const
	{
		return abs((X1 - X2) * (Y1 - Y2));
	}

	double GetDiagnoalDistance()
	{
		return sqrt((X1 - X2)*(X1 - X2) + (Y1 - Y2)*(Y1 - Y2));
	}

	//
	// -1: outside, 0: border, 1: inside 
	//
	int ResolvePosition(int x, int y) const
	{
		GeoRectangleInt rect(*this);

		rect.Normalize();

		if (x < rect.X1 || x > rect.X2 || y < rect.Y1 || y > rect.Y2)
		{
			return -1;
		}

		if (x > rect.X1 && x < rect.X2 && y > rect.Y1 && y < rect.Y2)
		{
			return 1;
		}

		return 0;
	}

	static GeoRectangleInt GetIntersection(const GeoRectangleInt& r1, const GeoRectangleInt& r2)
	{
		//
		// Note this is assume both are already normalized
		//

		int x5 = max(r1.X1, r2.X1);
		int y5 = max(r1.Y1, r2.Y1);
		int x6 = min(r1.X2, r2.X2);
		int y6 = min(r1.Y2, r2.Y2);

		if (x5 > x6 || y5 > y6)
		{
			return GeoRectangleInt();
		}

		return GeoRectangleInt(x5, y5, x6, y6);
	}

	static bool IfOverlap(const GeoRectangleInt &r1, const GeoRectangleInt &r2)
	{
		// The rectangles don't overlap if
		// one rectangle's minimum in some dimension 
		// is greater than the other's maximum in
		// that dimension.

		bool noOverlap = 
			r1.X1 > r2.X2 ||
			r2.X1 > r1.X2 ||
			r1.Y1 > r2.Y2 ||
			r2.Y1 > r1.Y2;

		return !noOverlap;
	}

	static bool IsNumInRect(vector< vector<int> >& grid, const GeoRectangleInt& r, int num)
	{
		for (int x = r.X1; x <= r.X2; x++)
		{
			for (int y = r.Y1; y <= r.Y2; y++)
			{
				if (grid[x][y] == num)
				{
					return true;
				}
			}
		}

		return false;
	}

	static bool IfContain(const GeoRectangleInt &r1, const GeoRectangleInt &r2)
	{
		return r1.X1 <= r2.X1 && r1.Y1 <= r2.Y1 && r1.X2 >= r2.X2 && r1.Y2 >= r2.Y2;
	}

	bool operator<(const GeoRectangleInt& other) const
	{
		if (this->X1 != other.X1)
		{
			return this->X1 < other.X1;
		}

		if (this->Y1 != other.Y1)
		{
			return this->Y1 < other.Y1;
		}

		if (this->X2 != other.X2)
		{
			return this->X2 < other.X2;
		}

		return this->Y2 < other.Y2;
	}

	string ToString() const
	{
		string str;

		str += to_string(this->X1);
		str += " ";
		str += to_string(this->Y1);
		str += " ";
		str += to_string(this->X2);
		str += " ";
		str += to_string(this->Y2);

		return str;
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class GeoRectangle
{
public:

	double X1, Y1; // top-left
	double X2, Y2; // bottom-right

public:

	GeoRectangle(double x1 = 0, double y1 = 0, double x2 = 0, double y2 = 0)
	{
		X1 = x1;
		Y1 = y1;
		X2 = x2;
		Y2 = y2;
	}

	//
	// Normalize the coordinates so that (X1,Y1) is always the top-left, and (X2, Y2) is the bottom-right
	//
	void Normalize()
	{
		if (X1 > X2)
		{
			double x = X1;
			double y = Y1;

			X1 = X2;
			Y1 = Y2;

			X2 = x;
			Y2 = y;
		}

		if (Y1 > Y2)
		{
			double height = abs(Y1 - Y2);

			Y1 -= height;
			Y2 += height;
		}
	}

	double GetWidth() const
	{
		return abs(X1 - X2);
	}

	double GetHeight() const
	{
		return abs(Y1 - Y2);
	}

	double GetArea() const
	{
		return abs((X1 - X2) * (Y1 - Y2));
	}

	double GetDiagnoalDistance()
	{
		return sqrt((X1 - X2)*(X1 - X2) + (Y1 - Y2)*(Y1 - Y2));
	}

	//
	// -1: outside, 0: border, 1: inside  
	//
	int ResolvePosition(double x, double y) const
	{
		GeoRectangle rect(*this);

		rect.Normalize();

		if (x < rect.X1 || x > rect.X2 || y < rect.Y1 || y > rect.Y2)
		{
			return -1;
		}

		if (x > rect.X1 && x < rect.X2 && y > rect.Y1 && y < rect.Y2)
		{
			return 1;
		}

		return 0;
	}

	static GeoRectangle GetIntersection(const GeoRectangle& r1, const GeoRectangle& r2)
	{
		//
		// Note this is assume both are already normalized
		//

		double x5 = max(r1.X1, r2.X1);
		double y5 = max(r1.Y1, r2.Y1);
		double x6 = min(r1.X2, r2.X2);
		double y6 = min(r1.Y2, r2.Y2);

		if (x5 > x6 || y5 > y6)
		{
			return GeoRectangle();
		}

		return GeoRectangle(x5, y5, x6, y6);
	}

	static bool IfOverlap(const GeoRectangle &r1, const GeoRectangle &r2)
	{
		// The rectangles don't overlap if
		// one rectangle's minimum in some dimension 
		// is greater than the other's maximum in
		// that dimension.

		bool noOverlap =
			r1.X1 > r2.X2 ||
			r2.X1 > r1.X2 ||
			r1.Y1 > r2.Y2 ||
			r2.Y1 > r1.Y2;

		return !noOverlap;
	}

	static bool IfContain(const GeoRectangle &r1, const GeoRectangle &r2)
	{
		return r1.X1 <= r2.X1 && r1.Y1 <= r2.Y1 && r1.X2 >= r2.X2 && r1.Y2 >= r2.Y2;
	}

	bool operator==(const GeoRectangle& other) const
	{
		return FLOAT_EQUAL(this->X1, other.X1) &&
			   FLOAT_EQUAL(this->Y1, other.Y1) &&
			   FLOAT_EQUAL(this->X2, other.X2) &&
			   FLOAT_EQUAL(this->Y2, other.Y2);
	}

	bool operator<(const GeoRectangle& other) const
	{
		if (!FLOAT_EQUAL(this->X1,other.X1))
		{
			return this->X1 < other.X1;
		}

		if (!FLOAT_EQUAL(this->Y1,other.Y1))
		{
			return this->Y1 < other.Y1;
		}

		if (!FLOAT_EQUAL(this->X2, other.X2))
		{
			return this->X2 < other.X2;
		}

		return this->Y2 < other.Y2;
	}

	string ToString() const
	{
		string str;

		str += to_string(this->X1);
		str += " ";
		str += to_string(this->Y1);
		str += " ";
		str += to_string(this->X2);
		str += " ";
		str += to_string(this->Y2);

		return str;
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class GeoCircleInt 
{
public:

	int X, Y, Radius;

private:
	
	double Area;

public:
	
	GeoCircleInt(int x = 0, int y = 0, int radius = 0)
	{
		this->X = x;
		this->Y = y;
		this->Radius = max(0, radius);

		this->Area = 0;
	}

	void ComputeArea()
	{
		this->Area = (this->Radius < 0) ? 0 : PI * this->Radius * this->Radius;
	}

	//
	// -1: intersect, 0: tangant, 1: apart 
	//
	int ResolvePosition(const GeoCircleInt& other) const
	{
		int dist = (this->X - other.X)*(this->X - other.X) + (this->Y - other.Y)*(this->Y - other.Y);
		int radius = (this->Radius + other.Radius) * (this->Radius + other.Radius);

		return (dist == radius) ? 0 : ( (dist < radius) ? -1 : 1 );
	}

	//
	// -1: inside, 0: tangant, 1: outside 
	//
	int ResolvePosition(int x, int y) const
	{
		int dist = (this->X - x)*(this->X - x) + (this->Y - y)*(this->Y - y);
		int radius = this->Radius * this->Radius;

		return (dist == radius) ? 0 : ((dist < radius) ? -1 : 1);
	}

	bool operator < (const GeoCircleInt& other) const
	{
		if (this->X != other.X)
		{
			return this->X < other.X;
		}

		if ( this->Y != other.Y)
		{
			return this->Y < other.Y;
		}

		return this->Radius < other.Radius;
	}

	string ToString() const
	{
		string str;

		str += to_string(this->X);
		str += " ";
		str += to_string(this->Y);
		str += " ";
		str += to_string(this->Radius);

		return str;
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class GeoCircle
{
public:

	double X, Y, Radius;

private:

	double Area;

public:

	GeoCircle(double x = 0, double y = 0, double radius = 0)
	{
		this->X = x;
		this->Y = y;
		this->Radius = max(0.0, radius);

		this->Area = 0;
	}

	void ComputeArea()
	{
		this->Area = (this->Radius < 0) ? 0 : PI * this->Radius * this->Radius;
	}

	//
	// -1: intersect, 0: tangant, 1: apart 
	//
	int ResolvePosition(const GeoCircle& other) const
	{
		double dist = (this->X - other.X)*(this->X - other.X) + (this->Y - other.Y)*(this->Y - other.Y);
		double radius = (this->Radius + other.Radius) * (this->Radius + other.Radius);

		return FLOAT_EQUAL(dist, radius) ? 0 : ((dist < radius) ? -1 : 1);
	}

	//
	// -1: inside, 0: tangant, 1: outside 
	//
	int ResolvePosition(double x, double y) const
	{
		double dist = (this->X - x)*(this->X - x) + (this->Y - y)*(this->Y - y);
		double radius = this->Radius * this->Radius;

		return FLOAT_EQUAL(dist,radius) ? 0 : ((dist < radius) ? -1 : 1);
	}

	bool operator == (const GeoCircle& other) const
	{
		return FLOAT_EQUAL(this->X, other.X) && FLOAT_EQUAL(this->Y, other.Y) && FLOAT_EQUAL(this->Radius, other.Radius);
	}

	bool operator < (const GeoCircle& other) const
	{
		if (!FLOAT_EQUAL(this->X, other.X))
		{
			return this->X < other.X;
		}

		if (!FLOAT_EQUAL(this->Y, other.Y))
		{
			return this->Y < other.Y;
		}

		return this->Radius < other.Radius;
	}

	string ToString() const
	{
		string str;

		str += to_string(this->X);
		str += " ";
		str += to_string(this->Y);
		str += " ";
		str += to_string(this->Radius);

		return str;
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class GeoTriangle
{
public:

	GeoPoint P1, P2, P3;

public:

	GeoTriangle(GeoPoint p1 = GeoPoint(), GeoPoint p2 = GeoPoint(), GeoPoint p3 = GeoPoint())
	{
		this->P1 = p1;
		this->P2 = p2;
		this->P3 = p3;
	}

	bool operator==(const GeoLineSegment& other) const
	{
		return this->P1 == other.P1 && this->P2 == other.P2;
	}

	bool operator<(const GeoLineSegment& other) const
	{
		if (!(this->P1 == other.P1))
		{
			return this->P1 < other.P1;
		}

		return this->P2 < other.P2;
	}

	string ToString() const
	{
		string str;

		str += to_string(this->P1.X);
		str += " ";
		str += to_string(this->P1.Y);
		str += " ";
		str += to_string(this->P2.X);
		str += " ";
		str += to_string(this->P2.Y);

		return str;
	}

	double GetPerimeter() const
	{
		double p1p2 = GeoPoint::GetDistance(this->P1, this->P2);
		double p1p3 = GeoPoint::GetDistance(this->P1, this->P3);
		double p2p3 = GeoPoint::GetDistance(this->P2, this->P3);

		double perimeter = (p1p2 + p1p3 + p2p3);

		return perimeter;
	}

	double GetArea() const
	{
		double p1p2 = GeoPoint::GetDistance(this->P1, this->P2);
		double p1p3 = GeoPoint::GetDistance(this->P1, this->P3);
		double p2p3 = GeoPoint::GetDistance(this->P2, this->P3);

		double perimeter = (p1p2 + p1p3 + p2p3);

		double s = perimeter / 2;

		double area = sqrt(s*(s - p1p2)*(s - p1p3)*(s - p2p3));

		return area;
	}

	double GetInscribedCircleRadius() const
	{
		double p1p2 = GeoPoint::GetDistance(this->P1, this->P2);
		double p1p3 = GeoPoint::GetDistance(this->P1, this->P3);
		double p2p3 = GeoPoint::GetDistance(this->P2, this->P3);

		double perimeter = (p1p2 + p1p3 + p2p3);

		double s = perimeter / 2;

		double area = sqrt(s*(s - p1p2)*(s - p1p3)*(s - p2p3));

		double radius = area / (0.5 * perimeter);

		return radius;
	}

	bool GetInscribedCircle(GeoPoint& center, double& radius) const
	{
		double p1p2 = GeoPoint::GetDistance(this->P1, this->P2);
		double p1p3 = GeoPoint::GetDistance(this->P1, this->P3);
		double p2p3 = GeoPoint::GetDistance(this->P2, this->P3);

		double perimeter = (p1p2 + p1p3 + p2p3);

		double s = perimeter / 2;

		double area = sqrt(s*(s - p1p2)*(s - p1p3)*(s - p2p3));

		radius = area / (0.5 * perimeter);

		if (FLOAT_EQUAL_ZERO(radius))
		{
			return false;
		}

		double ratio = p1p2 / p1p3;

		GeoPoint p = GeoVector::Translate(GeoVector::Scale(GeoVector::FromPoints(this->P2, this->P3), ratio / (1 + ratio)), this->P2);

		GeoLine l1 = GeoLine::FromPoints(this->P1, p);

		ratio = p1p2 / p2p3;

		p = GeoVector::Translate(GeoVector::Scale(GeoVector::FromPoints(this->P1, this->P3), ratio / (1 + ratio)), this->P1);

		GeoLine l2 = GeoLine::FromPoints(this->P2, p);

		GeoLine::AreIntersected(l1, l2, &center); // get their intersection point

		return true;
	}

	double GetCircumscribedCircleRadius() const
	{
		double p1p2 = GeoPoint::GetDistance(this->P1, this->P2);
		double p1p3 = GeoPoint::GetDistance(this->P1, this->P3);
		double p2p3 = GeoPoint::GetDistance(this->P2, this->P3);

		double perimeter = (p1p2 + p1p3 + p2p3);

		double s = perimeter / 2;

		double area = sqrt(s*(s - p1p2)*(s - p1p3)*(s - p2p3));

		double radius = (p1p2 * p1p3 * p2p3) / (4 * area);

		return radius;
	}

	bool IsPointInside(const GeoPoint& p) const
	{
		// find another point that's inside (the average of three vertices)
		GeoPoint p2 = GeoPoint((P1.X + P2.X + P3.X) / 3.0, (P1.Y + P2.Y + P3.Y) / 3.0);

		GeoLine l1 = GeoLine::FromPoints(P1, P2);
		if (!GeoLine::AreOnSameSide(l1, p, p2)) {
			return false;
		}

		GeoLine l2 = GeoLine::FromPoints(P1, P3);
		if (!GeoLine::AreOnSameSide(l1, p, p2)) {
			return false;
		}

		GeoLine l3 = GeoLine::FromPoints(P2, P3);
		return GeoLine::AreOnSameSide(l1, p, p2);
	}

public:

	static void Test()
	{
		GeoTriangle tri(GeoPoint(1, 1), GeoPoint(10, 2), GeoPoint(5, 5));

		bool b = tri.IsPointInside(GeoPoint(3, 3));

		printf("GeoTriangle::IsPointInside, %s, the answer is %s \n", (b == true) ? "Correct" : "Wrong", b ? "Yes" : "No");

		b = tri.IsPointInside(GeoPoint(20, 3));

		printf("GeoTriangle::IsPointInside, %s, the answer is %s \n", (b == false) ? "Correct" : "Wrong", b ? "Yes" : "No");
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


struct GeoAngleLessThanCompare
{
	GeoPoint Pivot;

	GeoAngleLessThanCompare(const GeoPoint& pivot)
	{
		Pivot = pivot;
	}

	bool operator()(const GeoPoint& a, const GeoPoint& b)
	{
		//if (GeoPoint::Orientation(this->Pivot, a, b) == 0) // special case
		if (GeoLine::AreCollinear(this->Pivot, a, b)) // special case
		{
			return GeoPoint::GetDistance(this->Pivot, a) < GeoPoint::GetDistance(this->Pivot, b); // check which one is closer
		}

		double d1x = a.X - this->Pivot.X;
		double d1y = a.Y - this->Pivot.Y;
		double d2x = b.X - this->Pivot.X;
		double d2y = b.Y - this->Pivot.Y;

		return (atan2(d1y, d1x) - atan2(d2y, d2x)) < 0;  // compare two angles
	} 

	bool IsLess(const GeoPoint& a, const GeoPoint& b)
	{
		return (*this)(a, b);
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class GeoPolygon
{
public:

	vector<GeoPoint> Vertices; // the input must be such that Vertices[0] = Vertices[n-1]

public:

	GeoPolygon()
	{
	}

	GeoPolygon(const GeoPolygon& other)
	{
		*this = other;
	}

	GeoPolygon& operator=(const GeoPolygon& other)
	{
		if (&other != this)
		{
			this->Vertices = other.Vertices;
		}

		return *this;
	}

	bool operator == (const GeoPolygon& other) const
	{
		if (&other == this)
		{
			return true;
		}

		return this->Vertices == other.Vertices;
	}

	bool operator < (const GeoPolygon& other) const
	{
		return (long long)this < (long long)&other;
	}

	string ToString() const
	{
		string str;

		for (size_t i = 0; i < Vertices.size(); i++)
		{
			str += "(";
			str += Vertices[i].ToString();
			str += ")";
		}

		return str;
	}

public:

	void Clear()
	{
		this->Vertices.clear();
	}

	//
	// returns the perimeter, which is the sum of Euclidian distances
	// of consecutive line segments (polygon edges)
	//
	double GetPerimeter() const
	{
		double result = 0.0;
		for (int i = 0; i < (int)Vertices.size() - 1; i++)  
		{
			result += GeoPoint::GetDistance(Vertices[i], Vertices[i + 1]);
		}

		return result;
	}

	//
	// returns the area, which is half the determinant
	//
	double GetArea() const
	{
		double result = 0.0, x1, y1, x2, y2;

		for (int i = 0; i < (int)Vertices.size() - 1; i++)
		{
			x1 = Vertices[i].X; x2 = Vertices[i + 1].X;
			y1 = Vertices[i].Y; y2 = Vertices[i + 1].Y;
		
			result += (x1 * y2 - x2 * y1);
		}

		return fabs(result) / 2.0;
	}

	//
	// check whether all three consecutive vertices of the polygon form the same turns.
	// returns true if all three consecutive vertices form the same turns.
	//
	bool IsConvex() const
	{ 
		int sz = (int)this->Vertices.size();  

		if (sz <= 3) 
		{ 
			return false; 
		}
		
		bool ccw = GeoVector::IsCounterClockwise(this->Vertices[0], this->Vertices[1], this->Vertices[2]);
		
		for (int i = 1; i < sz - 1; i++)
		{
			if (GeoVector::IsCounterClockwise(this->Vertices[i], this->Vertices[i + 1], this->Vertices[(i + 2) == sz ? 1 : i + 2]) != ccw)
			{
				return false; // different turn -> this polygon is concave
			}
		}
		
		return true;
	}

	bool IsConcave() const
	{
		int sz = (int)this->Vertices.size();

		if (sz <= 3)
		{
			return false;
		}

		return !IsConvex();
	}

	//
	// Winding number algorithm: It works by computing the sum of angles between
	// three points : {P[i], pt, P[i + 1]} where(P[i] - P[i + 1]) are consecutive sides of polygon.
	// If the final sum is 2*PI(360 degrees), then the point is inside polygon.
	//
	bool IsPointInside(const GeoPoint& p) const
	{
		if ((int)this->Vertices.size() == 0)
		{
			return false;
		}

		double sum = 0; // assume the first vertex is equal to the last vertex
		for (int i = 0; i < (int)this->Vertices.size() - 1; i++)
		{
			double angle = GeoVector::GetAngle(this->Vertices[i], p, this->Vertices[i + 1]);

			if (GeoVector::IsCounterClockwise(p, this->Vertices[i], this->Vertices[i + 1]))
			{
				sum += angle; // left turn/ccw
			}
			else 
			{
				sum -= angle; // right turn/cw
			}
		}

		sum = fabs(sum);

		return FLOAT_EQUAL(sum, 2 * PI);
	}	

	static void CutByLine(const GeoPolygon& _poly1, const GeoPoint& a, const GeoPoint& b, GeoPolygon& poly2)
	{
		GeoPolygon poly1 = _poly1;

		poly2.Clear();

		for (int i = 0; i < (int)poly1.Vertices.size(); i++) 
		{
			double left1 = GeoVector::CrossProduct(GeoVector::FromPoints(a, b), GeoVector::FromPoints(a, poly1.Vertices[i]));
			double left2 = 0;

			if (i != (int)poly1.Vertices.size() - 1)
			{ 
				left2 = GeoVector::CrossProduct(GeoVector::FromPoints(a, b), GeoVector::FromPoints(a, poly1.Vertices[i + 1]));
			}
			
			if (left1 > -EPS)
			{
				poly2.Vertices.push_back(poly1.Vertices[i]); // poly1[i] is on the left of ab
			}
			
			if (left1 * left2 < -EPS) // edge (poly1[i], poly1[i+1]) crosses line ab
			{
				GeoPoint intersectPoint;
				GeoLineSegment::AreIntersected(GeoLineSegment(poly1.Vertices[i], poly1.Vertices[i + 1]), GeoLine::FromPoints(a, b), &intersectPoint);

				poly2.Vertices.push_back(intersectPoint);
			}
		}

		if (!poly2.Vertices.empty() && !(poly2.Vertices.back() == poly2.Vertices.front()))
		{
			poly2.Vertices.push_back(poly2.Vertices.front()); // make poly2’s first point = poly2’s last point
		}
	}

	//
	// Ronald Graham’s Scan algorithm O(n log n )
	//
	static void GetConvexHull(vector<GeoPoint>& vertices, GeoPolygon& polygon) // the content of vertices may be reshuffled
	{
		polygon.Clear();

		int i, j, n = (int)vertices.size();
		
		if (n <= 3) 
		{
			if (!(vertices[0] == vertices[n - 1]))
			{
				vertices.push_back(vertices[0]); // safeguard from corner case
			}

			polygon.Vertices = vertices;
			return;
		} 
		
		// special case, the CH is vertices itself
		 // first, find vertices0 = point with lowest Y and if tie: rightmost X
	
		int vertices0 = 0;
		
		for (i = 1; i < n; i++)
		{
			if (vertices[i].Y < vertices[vertices0].Y || (vertices[i].Y == vertices[vertices0].Y && vertices[i].X > vertices[vertices0].X))
			{
				vertices0 = i;
			}
		}
		
		// swap vertices[vertices0] with vertices[0]
		GeoPoint temp = vertices[0]; vertices[0] = vertices[vertices0]; vertices[vertices0] = temp;
		
		// second, sort points by angle w.r.t. pivot vertices0
		GeoAngleLessThanCompare GALTC(vertices[0]);

		sort(++vertices.begin(), vertices.end(), GALTC); // we do not sort vertices[0]
		
		// third, the ccw tests
		polygon.Vertices.push_back(vertices[n - 1]); 
		polygon.Vertices.push_back(vertices[0]); 
		polygon.Vertices.push_back(vertices[1]); // initial S
		
		i = 2; // then, we check the rest
		while (i < n)  // note: N must be >= 3 for this method to work
		{	
			j = (int)polygon.Vertices.size() - 1;
			
			if (GeoVector::IsCounterClockwise(polygon.Vertices[j - 1], polygon.Vertices[j], vertices[i]))
			//if (GeoPoint::Orientation(polygon.Vertices[j - 1], polygon.Vertices[j], vertices[i]) == 1)
			{
				polygon.Vertices.push_back(vertices[i++]); // left turn, accept
			}
			else // or pop the top of S until we have a left turn
			{ 
				polygon.Vertices.pop_back();
			} 
		}   
	} 

public:

	static void Test()
	{
		//	the positions of P6 and P7 w.r.t the polygon:
		//	
		//	7 P5--------------P4
		//	6 |                  \
		//	5 |                    \
		//	4 |   P7                P3
		//	3 |   P1___            /
		//	2 | / P6    \ ___    /
		//	1 P0              P2
		//	0 1 2 3 4 5 6 7 8 9 101112

		GeoPolygon P;

		P.Vertices.push_back(GeoPoint(1, 1));
		P.Vertices.push_back(GeoPoint(3, 3));
		P.Vertices.push_back(GeoPoint(9, 1));
		P.Vertices.push_back(GeoPoint(12, 4));
		P.Vertices.push_back(GeoPoint(9, 7));
		P.Vertices.push_back(GeoPoint(1, 7));
		P.Vertices.push_back(P.Vertices[0]); // loop back

		printf("Perimeter of polygon = %.2lf\n", P.GetPerimeter()); // 31.64
		printf("Area of polygon = %.2lf\n", P.GetArea()); // 49.00
		printf("Is convex = %d\n", P.IsConvex()); // false (P1 is the culprit)

		GeoPoint P6(3, 2); // outside this (concave) polygon
		printf("Point P6 is inside this polygon = %d\n", P.IsPointInside(P6)); // false
		GeoPoint P7(3, 4); // inside this (concave) polygon
		printf("Point P7 is inside this polygon = %d\n", P.IsPointInside(P7)); // true

		//	 cutting the original polygon based on line P[2] -> P[4] (get the left side)
		//	
		//	7 P5--------------P4
		//	6 |               |  \
		//	5 |               |    \
		//	4 |               |     P3
		//	3 |   P1___       |    /
		//	2 | /       \ ___ |  /
		//	1 P0              P2
		//	0 1 2 3 4 5 6 7 8 9 101112
		//	
		//	  new polygon (notice the index are different now):
		//	7 P4--------------P3
		//	6 |               |
		//	5 |               |
		//	4 |               |
		//	3 |   P1___       |
		//	2 | /       \ ___ |
		//	1 P0              P2
		//	0 1 2 3 4 5 6 7 8 9


		GeoPolygon::CutByLine(P, P.Vertices[2], P.Vertices[4], P);

		printf("Perimeter of polygon = %.2lf\n", P.GetPerimeter()); // smaller now 29.15
		printf("Area of polygon = %.2lf\n", P.GetArea()); // 40.00

		//	running convex hull of the resulting polygon (index changes again)
		//
		//	7 P3--------------P2
		//	6 |               |
		//	5 |               |
		//	4 |   P7          |
		//	3 |               |
		//	2 |               |
		//	1 P0--------------P1
		//	0 1 2 3 4 5 6 7 8 9

		vector<GeoPoint> vertices = P.Vertices;

		GeoPolygon::GetConvexHull(vertices, P); // now this is a rectangle
		
		printf("Perimeter of polygon = %.2lf\n", P.GetPerimeter()); // precisely 28.00
		printf("Area of polygon = %.2lf\n", P.GetArea()); // precisely 48.00
		printf("Is convex = %d\n", P.IsConvex()); // true
		printf("Point P6 is inside this polygon = %d\n", P.IsPointInside(P6)); // true
		printf("Point P7 is inside this polygon = %d\n", P.IsPointInside(P7)); // true
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


static void Geometry2DTester()
{
	GeoPoint::Test();

	//GeoTriangle::Test();

	GeoPolygon::Test();
}


///////////////////////////////////////////////////////////////////////////////////////////







