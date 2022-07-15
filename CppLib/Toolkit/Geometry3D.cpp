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


struct PointInt
{
	int x, y, z;

	PointInt(int x = 0, int y = 0, int z = 0)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	double Length()		// relative to origin
	{
		return sqrt(x*x + y*y + z*z);
	}

	static double Distance(const PointInt& u, const PointInt& v)
	{
		return sqrt( (u.x - v.x)*(u.x - v.x) + (u.y - v.y)*(u.y - v.y) + (u.z - v.z)*(u.z - v.z));
	}

	static PointInt CrossProduct(const PointInt& u, const PointInt& v)
	{
		int i = u.y * v.z - u.z * v.y;
		int j = u.z * v.x - u.x * v.z;
		int k = u.x * v.y - u.y * v.x;

		return PointInt(i, j, k);
	}

	static int DotProduct(const PointInt& u, const PointInt& v)
	{
		return u.x*v.x + u.y*v.y + u.z*v.z;
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


struct Point
{
	double x, y, z;

	Point(double x = 0, double y = 0, double z = 0)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	double Length()		// relative to origin
	{
		return sqrt(x*x + y*y + z*z);
	}

	static double Distance(const Point& p1, const Point& p2)
	{
		return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
	}

	static Point CrossProduct(const Point& u, const Point& v)
	{
		double i = u.y * v.z - u.z * v.y;
		double j = u.z * v.x - u.x * v.z;
		double k = u.x * v.y - u.y * v.x;

		return Point(i, j, k);
	}

	static double DotProduct(const Point& u, const Point& v)
	{
		return u.x*v.x + u.y*v.y + u.z*v.z;
	}
};



///////////////////////////////////////////////////////////////////////////////////////////


struct Segment
{
	Point p1, p2;

	Segment(Point p1, Point p2)
	{
		this->p1 = p1;
		this->p2 = p2;
	}

	double Length()		
	{
		return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


struct Line
{
	double a, b, c, d;  // ax + by + cz + d = 0

	Line(double a = 0, double b = 0, double c = 0, double d = 0)
	{
		this->a = a;
		this->b = b;
		this->c = c;
		this->d = d;
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


struct Plane
{
	Point p1, p2, p3;

	Plane(Point p1, Point p2, Point p3)
	{
		this->p1 = p1;
		this->p2 = p2;
		this->p3 = p3;
	}

	double AreaByHeron()
	{
		double dist12 = Point::Distance(p1, p2);
		double dist13 = Point::Distance(p1, p3);
		double dist23 = Point::Distance(p2, p3);

		double s = (dist12 + dist13 + dist23) * 0.5;

		return sqrt(s * (s - dist12) * (s - dist13) * (s - dist23));
	}

public:

	static double AreaByHeron(double dist1, double dist2, double dist3)
	{
		double s = (dist1 + dist2 + dist3) * 0.5;

		return sqrt(s * (s - dist1) * (s - dist2) * (s - dist3));
	}
};


///////////////////////////////////////////////////////////////////////////////////////////


class GeoUtil
{
public:

	static bool AreSegmentsParallel(const Point& u1, const Point& u2, const Point& v1, const Point& v2)
	{
		Point u(u2.x - u1.x, u2.y - u1.y, u2.z - u1.z);
		Point v(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z);

		double cross = Point::CrossProduct(u, v).Length();

		return FLOAT_EQUAL_ZERO(cross);
	}

	static bool ArePointsCoplanar(const Point& p, const Point& p1, const Point& p2, const Point& p3)
	{
		//
		// The four points are in the same plane if and only if one of the areas of the four triangles you can make has an area equal to the sum
		//
		
		double dist12 = Point::Distance(p1, p2);
		double dist13 = Point::Distance(p1, p3);
		double dist23 = Point::Distance(p2, p3);

		double area123 = Plane::AreaByHeron(dist12, dist13, dist23);

		double dist01 = Point::Distance(p, p1);
		double dist02 = Point::Distance(p, p2);
		double dist03 = Point::Distance(p, p3);

		double area012 = Plane::AreaByHeron(dist01, dist02, dist12);
		double area013 = Plane::AreaByHeron(dist01, dist03, dist13);
		double area023 = Plane::AreaByHeron(dist02, dist03, dist23);

		if (FLOAT_EQUAL(area012 + area013 - area023, area123)) {
			return true;
		}

		if (FLOAT_EQUAL(area123 + area013 - area023, area012)) {
			return true;
		}

		if (FLOAT_EQUAL(area012 + area123 - area023, area013)) {
			return true;
		}

		if (FLOAT_EQUAL(area012 + area013 - area123, area023)) {
			return true;
		}

		return false;
	}

	static bool ArePointsCoplanar(const vector<Point>& p)
	{
		int size = (int)p.size();

		if (size < 3) {
			return false;
		}

		if (size == 3) {
			return true;
		}

		for (int i = 3; i < size; i++)
		{
			if (!ArePointsCoplanar(p[i], p[0], p[1], p[2])) {
				return false;
			}
		}

		return true;
	}

public:

	static void Test()
	{
		bool b = AreSegmentsParallel(Point(1, 3, 5), Point(2, 4, 5), Point(2, 1, 5), Point(3, 2, 5));

		printf("GeoUtil::AreSegmentsParallel, %s, the answer is %s \n", (b == true) ? "Correct" : "Wrong", b ? "Yes" : "No");

		b = AreSegmentsParallel(Point(1, 3, 5), Point(2, 4, 5), Point(2, 1, 5), Point(3, 2, 6));

		printf("GeoUtil::AreSegmentsParallel, %s, the answer is %s \n", (b == false) ? "Correct" : "Wrong", b ? "Yes" : "No");

		b = ArePointsCoplanar(Point(4, 4, 0), Point(10, 0, 0), Point(0, 10, 0), Point(10, 10, 0));

		printf("GeoUtil::ArePointsCoplanar, %s, the answer is %s \n", (b == true) ? "Correct" : "Wrong", b ? "Yes" : "No");

		b = ArePointsCoplanar(Point(5, 5, 0), Point(10, 0, 0), Point(0, 10, 0), Point(10, 10, 0));

		printf("GeoUtil::ArePointsCoplanar, %s, the answer is %s \n", (b == true) ? "Correct" : "Wrong", b ? "Yes" : "No");

		b = ArePointsCoplanar(Point(5, 5, 5), Point(10, 0, 0), Point(0, 10, 0), Point(10, 10, 0));

		printf("GeoUtil::ArePointsCoplanar, %s, the answer is %s \n", (b == false) ? "Correct" : "Wrong", b ? "Yes" : "No");

		vector<Point> p;
		p.push_back(Point(10, 0, 0));
		p.push_back(Point(0, 10, 0));
		p.push_back(Point(10, 10, 0));
		p.push_back(Point(5, 5, 0));
		p.push_back(Point(4, 4, 0));
		//p.push_back(Point(3, 4, 0));

		b = ArePointsCoplanar(p);

		printf("GeoUtil::ArePointsCoplanar, %s, the answer is %s \n", (b == true) ? "Correct" : "Wrong", b ? "Yes" : "No");
	}
};


///////////////////////////////////////////////////////////////////////////////////////////



static void Geometry3DTester()
{
	GeoUtil::Test();
}


///////////////////////////////////////////////////////////////////////////////////////////








