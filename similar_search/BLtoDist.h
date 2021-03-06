#ifndef _BLTODIST_H
#define _BLTODIST_H

#define PI 3.14159265358979323846
#define R 6371.0
#define toRadians(degree) (degree * PI / 180)

const double eps = 1.0E-15;
double epsilon;

struct IntPoint {
	int latitude, longitude;
};

inline int Int_dist(IntPoint& i, IntPoint& j) {

	int lat = i.latitude - j.latitude;
	int lon = i.longitude - j.longitude;
	int c = lat*lat + lon*lon;
	//return sqrt(c);
	return c;
}

bool bool_Intdist(IntPoint& i, IntPoint& j) {
	//double c = (*i).x * (*j).x + (*i).y * (*j).y + (*i).z * (*j).z;
	//if (c > 1.0) c = c - eps;
	//if (c < -1.0) c = c + eps;
	//double a = acos(c);
	////return a*R; //km
	//if ((a*R) > epsilon) return 0;
	//return 1;
	double lat = i.latitude - j.latitude;
	double lon = i.longitude - j.longitude;
	double c = lat*lat + lon*lon;
	//if (sqrt(c) > epsilon) return 0;
	if (c > epsilon) return 0;
	return 1;
}
struct Point { 
	double x, y, z;
	double latitude, longitude;
};

void toRad(Point& p, double latitude, double longitude) {
	p.latitude = latitude;
	p.longitude = longitude;
	/*latitude = toRadians(latitude);
	longitude = toRadians(longitude);
	p.x = cos(latitude)*cos(longitude);
	p.y = cos(latitude)*sin(longitude);
	p.z = sin(latitude);*/
}

bool bool_dist(Point& i, Point& j) {
	//double c = (*i).x * (*j).x + (*i).y * (*j).y + (*i).z * (*j).z;
	//if (c > 1.0) c = c - eps;
	//if (c < -1.0) c = c + eps;
	//double a = acos(c);
	////return a*R; //km
	//if ((a*R) > epsilon) return 0;
	//return 1;
	double lat = i.latitude - j.latitude;
	double lon = i.longitude - j.longitude;
	double c = lat*lat + lon*lon;
	//if (sqrt(c) > epsilon) return 0;
	if (c > epsilon) return 0;
	return 1;
}

inline double double_dist(Point& i, Point& j) {
	//double c = i.x*j.x + i.y*j.y + i.z*j.z;
	//if (c > 1.0) c = c - eps;
	//if (c < -1.0)c = c + eps;
	//double a = acos(c);
	//return a*R; //km
	double lat = i.latitude - j.latitude;
	double lon = i.longitude - j.longitude;
	double c = lat*lat + lon*lon; 
	//return sqrt(c);
	return c;
}

#endif