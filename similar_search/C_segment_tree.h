#ifndef _C_SEGMENT_TREE_H
#define _C_SEGMENT_TREE_H

using namespace std;

struct Csegment_node {
	vector<Point> ConvexHull;
	int ls, rs;
};

vector<vector<Csegment_node>> All_Cseg;
vector<int> All_CsegRoot;
vector<Point> TempPoints;
const double Ceps = 1.0E-15;

double dis(Point* A, Point* B) {
	return ((*A).latitude - (*B).latitude)*((*A).latitude - (*B).latitude) + ((*A).longitude - (*B).longitude)*((*A).longitude - (*B).longitude);
}

double multi(Point* A, Point* B, Point* C) { // AC * BC
	return ((*B).latitude - (*A).latitude)*((*C).longitude - (*A).longitude) - ((*C).latitude - (*A).latitude)*((*B).longitude - (*A).longitude);
}

bool CMP(Point A, Point B) { 
	double x = multi(&TempPoints[0], &A, &B);
	if ((x > -Ceps) || ((abs(x-0.0)<Ceps) && (dis(&A, &TempPoints[0]) < dis(&B, &TempPoints[0])))) return 1;
	return 0;
}

void Merge_ConvexHull(vector<Point>* M, vector<Point>* A, vector<Point>* B) {
	TempPoints.clear();
	for (int i = 0; i < (*A).size(); i++) TempPoints.push_back((*A)[i]);
	for (int i = 0; i < (*B).size(); i++) TempPoints.push_back((*B)[i]);

	int k = 0;
	for (int i = 0; i < TempPoints.size(); i++) {
		if ((TempPoints[i].longitude < TempPoints[k].longitude) || ((TempPoints[i].longitude < (TempPoints[k].longitude + Ceps)) && (TempPoints[i].latitude < TempPoints[k].latitude))) k = i;
	}

	swap(TempPoints[0], TempPoints[k]);
	sort(TempPoints.begin() + 1, TempPoints.end(), CMP);
	(*M).push_back(TempPoints[0]);
	(*M).push_back(TempPoints[1]);
	for (int i = 2; i < TempPoints.size(); i++) {
		//multi(A,B,C) AC * BC < 0 C on the left of AB, B invalid
		while (((*M).size() > 1) && (multi(&(*M)[(*M).size() - 2], &(*M)[(*M).size() - 1], &TempPoints[i]) < Ceps)) (*M).pop_back();
		(*M).push_back(TempPoints[i]);
	}
	return;
}

bool InCircle(Point* M, Point* A) {
	//cout << (*M).latitude << ' ' << (*M).longitude <<' ' << (*M).x << ' ' << (*M).y << ' ' << (*M).z << endl;
	//cout << (*A).x << ' ' << (*A).y << ' ' << (*A).z << endl;
	if (!bool_dist(M, A)) return false;
	return true;
}

bool InConvexHull(vector<Point>* M, Point* A) {
	if ((*M).size() < 3) return false;
	for (int i = 1; i < (*M).size(); i++) {
		if (multi(&((*M)[i - 1]), &((*M)[i]), A) < Ceps) return false;
	}
	if (multi(&((*M)[(*M).size() - 1]), &((*M)[0]), A) < Ceps)  return false;
	return true;
}

double Dis_point_segment(Point* start_p, Point* end_p, Point* A) {
	Point AA;
	AA.latitude = (*A).latitude + ((*start_p).longitude - (*end_p).longitude);
	AA.longitude = (*A).longitude + ((*end_p).latitude - (*start_p).latitude);
	if ((multi(A, &AA, start_p))*(multi(A, &AA, end_p)) > Ceps) {
		return min(dis(A, start_p), dis(A, end_p));
	}
	return fabs(multi(start_p,end_p,A)/(dis(start_p,end_p)+Ceps));
}

int OverlapC(Csegment_node* now, Point* Qp) {
	int tot = 0;
	for (int i = 0; i < (*now).ConvexHull.size(); i++) {
		if (InCircle(&((*now).ConvexHull[i]), Qp)) tot++;
	}
	if (tot == (*now).ConvexHull.size()) return 1;// convex hull in circle
	if (tot != 0) return 2;// overlap
	if (InConvexHull(&((*now).ConvexHull), Qp)) return 2;// overlap
	if ((*now).ConvexHull.size() > 1) {
		for (int i = 1; i < (*now).ConvexHull.size(); i++) {
			if (Dis_point_segment(&((*now).ConvexHull[i - 1]), &((*now).ConvexHull[i]), Qp) <= epsilon) return 2;//overlap
		}
		if (Dis_point_segment(&((*now).ConvexHull[(*now).ConvexHull.size()-1]), &((*now).ConvexHull[0]), Qp) <= epsilon) return 2;//overlap
	}
	return 0;// no overlap
}

bool in_CsegSearch(vector<Csegment_node>* T, Point* Qp, int now, int l, int r, int ll) {
	if ((l > r) || (r < ll)) return false;
	int c = OverlapC(&((*T)[now]), Qp);
	if (c == 0) return false;// no overlap
	if (c == 1) { // inside and updade segLen
		segLen = r;
		return true;
	}

	//if (l == r) return true;

	bool flag = true;
	if (ll <= ((l + r) / 2)) flag = in_CsegSearch(T, Qp, (*T)[now].ls, l, (l + r) / 2, ll);
	if (flag) { // if ll>mid or l-mid is all inside
		flag = in_CsegSearch(T, Qp, (*T)[now].rs, (l + r) / 2 + 1, r, ll);
	}
	return flag;
}

bool out_CsegSearch(vector<Csegment_node>* T, Point* Qp, int now, int l, int r, int ll) {
	if ((l > r) || (r < ll)) return false;
	int c = OverlapC(&((*T)[now]), Qp);
	if (c == 1) return false;// inside
	if (c == 0) { // outside and updade segLen
		segLen = r;
		return true;
	}

	//if (l == r) return true;

	bool flag = true;
	if (ll <= ((l + r) / 2)) flag = out_CsegSearch(T, Qp, (*T)[now].ls, l, (l + r) / 2, ll);
	if (flag) { // if ll>mid or l-mid is all outside
		flag = out_CsegSearch(T, Qp, (*T)[now].rs, (l + r) / 2 + 1, r, ll);
	}
	return flag;
}

int build_CsegTree(vector<Csegment_node>* T, vector<Point>* A, int l, int r) {
	if (l > r) return -1;
	Csegment_node Tempnode;
	if (l == r) {
		Tempnode.ConvexHull.clear();
		Tempnode.ConvexHull.push_back((*A)[l]);
		Tempnode.ls = -1;
		Tempnode.rs = -1;
		(*T).push_back(Tempnode);
		return ((*T).size() - 1);
	}

	Tempnode.ls = build_CsegTree(T, A, l, (l + r) / 2);
	Tempnode.rs = build_CsegTree(T, A, (l + r) / 2 + 1, r);
	Tempnode.ConvexHull.clear();
	Merge_ConvexHull(&(Tempnode.ConvexHull), &((*T)[Tempnode.ls].ConvexHull), &((*T)[Tempnode.rs].ConvexHull));
	(*T).push_back(Tempnode);
	return ((*T).size() - 1);
}

void init_Csegment() {
	All_Cseg.clear();
	All_CsegRoot.clear();

	vector<Csegment_node> TempTree;
	for (int i = 0; i < All_Data.size(); i++) {
		TempTree.clear();
		All_CsegRoot.push_back(build_CsegTree(&TempTree, &(All_Data[i].Points), 0, All_Data[i].Points.size() - 1));
		All_Cseg.push_back(TempTree);
	}
}


#endif
