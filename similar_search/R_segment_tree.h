#ifndef _R_SEGMENT_TREE_H
#define _R_SEGMENT_TREE_H

using namespace std;

struct Rsegment_node {
	Point lt, rt, lb, rb;// left right top bottom
	int ls, rs;
};

vector<vector<Rsegment_node>> All_Rseg;
vector<int> All_RsegRoot;

int build_RsegTree(vector<Rsegment_node>* T, vector<Point>* A, int l, int r) {
	if (l > r) return -1;
	Rsegment_node Tempnode;
	if (l == r) {
		toRad(&(Tempnode.rt), (*A)[l].latitude, (*A)[l].longitude);
		toRad(&(Tempnode.rb), (*A)[l].latitude, (*A)[l].longitude);
		toRad(&(Tempnode.lb), (*A)[l].latitude, (*A)[l].longitude);
		toRad(&(Tempnode.lt), (*A)[l].latitude, (*A)[l].longitude);
		Tempnode.ls = -1;
		Tempnode.rs = -1;
		(*T).push_back(Tempnode);
		return ((*T).size()-1);
	}

	Tempnode.ls = build_RsegTree(T, A, l, (l + r) / 2);
	Tempnode.rs = build_RsegTree(T, A, (l + r) / 2 + 1, r);
	double Max_lat = max((*T)[Tempnode.ls].rt.latitude, (*T)[Tempnode.rs].rt.latitude);
	double Min_lat = min((*T)[Tempnode.ls].lb.latitude, (*T)[Tempnode.rs].lb.latitude);
	double Max_lon = max((*T)[Tempnode.ls].rt.longitude, (*T)[Tempnode.rs].rt.longitude);
	double Min_lon = min((*T)[Tempnode.ls].lb.longitude, (*T)[Tempnode.rs].lb.longitude);
	toRad(&(Tempnode.rt), Max_lat, Max_lon);
	toRad(&(Tempnode.rb), Max_lat, Min_lon);
	toRad(&(Tempnode.lb), Min_lat, Min_lon);
	toRad(&(Tempnode.lt), Min_lat, Max_lon);
	(*T).push_back(Tempnode);
	return ((*T).size()-1);
}

int OverlapR(Rsegment_node* now, Point* Qp) {
	int tot = 0;
	if (bool_dist(&((*now).rt), Qp)) tot++;
	if (bool_dist(&((*now).rb), Qp)) tot++;
	if (bool_dist(&((*now).lb), Qp)) tot++;
	if (bool_dist(&((*now).lt), Qp)) tot++;
	if (tot == 4) return 1;// rectangle in circle
	if (tot > 0) return 2;// overlap

	Point tempp; double lat, lon;
	double eps = 1.0E-10;
	if (((*Qp).latitude > ((*now).lb.latitude - eps)) && ((*Qp).latitude < ((*now).rt.latitude + eps))) { lat = (*Qp).latitude; }
	else {
		if ((*Qp).latitude < ((*now).lb.latitude + eps)) { lat = (*now).lb.latitude; }
		else lat = (*now).rt.latitude;
	}
	if (((*Qp).longitude >((*now).lb.longitude - eps)) && ((*Qp).longitude < ((*now).rt.longitude + eps))) { lon = (*Qp).longitude; }
	else {
		if ((*Qp).longitude < ((*now).lb.longitude + eps)) { lon = (*now).lb.longitude; }
		else lon = (*now).rt.longitude;
	}
	
	toRad(&tempp, lat, lon);
	if (bool_dist(&tempp, Qp)) return 2; // overlap
	else return 0; // no overlap
}

bool in_RsegSearch(vector<Rsegment_node>* T, Point* Qp, int now, int l, int r, int ll) {
	if (now == -1) {
		cout << now << endl;
	}
	if ((l > r) || (r < ll)) return false;
	int c = OverlapR(&((*T)[now]), Qp);
	if (c == 0) return false;// no overlap
	if (c == 1) { // inside and updade segLen
		segLen = r;
		return true;
	} 

	//if (l == r) return true;

	bool flag = true;
	if (ll <= ((l + r) / 2)) flag = in_RsegSearch(T, Qp, (*T)[now].ls, l, (l + r) / 2, ll);
	if (flag) { // if ll>mid or l-mid is all inside
		flag = in_RsegSearch(T, Qp, (*T)[now].rs, (l + r) / 2 + 1, r, ll);
	}
	return flag;
}

bool out_RsegSearch(vector<Rsegment_node>* T, Point* Qp, int now, int l, int r, int ll) {
	if (now == -1) {
		cout << now << endl;
	}
	if ((l > r) || (r < ll)) return false;
	int c = OverlapR(&((*T)[now]), Qp);
	if (c == 1) return false;// inside
	if (c == 0) { // outside and updade segLen
		segLen = r;
		return true;
	}

	//if (l == r) return true;

	bool flag = true;
	if (ll <= ((l + r) / 2)) flag = out_RsegSearch(T, Qp, (*T)[now].ls, l, (l + r) / 2, ll);
	if (flag) { // if ll>mid or l-mid is all outside
		flag = out_RsegSearch(T, Qp, (*T)[now].rs, (l + r) / 2 + 1, r, ll);
	}
	return flag;
}

void init_Rsegment() {
	All_Rseg.clear();
	All_RsegRoot.clear();

	vector<Rsegment_node> TempTree;
	for (int i = 0; i < All_Data.size(); i++) {
		TempTree.clear();
		All_RsegRoot.push_back(build_RsegTree(&TempTree, &(All_Data[i].Points), 0, All_Data[i].Points.size() - 1));
		All_Rseg.push_back(TempTree);
	}
}

#endif // !_R_SEGMENT_TREE_H
