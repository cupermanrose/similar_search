#ifndef _KDTREE_H
#define _KDTREE_H

#include <iostream>     
#include <algorithm>
#include <BLtoDist.h>

using namespace std;

struct KD_node
{
	//double MAX_lat, MAX_lon, MIN_lat, MIN_lon; 
	int split_D;// 0 latitude; 1 longitude
	int node_tra, node_pos; // which trajectory and which position in trajectory
	int ls, rs;
};

vector<KD_node> Start_KD;
vector<KD_node> End_KD;
vector<vector<KD_node>> All_KD;
vector<int> All_KDroot;

bool latitude_cmp(KD_node x, KD_node y) { return (All_Data[x.node_tra].Points[x.node_pos].latitude <All_Data[y.node_tra].Points[y.node_pos].latitude); }
bool longitude_cmp(KD_node x, KD_node y) { return (All_Data[x.node_tra].Points[x.node_pos].longitude <All_Data[y.node_tra].Points[y.node_pos].longitude); }

int build_KDtree(vector<KD_node>* name_KD, int l, int r, int split) {
	if (l > r) return -1;
	int now = (l + r) / 2;
	if (split == 0) { nth_element((*name_KD).begin() + l, (*name_KD).begin() + now, (*name_KD).begin() + r + 1, latitude_cmp); }
	else { nth_element((*name_KD).begin() + l, (*name_KD).begin() + now, (*name_KD).begin() + r + 1, longitude_cmp); }
	
	(*name_KD)[now].split_D = split;
	(*name_KD)[now].ls = build_KDtree(name_KD, l, now - 1, (split + 1) % 2);
	(*name_KD)[now].rs = build_KDtree(name_KD, now + 1, r, (split + 1) % 2);
	//if ((*name_KD)[now].ls != -1) (*name_KD)[(*name_KD)[now].ls].fa = now;
	//if ((*name_KD)[now].rs != -1) (*name_KD)[(*name_KD)[now].rs].fa = now;

	return now;
}

void init_KD() {

	Start_KD.clear(); End_KD.clear();

	for (int i = 0; i < All_Data.size(); i++) {
		KD_node Tempnode;
		Tempnode.node_tra = i; Tempnode.node_pos = 0;
		Start_KD.push_back(Tempnode);
		Tempnode.node_pos = All_Data[i].Points.size() - 1;
		End_KD.push_back(Tempnode);
	}

	root_Start = build_KDtree(&Start_KD, 0, Start_KD.size()-1, 0);
	root_End = build_KDtree(&End_KD, 0, End_KD.size()-1, 0);

	All_KD.clear(); All_KDroot.clear();

	for (int i = 0; i < All_Data.size(); i++) {
		vector<KD_node> TempTree; TempTree.clear();
		KD_node Tempnode;
		for (int j = 0; j < All_Data[i].Points.size(); j++) {
			Tempnode.node_tra = i; Tempnode.node_pos = j;
			TempTree.push_back(Tempnode);
		}
		All_KD.push_back(TempTree);
		All_KDroot.push_back(build_KDtree(&All_KD[i], 0, All_KD[i].size() - 1, 0));
	}
}

void Range_KDsearch(vector<KD_node>* name_KD,set<int>* answer, double Qlatitude, double Qlongitude, int root) {
	Point Qpoint;
	toRad(&Qpoint, Qlatitude, Qlongitude);
	vector<int> SearchQueue;
	SearchQueue.clear();
	SearchQueue.push_back(root);
	while (!SearchQueue.empty()) {
		int now = SearchQueue.back();
		SearchQueue.pop_back();
		if (now == -1) continue;
		int TempTra = (*name_KD)[now].node_tra;
		int TempPos = (*name_KD)[now].node_pos;

		if (bool_dist(&Qpoint,&(All_Data[TempTra].Points[TempPos]))) (*answer).insert(TempTra);

		if ((*name_KD)[now].split_D == 0) {
			Point TempPoint;
			toRad(&TempPoint, All_Data[TempTra].Points[TempPos].latitude, Qlongitude);
			if (bool_dist(&Qpoint, &TempPoint)) {
				SearchQueue.push_back((*name_KD)[now].ls);
				SearchQueue.push_back((*name_KD)[now].rs);
			}
			else {
				if (Qlatitude < All_Data[TempTra].Points[TempPos].latitude) SearchQueue.push_back((*name_KD)[now].ls);
				else SearchQueue.push_back((*name_KD)[now].rs);
			}
		}
		else {
			Point TempPoint;
			toRad(&TempPoint, Qlatitude, All_Data[TempTra].Points[TempPos].longitude);
			if (bool_dist(&Qpoint, &TempPoint)) {
				SearchQueue.push_back((*name_KD)[now].ls);
				SearchQueue.push_back((*name_KD)[now].rs);
			}
			else {
				if (Qlongitude < All_Data[TempTra].Points[TempPos].longitude) SearchQueue.push_back((*name_KD)[now].ls);
				else SearchQueue.push_back((*name_KD)[now].rs);
			}
		}
	}
}

bool Range_KDsearch_forband(vector<KD_node>* name_KD, double Qlatitude, double Qlongitude, int root) {
	Point Qpoint;
	toRad(&Qpoint, Qlatitude, Qlongitude);
	vector<int> SearchQueue;
	SearchQueue.clear();
	SearchQueue.push_back(root);
	while (!SearchQueue.empty()) {
		int now = SearchQueue.back();
		SearchQueue.pop_back();
		if (now == -1) continue;
		int TempTra = (*name_KD)[now].node_tra;
		int TempPos = (*name_KD)[now].node_pos;

		if (bool_dist(&Qpoint, &(All_Data[TempTra].Points[TempPos]))) return true;

		if ((*name_KD)[now].split_D == 0) {
			Point TempPoint;
			toRad(&TempPoint, All_Data[TempTra].Points[TempPos].latitude, Qlongitude);
			if (bool_dist(&Qpoint, &TempPoint)) {
				SearchQueue.push_back((*name_KD)[now].ls);
				SearchQueue.push_back((*name_KD)[now].rs);
			}
			else {
				if (Qlatitude < All_Data[TempTra].Points[TempPos].latitude) SearchQueue.push_back((*name_KD)[now].ls);
				else SearchQueue.push_back((*name_KD)[now].rs);
			}
		}
		else {
			Point TempPoint;
			toRad(&TempPoint, Qlatitude, All_Data[TempTra].Points[TempPos].longitude);
			if (bool_dist(&Qpoint, &TempPoint)) {
				SearchQueue.push_back((*name_KD)[now].ls);
				SearchQueue.push_back((*name_KD)[now].rs);
			}
			else {
				if (Qlongitude < All_Data[TempTra].Points[TempPos].longitude) SearchQueue.push_back((*name_KD)[now].ls);
				else SearchQueue.push_back((*name_KD)[now].rs);
			}
		}
	}

	return false;
}

#endif