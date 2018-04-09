#ifndef _DFD_H
#define _DFD_H
#include <vector>
#include <unordered_set>
#include <grouping.h>

using namespace std;

bool f[2][MaxLen];
double g[2][MaxLen];
bool f_col[MaxLen];	// f_col[i]=false if this column is all 0;
int f_near[MaxLen]; // the nearest 1 on the right;

//struct pairhash {
//public:
//	template <typename T, typename U>
//	std::size_t operator()(const std::pair<T, U> &x) const
//	{
//		return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
//	}
//};
//unordered_set<pair<int, int>, pairhash> rep;

double LB_cell(vector<Point>& A, vector<Point>& B) {
	return (max(double_dist(A[0], B[0]), double_dist(A[A.size() - 1], B[B.size() - 1])));
}

bool LB_band(int Q, int A) { // true: has a LB_band; false no lower bound
	bool flag = false;
	for (int p = 0; p < All_Query[Q].Points.size(); p++) {
		if (!Range_KDsearch_forband(All_KD[A], All_Query[Q].Points[p].latitude, All_Query[Q].Points[p].longitude, All_KDroot[A])) { flag = true; break; }
	}
	if (flag) return true;
	flag = false;
	for (int p = 0; p < All_Data[A].Points.size(); p++) {
		if (!Range_KDsearch_forband(All_KD[Q], All_Data[A].Points[p].latitude, All_Data[A].Points[p].longitude, All_KDroot[Q])) { flag = true; break; }
	}
	if (flag) return true;
	return false;
}

//double LB_cross(vector<Point>* A, vector<Point>* B) {
//	double LB = 0, LB_row, LB_col;
//	int LengthA = (*A).size();
//	int LengthB = (*B).size();
//	//start
//	//i=1 j=0..lengthB
//	LB_row = INFINITE;
//	for (int j = 0; j < LengthB; j++) {
//		LB_row = min(LB_row, double_dist((*A)[1], (*B)[j]));
//	}
//	//i=0..lengthA,j=1
//	LB_col = INFINITE;
//	for (int i = 0; i < LengthA; i++) {
//		LB_col = min(LB_col, double_dist((*A)[i], (*B)[1]));
//	}
//	LB = max(LB_row, LB_col);
//
//	//end
//	//i=LengthA-2 j=0..lengthB
//	LB_row = INFINITE;
//	for (int j = 0; j < LengthB; j++) {
//		LB_row = min(LB_row, double_dist((*A)[LengthA], (*B)[j]));
//	}
//	//i=0..lengthA,j=LengthB-2
//	LB_col = INFINITE;
//	for (int i = 0; i < LengthA; i++) {
//		LB_col = min(LB_col, double_dist((*A)[i], (*B)[LengthB]));
//	}
//
//	LB = max(LB, max(LB_row, LB_col));
//	return LB;
//}

//bool LB_corner(vector<Point>* A, vector<Point>*B) {
//	int LengthA = (*A).size();
//	int LengthB = (*B).size();
//	//start
//	// i=1 j=0..lengthB
//	for (int j = 0; j < LengthB; j++) {
//		if (bool_dist(&((*A)[1]), &((*B)[j]))) break;
//		if ((!bool_dist(&((*A)[0]), &((*B)[j]))) || (j == (LengthB - 1))) return 1;
//	}
//	// i=0..lengthA j=1
//	for (int i = 0; i < LengthA; i++) {
//		if (bool_dist(&((*A)[i]), &((*B)[1]))) break;
//		if ((!bool_dist(&((*A)[i]), &((*B)[0]))) || (i == (LengthA - 1))) return 1;
//	}
//
//	//end
//	// i=LengthA-2 j=0..lengthB
//	for (int j = LengthB - 1; j >= 0; j++) {
//		if (bool_dist(&((*A)[LengthA - 2]), &((*B)[j]))) break;
//		if ((!bool_dist(&((*A)[LengthA - 1]), &((*B)[j]))) || (j == 0)) return 1;
//	}
//	// i=0..lengthA j=1
//	for (int i = LengthA - 1; i >= 0; i++) {
//		if (bool_dist(&((*A)[i]), &((*B)[LengthB - 2]))) break;
//		if ((!bool_dist(&((*A)[i]), &((*B)[LengthB - 1])) || (i == 0))) return 1;
//	}
//	return 0;
//}

//bool DFD(vector<Point>* A, vector<Point>* B) { // standard DFD with bool
//	int LengthA = (*A).size();
//	int LengthB = (*B).size();
//
//	for (int i = 0; i < LengthA; i++) {
//		for (int j = 0; j < LengthB; j++) {
//			bool flag = bool_dist(&((*A)[i]), &((*B)[j]));
//			if (!flag) { f[i % 2][j] = flag; continue; }
//			if ((i == 0) && (j == 0)) { f[i % 2][j] = flag; continue; }
//			if (i == 0) { f[i % 2][j] = f[i % 2][j - 1] & flag; continue; }
//			if (j == 0) { f[i % 2][j] = f[(i - 1) % 2][j] & flag; continue; }
//			f[i % 2][j] = (f[(i - 1) % 2][j] | f[i % 2][j - 1] | f[(i - 1) % 2][j - 1]) & flag;
//		}
//	}
//
//	return f[(LengthA - 1) % 2][LengthB - 1];
//}

double double_DFD(vector<Point>& A, vector<Point>& B) { //	standard DFD
	int LengthA = A.size();
	int LengthB = B.size();
	for (int i = 0; i < LengthA; i++) {
		for (int j = 0; j < LengthB; j++) {
			double distemp = double_dist(A[i], B[j]);
			if ((i == 0) && (j == 0)) { g[i % 2][j] = distemp; continue; }
			if (i == 0) { g[i % 2][j] = max(g[i % 2][j - 1], distemp); continue; }
			if (j == 0) { g[i % 2][j] = max(g[(i - 1) % 2][j], distemp); continue; }
			g[i % 2][j] = max(min(min(g[(i - 1) % 2][j], g[i % 2][j - 1]), g[(i - 1) % 2][j - 1]), distemp);
		}
	}
	return g[(LengthA - 1) % 2][LengthB - 1];
}

double DFD_LBrow(vector<Point>& A, vector<Point>& B) { //	standard DFD with LBrow
	int LengthA = A.size();
	int LengthB = B.size();
	for (int i = 0; i < LengthA; i++) {
		double LBrow = INFINITE;
		for (int j = 0; j < LengthB; j++) {
			double distemp = double_dist(A[i], B[j]);
			LBrow = min(distemp, LBrow);
			if ((i == 0) && (j == 0)) { g[i % 2][j] = distemp; continue; }
			if (i == 0) { g[i % 2][j] = max(g[i % 2][j - 1], distemp); continue; }
			if (j == 0) { g[i % 2][j] = max(g[(i - 1) % 2][j], distemp); continue; }
			g[i % 2][j] = max(min(min(g[(i - 1) % 2][j], g[i % 2][j - 1]), g[(i - 1) % 2][j - 1]), distemp);
		}
		if (LBrow > epsilon) return LBrow;
	}
	dfd_flag = true;
	return g[(LengthA - 1) % 2][LengthB - 1];
}

double DFD_greedy(vector<Point>& A, vector<Point>& B) { // greedy solution
	int LengthA = A.size() - 1;
	int LengthB = B.size() - 1;
	int i = 0, j = 0; double ans = double_dist(A[0], B[0]);
	while ((i < LengthA) || (j < LengthB)) {
		double dis1 = INFINITE, dis2 = INFINITE, dis3 = INFINITE;
		if (i < LengthA) dis1 = double_dist(A[i + 1], B[j]);
		if (j < LengthB) dis2 = double_dist(A[i], B[j + 1]);
		if ((i < LengthA) && (j < LengthB)) dis3 = double_dist(A[i + 1], B[j + 1]);

		if (dis1 < dis2) {
			if (dis1 < dis3) { i = i + 1; ans = max(ans, dis1); }
			else { i = i + 1; j = j + 1; ans = max(ans, dis3); }
		}
		else {
			if (dis2 < dis3) { j = j + 1; ans = max(ans, dis2); }
			else { i = i + 1; j = j + 1; ans = max(ans, dis3); }
		}
	}
	return ans;
}

double EP_DFD(vector<Point>& A, vector<Point>& B) { // early pruning 0-corner or 0 row;
	int LengthA = A.size();
	int LengthB = B.size();

	bool f_col[MaxLen];// f_col[i]=false if this column is all > epsilon;
	memset(f_col, false, sizeof(f_col));

	for (int i = 0; i < LengthA; i++) {
		double lbrow = INFINITE;
		for (int j = 0; j < LengthB; j++) {
			double distemp = double_dist(A[i], B[j]);
			lbrow = min(lbrow, distemp);
			if ((i == 0) && (j == 0)) { g[i % 2][j] = distemp; continue; }
			if (i == 0) { g[i % 2][j] = max(g[i % 2][j - 1], distemp); continue; }
			if (j == 0) { g[i % 2][j] = max(g[(i - 1) % 2][j], distemp); continue; }
			g[i % 2][j] = max(min(min(g[(i - 1) % 2][j], g[i % 2][j - 1]), g[(i - 1) % 2][j - 1]), distemp);
		}

		lbj = max(lbj, lbrow); //update lbj

		for (int j = 0; j < LengthB; j++) {
			if (g[i % 2][j]<epsilon) break;// find g[i][j]<epsilon can't prune
			if ((!f_col[j]) || (j == (LengthB - 1))) return INFINITE; // if exist 0-corner or 0-row
		}

		for (int j = 0; j < LengthB; j++) {
			if (g[i % 2][j] < epsilon) f_col[j] = true;
		}
	}

	dfd_flag = true;
	return g[(LengthA - 1) % 2][LengthB - 1];
}

bool EPplus_DFD(vector<Point>& A, vector<Point>& B) {
	int LengthA = A.size();
	int LengthB = B.size();
	
	memset(f_col, false, sizeof(f_col));
	memset(f_near, LengthB, sizeof(f_near));


	for (int i = 0; i < LengthA; i++) {
		int Nowi = i % 2; // reduce % times
		int Prei = (i - 1) % 2;
		int invalid_region = -1; // init all region is valid

		for (int j = 0; j < LengthB; j++) {
			if (j <= invalid_region) { // in invalid_region
				f[Nowi][j] = false;
				continue;
			}
			//cnt++;
			bool flag = bool_dist(A[i], B[j]);
			if (!flag) {
				f[Nowi][j] = flag;
				invalid_region = f_near[j] - 1; // update invalid_region;
				continue;
			}
			if ((i == 0) && (j == 0)) { f[Nowi][j] = flag; continue; }
			if (i == 0) { f[Nowi][j] = f[Nowi][j - 1] & flag; continue; }
			if (j == 0) { f[Nowi][j] = f[Prei][j] & flag; continue; }
			f[Nowi][j] = (f[Prei][j] | f[Nowi][j - 1] | f[Prei][j - 1]) & flag;
		}

		for (int j = 0; j < LengthB; j++) {
			if (f[Nowi][j]) break;// find a 1 can't prune
			if ((!f_col[j]) || (j == (LengthB - 1))) return false; // if exist 0-corner or 0-row
		}

		for (int j = 0; j < LengthB; j++) { // update f_col
			f_col[j] = f_col[j] | f[Nowi][j];
		}

		int nearest = LengthB;
		for (int j = LengthB - 1; j >= 0; j--) { //update f_near
			if (f[Nowi][j]) nearest = j;
			f_near[j] = nearest;
		}
	}

	return f[(LengthA - 1) % 2][LengthB - 1];
}

bool EPplusRQ_DFD(vector<Point>& A, vector<Point>& B,int A_num,int B_num) {// A_num and B_num is two trajectory in All_Query and All_Data
	int LengthA = A.size();
	int LengthB = B.size();

	memset(f_col, false, sizeof(f_col));
	memset(f_near, LengthB, sizeof(f_near));

	int PrePos = 0; // lbband of Qi < lbband of Qi+1

	for (int i = 0; i < LengthA; i++) {
		int Nowi = i % 2; // reduce % times
		int Prei = (i - 1) % 2;
		int invalid_region = -1; // init all region is valid

		//  use range query to replace distance computing
		for (int j = 0; j < LengthB; j++) DisRecord[j] = false;
		grouping::LB_bandFindAll(A[i], B_num, PrePos);

		for (int j = 0; j < LengthB; j++) {
			if (j <= invalid_region) { // in invalid_region
				f[Nowi][j] = false;
				continue;
			}
			//cnt++;
			bool flag;
			if (DisRecord[j]) flag = bool_dist(A[i], B[j]);
			else flag = false;
			//bool flag = DisRecord[j];

			if (!flag) {
				f[Nowi][j] = flag;
				invalid_region = f_near[j] - 1; // update invalid_region;
				continue;
			}
			if ((i == 0) && (j == 0)) { f[Nowi][j] = flag; continue; }
			if (i == 0) { f[Nowi][j] = f[Nowi][j - 1] & flag; continue; }
			if (j == 0) { f[Nowi][j] = f[Prei][j] & flag; continue; }
			f[Nowi][j] = (f[Prei][j] | f[Nowi][j - 1] | f[Prei][j - 1]) & flag;
		}

		for (int j = 0; j < LengthB; j++) {
			if (f[Nowi][j]) break;// find a 1 can't prune
			if ((!f_col[j]) || (j == (LengthB - 1))) return false; // if exist 0-corner or 0-row
		}

		for (int j = 0; j < LengthB; j++) { // update f_col
			f_col[j] = f_col[j] | f[Nowi][j];
		}

		int nearest = LengthB;
		for (int j = LengthB - 1; j >= 0; j--) { //update f_near
			if (f[Nowi][j]) nearest = j;
			f_near[j] = nearest;
		}
	}

	return f[(LengthA - 1) % 2][LengthB - 1];
}
//bool EPplus_RSeg_DFD(vector<Point>* A, vector<Point>* B, int datai, int dataj) {
//	int LengthA = (*A).size();
//	int LengthB = (*B).size();
//
//	memset(f_col, false, sizeof(f_col));
//	memset(f_near, LengthB, sizeof(f_near));
//
//
//	for (int i = 0; i < LengthA; i++) {
//		int Nowi = i % 2; // reduce % times
//		int Prei = (i - 1) % 2;
//		int invalid_region = -1; // init all region is valid
//
//		segLen = -1;//init segLen and value
//		bool flag;// bool_dist 
//		for (int j = 0; j < LengthB; j++) {
//
//			if (j <= invalid_region) { // in invalid_region
//				f[Nowi][j] = false;
//				continue;
//			}
//			//cnt++;
//			
//			//bool flag = bool_dist(&((*A)[i]), &((*B)[j]));
//			if (j > segLen) {
//				flag = bool_dist(&((*A)[i]), &((*B)[j]));
//				if (flag) { bool ff = in_RsegSearch(&All_Rseg[dataj], &(*A)[i], All_RsegRoot[dataj], 0, LengthB - 1, j); }
//				else { bool ff = out_RsegSearch(&All_Rseg[dataj], &(*A)[i], All_RsegRoot[dataj], 0, LengthB - 1, j); }
//			}
//
//			if (!flag) {
//				f[Nowi][j] = flag;
//				invalid_region = f_near[j] - 1; // update invalid_region;
//				continue;
//			}
//			if ((i == 0) && (j == 0)) { f[Nowi][j] = flag; continue; }
//			if (i == 0) { f[Nowi][j] = f[Nowi][j - 1] & flag; continue; }
//			if (j == 0) { f[Nowi][j] = f[Prei][j] & flag; continue; }
//			f[Nowi][j] = (f[Prei][j] | f[Nowi][j - 1] | f[Prei][j - 1]) & flag;
//		}
//
//		for (int j = 0; j < LengthB; j++) {
//			if (f[Nowi][j]) break;// find a 1 can't prune
//			if ((!f_col[j]) || (j == (LengthB - 1))) return false; // if exist 0-corner or 0-row
//		}
//
//		for (int j = 0; j < LengthB; j++) { // update f_col
//			f_col[j] = f_col[j] | f[Nowi][j];
//		}
//
//		int nearest = LengthB;
//		for (int j = LengthB - 1; j >= 0; j--) { //update f_near
//			if (f[Nowi][j]) nearest = j;
//			f_near[j] = nearest;
//		}
//	}
//
//	return f[(LengthA - 1) % 2][LengthB - 1];
//}

//bool EPplus_CSeg_DFD(vector<Point>* A, vector<Point>* B, int datai, int dataj) {
//	int LengthA = (*A).size();
//	int LengthB = (*B).size();
//
//	memset(f_col, false, sizeof(f_col));
//	memset(f_near, LengthB, sizeof(f_near));
//
//
//	for (int i = 0; i < LengthA; i++) {
//		int Nowi = i % 2; // reduce % times
//		int Prei = (i - 1) % 2;
//		int invalid_region = -1; // init all region is valid
//
//		segLen = -1;//init segLen and value
//		bool flag;// bool_dist 
//		for (int j = 0; j < LengthB; j++) {
//
//			if (j <= invalid_region) { // in invalid_region
//				f[Nowi][j] = false;
//				continue;
//			}
//			//cnt++;
//
//			//bool flag = bool_dist(&((*A)[i]), &((*B)[j]));
//			if (j > segLen) {
//				flag = bool_dist(&((*A)[i]), &((*B)[j]));
//				if (flag) { bool ff = in_CsegSearch(&All_Cseg[dataj], &(*A)[i], All_CsegRoot[dataj], 0, LengthB - 1, j); }
//				else { bool ff = out_CsegSearch(&All_Cseg[dataj], &(*A)[i], All_CsegRoot[dataj], 0, LengthB - 1, j); }
//			}
//
//			if (!flag) {
//				f[Nowi][j] = flag;
//				invalid_region = f_near[j] - 1; // update invalid_region;
//				continue;
//			}
//			if ((i == 0) && (j == 0)) { f[Nowi][j] = flag; continue; }
//			if (i == 0) { f[Nowi][j] = f[Nowi][j - 1] & flag; continue; }
//			if (j == 0) { f[Nowi][j] = f[Prei][j] & flag; continue; }
//			f[Nowi][j] = (f[Prei][j] | f[Nowi][j - 1] | f[Prei][j - 1]) & flag;
//		}
//
//		for (int j = 0; j < LengthB; j++) {
//			if (f[Nowi][j]) break;// find a 1 can't prune
//			if ((!f_col[j]) || (j == (LengthB - 1))) return false; // if exist 0-corner or 0-row
//		}
//
//		for (int j = 0; j < LengthB; j++) { // update f_col
//			f_col[j] = f_col[j] | f[Nowi][j];
//		}
//
//		int nearest = LengthB;
//		for (int j = LengthB - 1; j >= 0; j--) { //update f_near
//			if (f[Nowi][j]) nearest = j;
//			f_near[j] = nearest;
//		}
//	}
//
//	return f[(LengthA - 1) % 2][LengthB - 1];
//}
//void DFS_DFD(vector<Point>* A, vector<Point>* B,int i,int j) {
//	if (suc) return;
//	int LengthA = (*A).size();
//	int LengthB = (*B).size();
//	if ((i >= LengthA) || (j >= LengthB)) return;
//	pair<int, int> temp_pair;
//	temp_pair = make_pair(i, j);
//	if (rep.find(temp_pair)!=rep.end()) return;
//	rep.insert(temp_pair);
//
//	if ((i == (LengthA-1)) && (j == (LengthB-1))) {
//		suc = true;
//		return;
//	}
//
//	if (!bool_dist((*A)[i], (*B)[j])) return;
//
//	DFS_DFD(A, B, i + 1, j + 1);
//	DFS_DFD(A, B, i, j + 1);
//	DFS_DFD(A, B, i + 1, j);
//	return;
//}

//bool DFS_DFD(vector<Point>* A, vector<Point>* B) {
//	int LengthA = (*A).size();
//	int LengthB = (*B).size();
//	
//	rep.clear(); Vstack.clear();
//	pair<int, int> temppair = make_pair(0, 0);
//	Vstack.push_back(temppair);
//	rep.insert(temppair);
//	while (Vstack.size()!=0) {
//		temppair = Vstack.back();
//		Vstack.pop_back();
//		int i = temppair.first; int j = temppair.second;
//		if ((i >= LengthA) || (j >= LengthB)) continue;
//		cnt++;
//		if (!bool_dist((*A)[i], (*B)[j])) continue;
//		if ((i == (LengthA - 1)) && (j == (LengthB - 1))) return true;
//
//		temppair = make_pair(i + 1, j + 1);
//		if (rep.find(temppair) == rep.end()) {
//			rep.insert(temppair);
//			Vstack.push_back(temppair);
//		}
//		temppair = make_pair(i + 1, j);
//		if (rep.find(temppair) == rep.end()) {
//			rep.insert(temppair);
//			Vstack.push_back(temppair);
//		}
//		temppair = make_pair(i, j + 1);
//		if (rep.find(temppair) == rep.end()) {
//			rep.insert(temppair);
//			Vstack.push_back(temppair);
//		}
//	}
//
//	return false;
//}

//bool DFS_DFD_Matrix(vector<Point>* A, vector<Point>* B) {
//	int LengthA = (*A).size();
//	int LengthB = (*B).size();
//	
//	Vstack.clear();
//	pair<int, int> temppair = make_pair(0, 0);
//	Vstack.push_back(temppair);
//	RepMatrix[0][0] = true;
//	tot = 0;// searched points
//	Searchedi[tot] = 0; Searchedj[tot] = 0;
//	while (Vstack.size() != 0) {
//		temppair = Vstack.back();
//		Vstack.pop_back();
//		int i = temppair.first; int j = temppair.second;
//		if ((i >= LengthA) || (j >= LengthB)) continue;
//		if (!bool_dist(&((*A)[i]), &((*B)[j]))) continue;
//		if ((i == (LengthA - 1)) && (j == (LengthB - 1))) return true;
//
//		temppair = make_pair(i + 1, j + 1);
//		if (!RepMatrix[i + 1][j + 1]) {
//			RepMatrix[i + 1][j + 1] = true; tot++; Searchedi[tot] = i + 1; Searchedj[tot] = j + 1;
//			Vstack.push_back(temppair);
//		}
//		temppair = make_pair(i, j + 1);
//		if (!RepMatrix[i][j + 1]) {
//			RepMatrix[i][j + 1] = true; tot++; Searchedi[tot] = i; Searchedj[tot] = j + 1;
//			Vstack.push_back(temppair);
//		}
//		temppair = make_pair(i + 1, j);
//		if (!RepMatrix[i + 1][j]) {
//			RepMatrix[i + 1][j] = true; tot++; Searchedi[tot] = i + 1; Searchedj[tot] = j;
//			Vstack.push_back(temppair);
//		}
//	}
//
//	return false;
//}

//bool EPplus_DFD(vector<Point>* A, vector<Point>* B) { // space N*N
//	int LengthA = (*A).size();
//	int LengthB = (*B).size();
//
//	bool f_col[MaxLen];	// f_col[i]=false if this column is all 0;
//	memset(f_col, false, sizeof(f_col));
//
//	int f_near[MaxLen]; // the nearest 1 on the right;
//	memset(f_near, LengthB, sizeof(f_near));
//
//	
//
//	for (int i = 0; i < LengthA; i++) {
//
//		int invalid_region = -1; // init all region is valid
//
//		for (int j = 0; j < LengthB; j++) {
//			if (j <= invalid_region) { // in invalid_region
//				f[i][j] = false;
//				continue;
//			}
//			bool flag = bool_dist((*A)[i], (*B)[j]);
//			if (!flag) { 
//				f[i][j] = flag; 
//				invalid_region = f_near[j] - 1; // update invalid_region;
//				continue; 
//			}
//			if ((i == 0) && (j == 0)) { f[i][j] = flag; continue; }
//			if (i == 0) { f[i][j] = f[i][j - 1] & flag; continue; }
//			if (j == 0) { f[i][j] = f[i - 1][j] & flag; continue; }
//			f[i][j] = (f[i - 1][j] | f[i][j - 1] | f[i - 1][j - 1]) & flag;
//		}
//
//		for (int j = 0; j < LengthB; j++) {
//			if (f[i][j]) break;// find a 1 can't prune
//			if ((!f_col[j]) || (j == (LengthB - 1))) return false; // if exist 0-corner or 0-row
//		}
//
//		for (int j = 0; j < LengthB; j++) { // update f_col
//			f_col[j] = f_col[j] | f[i][j];
//		}
//
//		int nearest = LengthB;
//		for (int j = LengthB - 1; j >= 0; j--) { //update f_near
//			if (f[i][j]) nearest = j;
//			f_near[j] = nearest;
//		}
//	}
//
//	return f[LengthA - 1][LengthB - 1];
//}
#endif