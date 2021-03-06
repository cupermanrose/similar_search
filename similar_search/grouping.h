#ifndef _GROUPING_H
#define _GROUPING_H

#include <vector>
#include <iostream> 
#include <algorithm>
#include <BLtoDist.h>
#include <Init.h>

using namespace std;

namespace grouping {

	const int tau = 32;
	const int dimension = 2;
	const double eps = 1.0E-10;

	double g[2][tau];

	struct MBR {
		double U[dimension], L[dimension];
	};

	struct GroupTra {
		int length,number;
		int position[tau]; // store ending position of each position
		MBR MBR[tau];
	};

	vector<GroupTra> AllGTra;

	void ClearMBR(MBR& MBR) {
		for (int i = 0; i < dimension; i++) {
			MBR.L[i] = INFINITY;
			MBR.U[i] = -INFINITY;
		}
		return;
	}

	void CreateMBR(MBR& MBR, Trajectory& T, int left, int right) { // create MBR from left to right-1
		ClearMBR(MBR);
		for (int i = 0; i < dimension; i++) {
			for (int j = left; j < right; j++) {
				if (i == 0) {
					MBR.L[0] = min(MBR.L[0], T.Points[j].latitude);
					MBR.U[0] = max(MBR.U[0], T.Points[j].latitude);
				}
				else {
					MBR.L[1] = min(MBR.L[1], T.Points[j].longitude);
					MBR.U[1] = max(MBR.U[1], T.Points[j].longitude);
				}
			}
		}
		return;
	}

	double UDis(MBR& M1, MBR& M2) {
		double res = 0;
		for (int i = 0; i < dimension; i++) {
			double dis1 = abs(M1.U[i] - M2.L[i]);
			double dis2 = abs(M1.L[i] - M2.U[i]);
			dis1 = max(dis1, dis2);
			res = res + dis1*dis1;
		}
		//return sqrt(res);
		return res;
	}

	double LDis(MBR& M1, MBR& M2) {
		double res = 0;
		for (int i = 0; i < dimension; i++) {
			if ((M1.U[i] + eps > M2.U[i]) && (M1.L[i] - eps < M2.U[i])) continue;
			if ((M2.U[i] + eps > M1.U[i]) && (M2.L[i] - eps < M1.U[i])) continue;
			/*if ((M1.U[i] + eps > M2.L[i]) && (M1.U[i] - eps < M2.U[i])) continue;
			if ((M1.L[i] + eps > M2.L[i]) && (M1.L[i] - eps < M2.U[i])) continue;
			if ((M2.U[i] + eps > M1.L[i]) && (M2.U[i] - eps < M1.U[i])) continue;
			if ((M2.L[i] + eps > M1.L[i]) && (M2.L[i] - eps < M1.U[i])) continue;*/
			double dis1 = abs(M1.U[i] - M2.L[i]);
			double dis2 = abs(M1.L[i] - M2.U[i]);
			dis1 = min(dis1, dis2);
			res = res + dis1*dis1;
		}
		//return sqrt(res);
		return res;
	}
	
	double GlbDFD(GroupTra& A, GroupTra& B) { //	standard DFD lb
		int LengthA = A.length;
		int LengthB = B.length;
		for (int i = 0; i < LengthA; i++) {
			for (int j = 0; j < LengthB; j++) {
				double distemp = LDis(A.MBR[i], B.MBR[j]);
				if ((i == 0) && (j == 0)) { g[i % 2][j] = distemp; continue; }
				if (i == 0) { g[i % 2][j] = max(g[i % 2][j - 1], distemp); continue; }
				if (j == 0) { g[i % 2][j] = max(g[(i - 1) % 2][j], distemp); continue; }
				g[i % 2][j] = max(min(min(g[(i - 1) % 2][j], g[i % 2][j - 1]), g[(i - 1) % 2][j - 1]), distemp);
			}
		}
		return g[(LengthA - 1) % 2][LengthB - 1];
	}

	double GubDFD(GroupTra& A, GroupTra& B) { //	standard DFD ub
		int LengthA = A.length;
		int LengthB = B.length;
		for (int i = 0; i < LengthA; i++) {
			for (int j = 0; j < LengthB; j++) {
				double distemp = UDis(A.MBR[i], B.MBR[j]);
				if ((i == 0) && (j == 0)) { g[i % 2][j] = distemp; continue; }
				if (i == 0) { g[i % 2][j] = max(g[i % 2][j - 1], distemp); continue; }
				if (j == 0) { g[i % 2][j] = max(g[(i - 1) % 2][j], distemp); continue; }
				g[i % 2][j] = max(min(min(g[(i - 1) % 2][j], g[i % 2][j - 1]), g[(i - 1) % 2][j - 1]), distemp);
			}
		}
		return g[(LengthA - 1) % 2][LengthB - 1];
	}

	double GlbDFD_LBrow(GroupTra& A, GroupTra& B) {
		int LengthA = A.length;
		int LengthB = B.length;
		for (int i = 0; i < LengthA; i++) {
			double LBrow = INFINITE;
			for (int j = 0; j < LengthB; j++) {
				double distemp = LDis(A.MBR[i], B.MBR[j]);
				LBrow = min(distemp, LBrow);
				if ((i == 0) && (j == 0)) { g[i % 2][j] = distemp; continue; }
				if (i == 0) { g[i % 2][j] = max(g[i % 2][j - 1], distemp); continue; }
				if (j == 0) { g[i % 2][j] = max(g[(i - 1) % 2][j], distemp); continue; }
				g[i % 2][j] = max(min(min(g[(i - 1) % 2][j], g[i % 2][j - 1]), g[(i - 1) % 2][j - 1]), distemp);
			}
			if (LBrow > epsilon) return LBrow;
		}
		return g[(LengthA - 1) % 2][LengthB - 1];
	}

	double GubDFD_LBrow(GroupTra& A, GroupTra& B) {
		int LengthA = A.length;
		int LengthB = B.length;
		for (int i = 0; i < LengthA; i++) {
			double LBrow = INFINITE;
			for (int j = 0; j < LengthB; j++) {
				double distemp = UDis(A.MBR[i], B.MBR[j]);
				LBrow = min(distemp, LBrow);
				if ((i == 0) && (j == 0)) { g[i % 2][j] = distemp; continue; }
				if (i == 0) { g[i % 2][j] = max(g[i % 2][j - 1], distemp); continue; }
				if (j == 0) { g[i % 2][j] = max(g[(i - 1) % 2][j], distemp); continue; }
				g[i % 2][j] = max(min(min(g[(i - 1) % 2][j], g[i % 2][j - 1]), g[(i - 1) % 2][j - 1]), distemp);
			}
			if (LBrow > epsilon) return LBrow;
		}
		return g[(LengthA - 1) % 2][LengthB - 1];
	}

	double GetDistPM(Point& A, MBR& M) { //  min dist
		double res = 0;
		for (int i = 0; i < dimension; i++) {
			if (i == 0) {
				if (A.latitude < M.L[0]) res = res + (M.L[0] - A.latitude)*(M.L[0] - A.latitude);
				if (A.latitude > M.U[0]) res = res + (A.latitude - M.U[0])*(A.latitude - M.U[0]);
			}
			else {
				if (A.longitude < M.L[1]) res = res + (M.L[1] - A.longitude)*(M.L[1] - A.longitude);
				if (A.longitude > M.U[1]) res = res + (A.longitude - M.U[1])*(A.longitude - M.U[1]);
			}
		}
		//return sqrt(res);
		return res;
	}

	double GetPdist(Point&A, double lat, double lon) {
		Point B;
		B.latitude = lat; B.longitude = lon;
		return double_dist(A, B);
	}

	bool LB_band(int Q, int A) {// true: has a LB_band; false no lower bound
		vector<Point>& TempQ = All_Query[Q].Points;
		GroupTra& TempGT = AllGTra[A];
		double d1, d2, d3, d4, d;
		int PrePos = 0; // lbband of Q_i must < lbband of Q_i+1
		for (int i = 0; i < All_Query[Q].Points.size(); i++) {
			bool flag = false;
			for (int j = 0; j < TempGT.length; j++) {
				if (flag) continue;
				d = GetDistPM(TempQ[i], TempGT.MBR[j]);
				d1 = GetPdist(TempQ[i], TempGT.MBR[j].L[0], TempGT.MBR[j].L[1]);
				d2 = GetPdist(TempQ[i], TempGT.MBR[j].L[0], TempGT.MBR[j].U[1]);
				d3 = GetPdist(TempQ[i], TempGT.MBR[j].U[0], TempGT.MBR[j].L[1]);
				d4 = GetPdist(TempQ[i], TempGT.MBR[j].U[0], TempGT.MBR[j].U[1]);
				if ((d1 < epsilon) && (d2 < epsilon)) { flag = true; continue; }
				if ((d1 < epsilon) && (d3 < epsilon)) { flag = true; continue; }
				if ((d2 < epsilon) && (d4 < epsilon)) { flag = true; continue; }
				if ((d3 < epsilon) && (d4 < epsilon)) { flag = true; continue; }
				if (d < epsilon) {
					int st, en; // position of Tra
					if (j == 0) { st = 0; }
					else { st = AllGTra[A].position[j - 1]; }
					en = AllGTra[A].position[j];
					for (int k = st; k < en; k++) {
						if (k < PrePos) continue;
						if (double_dist(TempQ[i], All_Data[A].Points[k]) < epsilon) flag = true;
						if (flag) PrePos = k;
						if (flag) break;
					}
				}
			}
			if (!flag) return true; // can find any point is close to All_Query[Q].Points[i]
		}
		return false;
	}

	void LB_bandFindAll(Point& Qi, int A, int& PrePos) {
		GroupTra& TempGT = AllGTra[A];
		for (int i = 0; i < TempGT.length; i++) {
			double d = GetDistPM(Qi, TempGT.MBR[i]);
			int st, en; // position of Tra
			if (i == 0) { st = 0; }
			else { st = AllGTra[A].position[i - 1]; }
			en = AllGTra[A].position[i];

			if (d < epsilon) {
				//int st, en; // position of Tra
				//if (i == 0) { st = 0; }
				//else { st = AllGTra[A].position[i - 1]; }
				//en = AllGTra[A].position[i];

				//for (int k = st; k < en; k++) {
				//	if (k < PrePos) continue;
				//	DisRecord[k] = bool_dist(Qi, All_Data[A].Points[k]);
				//}
				for (int k = st; k < en; k++) DisRecord[k] = true;
			}
			else for (int k = st; k < en; k++) DisRecord[k] = false;
		}

		/*for (int i = 0; i < All_Data[A].Points.size(); i++) {
			if (DisRecord[i]) {
				PrePos = i;
				return;
			}
		}*/
		
		return;
	}

	void DivideTrajectory(Trajectory& T, int TID, GroupTra& NewT) {
		NewT.number = TID;
		if (T.Points.size() < tau) {
			NewT.length = 1;
			NewT.position[0] = T.Points.size();
			CreateMBR(NewT.MBR[0], T, 0, NewT.position[0]);
			return;
		}

		NewT.length = tau;
		int unit = T.Points.size() / tau;
		for (int i = 0; i < tau; i++) {
			int left, right;
			left = i*unit;
			if (i == tau - 1) right = T.Points.size();
			else right = left + unit;
			NewT.position[i] = right;
			//cout << left << ' ' << right << ' ' << T.Points.size() << endl;
			CreateMBR(NewT.MBR[i], T, left, right);
		}
		return;
	}

	void Init() {
		AllGTra.clear();
		for (int i = 0; i < All_Data.size(); i++) {
			GroupTra NewT;
			DivideTrajectory(All_Data[i], i,  NewT);
			AllGTra.push_back(NewT);
		}
	}
}

#endif