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

	void CreateMBR(MBR& MBR, int TID, int left, int right) { // create MBR from left to right-1
		ClearMBR(MBR);
		for (int i = 0; i < dimension; i++) {
			for (int j = left; j < right; j++) {
				if (i == 0) {
					MBR.L[0] = min(MBR.L[0], All_Data[TID].Points[j].latitude);
					MBR.U[0] = max(MBR.U[0], All_Data[TID].Points[j].latitude);
				}
				else {
					MBR.L[1] = min(MBR.L[1], All_Data[TID].Points[j].longitude);
					MBR.U[1] = max(MBR.U[1], All_Data[TID].Points[j].longitude);
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
		return sqrt(res);
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
		return sqrt(res);
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

	void DivideTrajectory(Trajectory& T, int TID, GroupTra& NewT) {
		NewT.number = TID;
		if (T.Points.size() < tau) {
			NewT.length = 1;
			NewT.position[0] = T.Points.size();
			CreateMBR(NewT.MBR[0], NewT.number, 0, NewT.position[0]);
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
			CreateMBR(NewT.MBR[i], NewT.number, left, right);
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