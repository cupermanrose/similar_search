#ifndef _SIMILARITY_H
#define _SIMILARITY_H

#include <algorithm>
#include <iterator>
#include <iomanip>
#include <M_tree.h>

using namespace std;

set<int> Start_set, End_set, Intersect_set;
double lower[20000], upper[20000];
vector<int> answer;

void similarity_search() {
	int anssum = 0;
	for (int i = 0; i < All_Query.size(); i++) {
		answer.clear();
		for (int j = 0; j < All_Data.size(); j++) {
			double bdis = double_DFD(&All_Query[i].Points, &All_Data[j].Points);
			if (bdis < epsilon) answer.push_back(j);
		}
		anssum = anssum + answer.size();
	}
	fout << "Brute force answer: " << anssum << endl;
}

void similarity_search_mtree() {
	int anssum = 0;
	for (int i = 0; i < All_Query.size(); i++) {
	//for (int i = 0; i < 30; i++) {
		answer.clear();
		results.clear();
		RangeQuery_mtree(root_mtree, i, epsilon);
		//cout << i << ": ";
		/*for (int j = 0; j < results.size(); j++) {
			cout << results[j] << " " ;
		}
		cout << endl;*/
		anssum = anssum + results.size();
	}
	fout << "Mtree answer: " << anssum << endl;
}

void similarity_search_baseline() {
	int anssum = 0, lbcell = 0;
	for (int i = 0; i < All_Query.size(); i++) {
		answer.clear();
		for (int j = 0; j < All_Data.size(); j++) {
			if (LB_cell(&All_Query[i].Points, &All_Data[j].Points) > epsilon) continue;
			lbcell++;
			double bdis = DFD_LBrow(&All_Query[i].Points, &All_Data[j].Points);
			if (bdis < epsilon) answer.push_back(j);
		}
		anssum = anssum + answer.size();
	}
	fout << "lbcell: " << lbcell << endl;
	fout << "baseline answer: " << anssum << endl;
}

void update_lbub(int x, double lb, double ub, double exact_dfd) {
	for (int p = 0; p < ExactDFD[x].size(); p++) {
		int y = ExactDFD[x][p].first;
		if (y >= All_Data.size()) continue;
		lower[y] = max(lower[y], lb - ExactDFD[x][p].second);
		lower[y] = max(lower[y], ExactDFD[x][p].second - ub);
		if (dfd_flag) lower[y] = max(lower[y], abs(exact_dfd - ExactDFD[x][p].second));
		if (dfd_flag) upper[y] = min(upper[y], ExactDFD[x][p].second + exact_dfd);
		upper[y] = min(upper[y], ExactDFD[x][p].second + ub);
	}
	return;
}

void similarity_search_triangle() {

	int anssum = 0;
	int lbcell = 0, trilb = 0, triub = 0, lbband = 0, ubgreedy = 0, finaldfd = 0;

	for (int i = 0; i < All_Query.size(); i++) {

		answer.clear();
		for (int j = 0; j < All_Data.size(); j++) {
			lower[j] = 0;upper[j] = INFINITE;
		}
		
		//LBcell
		Start_set.clear(); End_set.clear(); Intersect_set.clear();
		Range_KDsearch(&Start_KD, &Start_set, All_Query[i].Points[0].latitude, All_Query[i].Points[0].longitude, root_Start);
		Range_KDsearch(&End_KD, &End_set, All_Query[i].Points.back().latitude, All_Query[i].Points.back().longitude, root_End);	
		set_intersection(Start_set.begin(), Start_set.end(), End_set.begin(), End_set.end(), inserter(Intersect_set, Intersect_set.end()));
		
		set<int>::iterator it;
		for (it = Intersect_set.begin(); it != Intersect_set.end(); it++) {
			int j = *it;

			lbj = lower[j];
			ubj = upper[j];
			double dfdj;
			dfd_flag = false;

			lbcell++;

			if (lower[j] > epsilon) { update_lbub(j, lbj, ubj, dfdj); continue; }
			trilb++;
			if (upper[j] < epsilon) { update_lbub(j, lbj, ubj, dfdj); answer.push_back(j); continue; }
			triub++;

			//LB_band
			bool flag = false;
			for (int p = 0; p < All_Query[i].Points.size(); p++) {
				if (!Range_KDsearch_forband(&All_KD[j], All_Query[i].Points[p].latitude, All_Query[i].Points[p].longitude, All_KDroot[j])) { flag = true; break; }
			}
			if (flag) { update_lbub(j, lbj, ubj, dfdj); continue; }
			flag = false;
			for (int p = 0; p < All_Data[j].Points.size(); p++) {
				if (!Range_KDsearch_forband(&All_KD[i], All_Data[j].Points[p].latitude, All_Data[j].Points[p].longitude, All_KDroot[i])) { flag = true; break; }
			}
			if (flag) { update_lbub(j, lbj, ubj, dfdj); continue; }
			lbband++;

			//UB_greedy
			double ub_temp = DFD_greedy(&All_Query[i].Points, &All_Data[j].Points);
			ubj = min(ubj, ub_temp);
			if (ub_temp < epsilon) { update_lbub(j, lbj, ubj, dfdj); answer.push_back(j); }
			else {
				ubgreedy++;
				dfdj= EP_DFD(&All_Query[i].Points, &All_Data[j].Points);
				if (dfdj < epsilon) { finaldfd++; answer.push_back(j); }
				update_lbub(j, lbj, ubj, dfdj);
			}
		}

		anssum = anssum + answer.size(); // the similarity amount of all query
	}
	fout << "lbcell: " << lbcell << endl;
	fout << "trilb: " << trilb << endl;
	fout << "triub: " << triub << endl;
	fout << "lbband: " << lbband << endl;
	fout << "ubgreedy: " << ubgreedy << endl;
	fout << "finaldfd: " << finaldfd << endl;
	fout << "triangle answer: " << anssum << endl;
}

//void similarity_search_Index_EPplus() {
//	int anssum = 0;
//	int lbcell = 0, lbcorner = 0, lbband = 0;
//	memset(RepMatrix, false, sizeof(RepMatrix));
//	set<int> Start_set, End_set, Intersect_set, All_set, res_set;
//	for (int i = 0; i < All_Query.size(); i++) {
//		answer.clear();
//		//LB_cell
//		Start_set.clear(); End_set.clear(); Intersect_set.clear();
//		Range_KDsearch(&Start_KD, &Start_set, All_Query[i].Points[0].latitude, All_Query[i].Points[0].longitude, root_Start);
//		Range_KDsearch(&End_KD, &End_set, All_Query[i].Points.back().latitude, All_Query[i].Points.back().longitude, root_End);
//		set_intersection(Start_set.begin(), Start_set.end(), End_set.begin(), End_set.end(), inserter(Intersect_set, Intersect_set.end()));
//
//		set<int>::iterator it;
//		for (it = Intersect_set.begin(); it != Intersect_set.end(); it++) {
//			int j = *it;
//			lbcell++;
//			//if (LB_cell(&All_Query[i].Points, &All_Data[j].Points) > epsilon) continue;
//			if (LB_corner(&All_Query[i].Points, &All_Data[j].Points)) continue;
//			
//			lbcorner++;
//			//LB_band
//			bool flag = false;
//			//for (int p = 0; p < All_Query[i].Points.size(); p = p + (All_Query[i].Points.size() / 50 + 1 )) {
//			for (int p = 0; p < All_Query[i].Points.size(); p++) {
//				All_set.clear();
//				if (!Range_KDsearch_forband(&All_KD[j], All_Query[i].Points[p].latitude, All_Query[i].Points[p].longitude, All_KDroot[j])) { flag = true; break; }
//			}
//			if (flag) continue;
//			flag = false;
//			for (int p = 0; p < All_Data[j].Points.size(); p++) {
//				All_set.clear();
//				if (!Range_KDsearch_forband(&All_KD[i], All_Data[j].Points[p].latitude, All_Data[j].Points[p].longitude, All_KDroot[i])) { flag = true; break; }
//			}
//			if (flag) continue;
//			lbband++;
//			bool bdis = EPplus_DFD(&All_Query[i].Points, &All_Data[j].Points);
//			if (bdis) answer.push_back(j);
//		}
//		anssum = anssum + answer.size(); // the similarity amount of all query
//	}
//	//fout << "LB_cell: " << lbcell << endl;
//	//fout << "LB_corner: " << lbcorner << endl;
//	//fout << "LB_band: " << lbband << endl;
//	fout << "Index answer: " << anssum << endl;
//}

//void similarity_search_LB_cellAndcross() {
//	int anssum = 0;
//	for (int i = 0; i < All_Query.size(); i++) {
//		answer.clear();
//		for (int j = 0; j < All_Data.size(); j++) {
//			if (LB_cell(&All_Query[i].Points, &All_Data[j].Points) > epsilon) continue;
//			if (LB_cross(&All_Query[i].Points, &All_Data[j].Points) > epsilon) continue;
//			bool bdis = DFD(&All_Query[i].Points, &All_Data[j].Points);
//			if (bdis) answer.push_back(j);
//		}
//		anssum = anssum + answer.size(); // the similarity amount of all query
//	}
//	fout << "baseline answer: " << anssum << endl;
//}

//void similarity_search_LB_cellAndcorner() {
//	int anssum = 0;
//	for (int i = 0; i < All_Query.size(); i++) {
//		answer.clear();
//		for (int j = 0; j < All_Data.size(); j++) {
//			if (LB_cell(&All_Query[i].Points, &All_Data[j].Points) > epsilon) continue;
//			if (LB_corner(&All_Query[i].Points, &All_Data[j].Points)) continue;
//			bool bdis = DFD(&All_Query[i].Points, &All_Data[j].Points);
//			if (bdis) answer.push_back(j);
//		}
//		anssum = anssum + answer.size(); 
//	}
//	fout << "LB_cell&corner: " << anssum << endl;
//}

//void similarity_search_EP() {
//	int anssum = 0;
//	for (int i = 0; i < All_Query.size(); i++) {
//		answer.clear();
//		for (int j = 0; j < All_Data.size(); j++) {
//			if (LB_cell(&All_Query[i].Points, &All_Data[j].Points) > epsilon) continue;
//			if (LB_corner(&All_Query[i].Points, &All_Data[j].Points)) continue;
//			bool bdis = EP_DFD(&All_Query[i].Points, &All_Data[j].Points);
//			if (bdis) answer.push_back(j);
//		}
//		anssum = anssum + answer.size(); // the similarity amount of all query
//	}
//	fout << "EP_DFD answer: " << anssum << endl;
//}
//
//void similarity_search_EPplus() {
//	int anssum = 0; cnt = 0;
//	for (int i = 0; i < All_Query.size(); i++) {
//		answer.clear();
//		for (int j = 0; j < All_Data.size(); j++) {
//			if (LB_cell(&All_Query[i].Points, &All_Data[j].Points) > epsilon) continue;
//			if (LB_corner(&All_Query[i].Points, &All_Data[j].Points)) continue;
//			bool bdis = EPplus_DFD(&All_Query[i].Points, &All_Data[j].Points);
//			if (bdis) answer.push_back(j);
//		}
//		anssum = anssum + answer.size(); // the similarity amount of all query
//	}
//	fout << "EPplus_DFD distance calculate: " << cnt << endl;
//	fout << "EPplus_DFD answer: " << anssum << endl;
//}

//void similarity_search_DFS_DFD() {
//	int anssum = 0; cnt = 0;
//	for (int i = 0; i < All_Query.size(); i++) {
//		answer.clear();
//		for (int j = 0; j < All_Data.size(); j++) {
//			if (LB_cell(&All_Query[i].Points, &All_Data[j].Points) > epsilon) continue;
//			//if (LB_corner(&All_Query[i].Points, &All_Data[j].Points)) continue;
//			//rep.clear();
//			//suc = false;
//			//DFS_DFD(&All_Query[i].Points, &All_Data[j].Points,0,0); // search from (0,0)
//			//if (suc) answer.push_back(j);
//			if (DFS_DFD(&All_Query[i].Points, &All_Data[j].Points)) answer.push_back(j);
//		}
//		anssum = anssum + answer.size(); // the similarity amount of all query
//	}
//	fout << "DFS_DFD distance calculate: " << cnt << endl;
//	fout << "DFS_DFD answer: " << anssum << endl;
//}

//void similarity_search_DFS_DFD_Matrix() {
//	int anssum = 0; cnt = 0;
//	memset(RepMatrix, false, sizeof(RepMatrix));
//	for (int i = 0; i < All_Query.size(); i++) {
//		answer.clear();
//		for (int j = 0; j < All_Data.size(); j++) {
//			if (LB_cell(&All_Query[i].Points, &All_Data[j].Points) > epsilon) continue;
//			if (LB_corner(&All_Query[i].Points, &All_Data[j].Points)) continue;
//			if (DFS_DFD_Matrix(&All_Query[i].Points, &All_Data[j].Points)) answer.push_back(j);
//			for (int k = 0; k <= tot; k++) RepMatrix[Searchedi[k]][Searchedj[k]] = false;
//		}
//		anssum = anssum + answer.size(); // the similarity amount of all query
//	}
//	fout << "DFS_DFD_Matrix distance calculate: " << cnt << endl;
//	fout << "DFS_DFD_Matrix answer: " << anssum << endl;
//}

//void similarity_search_Combine() {
//	int anssum = 0; cnt = 0;
//	memset(RepMatrix, false, sizeof(RepMatrix));
//	for (int i = 0; i < All_Query.size(); i++) {
//		answer.clear();
//		for (int j = 0; j < All_Data.size(); j++) {
//			if (LB_cell(&All_Query[i].Points, &All_Data[j].Points) > epsilon) continue;
//			if (LB_corner(&All_Query[i].Points, &All_Data[j].Points)) continue;
//			
//			if ((All_Query[i].Points.size() < MatrixLen) && (All_Query[j].Points.size() < MatrixLen)) {
//				if (DFS_DFD_Matrix(&All_Query[i].Points, &All_Data[j].Points)) answer.push_back(j);
//				for (int k = 0; k <= tot; k++) RepMatrix[Searchedi[k]][Searchedj[k]] = false;
//			}
//			else {
//				bool bdis = EPplus_DFD(&All_Query[i].Points, &All_Data[j].Points);
//				if (bdis) answer.push_back(j);
//			}
//		}
//		anssum = anssum + answer.size(); // the similarity amount of all query
//	}
//	fout << "Combine calculate: " << cnt << endl;
//	fout << "Combine answer: " << anssum << endl;
//}

//void similarity_search_Index() {
//	int anssum = 0;
//	int lbcell = 0, lbcorner = 0, lbband = 0;
//	memset(RepMatrix, false, sizeof(RepMatrix));
//	set<int> Start_set, End_set, Intersect_set, All_set, res_set;
//	for (int i = 0; i < All_Query.size(); i++) {
//		answer.clear();
//		//LB_cell
//		Start_set.clear(); End_set.clear(); Intersect_set.clear();
//		Range_KDsearch(&Start_KD, &Start_set, All_Query[i].Points[0].latitude, All_Query[i].Points[0].longitude, root_Start);
//		Range_KDsearch(&End_KD, &End_set, All_Query[i].Points.back().latitude, All_Query[i].Points.back().longitude, root_End);
//		set_intersection(Start_set.begin(), Start_set.end(), End_set.begin(), End_set.end(), inserter(Intersect_set, Intersect_set.end()));
//
//		set<int>::iterator it;
//		for (it = Intersect_set.begin(); it != Intersect_set.end(); it++) {
//			int j = *it;
//			lbcell++;
//			//if (LB_cell(&All_Query[i].Points, &All_Data[j].Points) > epsilon) continue;
//			if (LB_corner(&All_Query[i].Points, &All_Data[j].Points)) continue;
//			
//			lbcorner++;
//			//LB_band
//			bool flag = false;
//			for (int p = 0; p < All_Query[i].Points.size(); p = p + (All_Query[i].Points.size() / 50) ) {
//				All_set.clear();
//				Range_KDsearch(&All_KD[j], &All_set, All_Query[i].Points[p].latitude, All_Query[i].Points[p].longitude, All_KDroot[j]);
//				if (All_set.empty()) { flag = true; break; }
//			}
//			if (flag) continue;
//			lbband++;
//
//			if ((All_Query[i].Points.size() < MatrixLen) && (All_Query[j].Points.size() < MatrixLen)) {
//				if (DFS_DFD_Matrix(&All_Query[i].Points, &All_Data[j].Points)) answer.push_back(j);
//				for (int k = 0; k <= tot; k++) RepMatrix[Searchedi[k]][Searchedj[k]] = false;
//			}
//			else {
//				bool bdis = EPplus_DFD(&All_Query[i].Points, &All_Data[j].Points);
//				if (bdis) answer.push_back(j);
//			}
//		}
//		anssum = anssum + answer.size(); // the similarity amount of all query
//	}
//	fout << "LB_cell: " << lbcell << endl;
//	fout << "LB_corner: " << lbcorner << endl;
//	fout << "LB_band: " << lbband << endl;
//	fout << "Index answer: " << anssum << endl;
//}

//void similarity_search_Index_EPRseg() {
//	int anssum = 0;
//	int lbcell = 0, lbcorner = 0, lbband = 0;
//	memset(RepMatrix, false, sizeof(RepMatrix));
//	set<int> Start_set, End_set, Intersect_set, All_set, res_set;
//	for (int i = 0; i < All_Query.size(); i++) {
//		answer.clear();
//		//LB_cell
//		Start_set.clear(); End_set.clear(); Intersect_set.clear();
//		Range_KDsearch(&Start_KD, &Start_set, All_Query[i].Points[0].latitude, All_Query[i].Points[0].longitude, root_Start);
//		Range_KDsearch(&End_KD, &End_set, All_Query[i].Points.back().latitude, All_Query[i].Points.back().longitude, root_End);
//		set_intersection(Start_set.begin(), Start_set.end(), End_set.begin(), End_set.end(), inserter(Intersect_set, Intersect_set.end()));
//
//		set<int>::iterator it;
//		for (it = Intersect_set.begin(); it != Intersect_set.end(); it++) {
//			int j = *it;
//			lbcell++;
//			//if (LB_cell(&All_Query[i].Points, &All_Data[j].Points) > epsilon) continue;
//			if (LB_corner(&All_Query[i].Points, &All_Data[j].Points)) continue;
//
//			lbcorner++;
//			//LB_band
//			bool flag = false;
//			//for (int p = 0; p < All_Query[i].Points.size(); p = p + (All_Query[i].Points.size() / 50 + 1)) {
//			for (int p = 0; p < All_Query[i].Points.size(); p++) {
//				All_set.clear();
//				if (!Range_KDsearch_forband(&All_KD[j], All_Query[i].Points[p].latitude, All_Query[i].Points[p].longitude, All_KDroot[j])) { flag = true; break; }
//			}
//			if (flag) continue;
//			flag = false;
//			for (int p = 0; p < All_Data[j].Points.size(); p++) {
//				All_set.clear();
//				if (!Range_KDsearch_forband(&All_KD[i], All_Data[j].Points[p].latitude, All_Data[j].Points[p].longitude, All_KDroot[i])) { flag = true; break; }
//			}
//			if (flag) continue;
//			lbband++;
//			bool bdis = EPplus_RSeg_DFD(&All_Query[i].Points, &All_Data[j].Points, i, j);
//			if (bdis) answer.push_back(j);
//		}
//		anssum = anssum + answer.size(); // the similarity amount of all query
//	}
//	//fout << "LB_cell: " << lbcell << endl;
//	//fout << "LB_corner: " << lbcorner << endl;
//	//fout << "LB_band: " << lbband << endl;
//	fout << "EP_Rseg answer: " << anssum << endl;
//}

//void similarity_search_Index_EPCseg() {
//	int anssum = 0;
//	int lbcell = 0, lbcorner = 0, lbband = 0;
//	memset(RepMatrix, false, sizeof(RepMatrix));
//	set<int> Start_set, End_set, Intersect_set, All_set, res_set;
//	for (int i = 0; i < All_Query.size(); i++) {
//		answer.clear();
//		//LB_cell
//		Start_set.clear(); End_set.clear(); Intersect_set.clear();
//		Range_KDsearch(&Start_KD, &Start_set, All_Query[i].Points[0].latitude, All_Query[i].Points[0].longitude, root_Start);
//		Range_KDsearch(&End_KD, &End_set, All_Query[i].Points.back().latitude, All_Query[i].Points.back().longitude, root_End);
//		set_intersection(Start_set.begin(), Start_set.end(), End_set.begin(), End_set.end(), inserter(Intersect_set, Intersect_set.end()));
//
//		set<int>::iterator it;
//		for (it = Intersect_set.begin(); it != Intersect_set.end(); it++) {
//			int j = *it;
//			lbcell++;
//			//if (LB_cell(&All_Query[i].Points, &All_Data[j].Points) > epsilon) continue;
//			if (LB_corner(&All_Query[i].Points, &All_Data[j].Points)) continue;
//
//			lbcorner++;
//			//LB_band
//			bool flag = false;
//			//for (int p = 0; p < All_Query[i].Points.size(); p = p + (All_Query[i].Points.size() / 50 + 1)) {
//			for (int p = 0; p < All_Query[i].Points.size(); p++) {
//				All_set.clear();
//				if (!Range_KDsearch_forband(&All_KD[j], All_Query[i].Points[p].latitude, All_Query[i].Points[p].longitude, All_KDroot[j])) { flag = true; break; }
//			}
//			if (flag) continue;
//			flag = false;
//			for (int p = 0; p < All_Data[j].Points.size(); p++) {
//				All_set.clear();
//				if (!Range_KDsearch_forband(&All_KD[i], All_Data[j].Points[p].latitude, All_Data[j].Points[p].longitude, All_KDroot[i])) { flag = true; break; }
//			}
//			if (flag) continue;
//			lbband++;
//			bool bdis = EPplus_CSeg_DFD(&All_Query[i].Points, &All_Data[j].Points, i, j);
//			if (bdis) answer.push_back(j);
//		}
//		anssum = anssum + answer.size(); // the similarity amount of all query
//	}
//	//fout << "LB_cell: " << lbcell << endl;
//	//fout << "LB_corner: " << lbcorner << endl;
//	//fout << "LB_band: " << lbband << endl;
//	fout << "EP_Cseg answer: " << anssum << endl;
//}

void Get_exact_DFD() {
	//set<int> Start_set, End_set, Intersect_set, All_set, res_set;
	//for (int i = 0; i < All_Query.size(); i++) {
	//	answer.clear();
	//	//LB_cell
	//	Start_set.clear(); End_set.clear(); Intersect_set.clear();
	//	Range_KDsearch(&Start_KD, &Start_set, All_Query[i].Points[0].latitude, All_Query[i].Points[0].longitude, root_Start);
	//	Range_KDsearch(&End_KD, &End_set, All_Query[i].Points.back().latitude, All_Query[i].Points.back().longitude, root_End);
	//	set_intersection(Start_set.begin(), Start_set.end(), End_set.begin(), End_set.end(), inserter(Intersect_set, Intersect_set.end()));

	//	set<int>::iterator it;
	//	for (it = Intersect_set.begin(); it != Intersect_set.end(); it++) {
	//		int j = *it;
	//		if (i == j) continue;
	//		double dfd = double_DFD(&All_Query[i].Points, &All_Data[j].Points);
	//		fout << fixed;
	//		fout << i << ' ' << j << ' ' << ' ' << setprecision(10) <<dfd << endl;
	//	}
	//}

	// all pairs
	for (int i = 0; i < All_Query.size(); i++) {
		for (int j = i + 1; j < All_Query.size(); j++) {
			double dfd= double_DFD(&All_Query[i].Points, &All_Data[j].Points);
			fout << fixed;
			fout << i << ' ' << j << ' ' << ' ' << setprecision(10) << dfd << endl;
		}
	}
}
#endif
