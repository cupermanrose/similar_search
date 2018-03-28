#ifndef _SIMILARITY_H
#define _SIMILARITY_H

#include <algorithm>
#include <iterator>
#include <iomanip>
#include <M_tree.h>
#include <M_tree_bulkload.h>
#include <R_tree.h>
#include <grouping.h>

using namespace std;

set<int> Start_set, End_set, Intersect_set;
double lower[20000], upper[20000];
vector<int> answer;

void similarity_search() {
	int anssum = 0;
	for (int i = 0; i < All_Query.size(); i++) {
		answer.clear();
		for (int j = 0; j < All_Data.size(); j++) {
			double bdis = double_DFD(All_Query[i].Points, All_Data[j].Points);
			if (bdis < epsilon) answer.push_back(j);
		}
		anssum = anssum + answer.size();
	}
	fout << "Brute force answer: " << anssum << endl;
}

void similarity_search_mtreeBLLoose() {

	init_time();
	init_KD();
	out_time("init KD: ");

	init_time();
	grouping::Init();
	out_time("grouping Init: ");
	init_time();

	string filename = "G:\\work\\DFD_convoy\\experimence_results\\similarity_search\\Geolife\\MTreeBulkLoad_" + to_string(TestNumber) + "_" + to_string(MtreeBulkLoad::Capacity) + ".txt";
	/*init_time();
	MtreeBulkLoad::Build(filename.c_str(), All_Data.size());
	out_time("MtreeBulkLoad building: ");
	init_time();*/
	MtreeBulkLoad::ReadFromDisk(filename.c_str());

	int anssum = 0, lbcell = 0, ubmtree = 0, lbmtree = 0, lbband = 0, ubgreedy = 0, Grouplb = 0, Groupub = 0;
	MtreeBulkLoad::DisNum = 0;
	for (int i = 0; i < All_Query.size(); i++) {

		// Mtree Index
		MtreeBulkLoad::Answer.clear();
		MtreeBulkLoad::Candidate.clear();
		MtreeBulkLoad::Entry Q;
		MtreeBulkLoad::CreateEntry(Q, i, -1, 0, epsilon);
		MtreeBulkLoad::RangeQueryLoose(MtreeBulkLoad::Tree[MtreeBulkLoad::root], Q, -INFINITY, INFINITY);
		ubmtree = ubmtree + MtreeBulkLoad::Candidate.size() + MtreeBulkLoad::Answer.size();
		lbmtree = lbmtree + MtreeBulkLoad::Candidate.size();

		for (int j = 0; j < MtreeBulkLoad::Candidate.size(); j++) {

			int jj = MtreeBulkLoad::Candidate[j];
			//LB_cell
			if (LB_cell(All_Query[i].Points, All_Data[jj].Points) > epsilon) continue;
			lbcell++;

			//LB_band
			if (LB_band(i,jj)) continue;
			lbband++;

			//UB_greedy
			
			if (DFD_greedy(All_Query[i].Points, All_Data[jj].Points) < epsilon) { 
				MtreeBulkLoad::Answer.push_back(jj);
				continue;
			}
			ubgreedy++;

			//grouping
			grouping::GroupTra NewQ;
			grouping::DivideTrajectory(All_Query[i], i, NewQ);

			if (grouping::GlbDFD_LBrow(NewQ, grouping::AllGTra[jj]) > epsilon) continue;
			Grouplb++;

			if (grouping::GubDFD_LBrow(NewQ, grouping::AllGTra[jj]) < epsilon) {
				MtreeBulkLoad::Answer.push_back(jj);
				continue;
			}
			Groupub++;
				
			//double bdis = DFD_LBrow(All_Query[i].Points, All_Data[MtreeBulkLoad::Candidate[j]].Points);
			//double bdis = EP_DFD(All_Query[i].Points, All_Data[jj].Points); 
			//if (bdis < epsilon) MtreeBulkLoad::Answer.push_back(jj);
			if (EPplus_DFD(All_Query[i].Points, All_Data[jj].Points)) MtreeBulkLoad::Answer.push_back(jj);
		}
		anssum = anssum + MtreeBulkLoad::Answer.size();
	}
	fout << "MtreeBulkLoad::DisNumLoose: " << MtreeBulkLoad::DisNum << endl;
	fout << "ubmtreeLoose: " << ubmtree << endl;
	fout << "lbmtreeLoose: " << lbmtree << endl;
	fout << "lbcellLoose: " << lbcell << endl;
	fout << "lbbandLoose: " << lbband << endl;
	fout << "ubgreedyLoose: " << ubgreedy << endl;
	fout << "GrouplbLoose: " << Grouplb << endl;
	fout << "GroupubLoose: " << Groupub << endl;
	fout << "MtreeBulkLoadLoose answer: " << anssum << endl;
}

void similarity_search_baseline() {
	int anssum = 0, lbcell = 0;
	for (int i = 0; i < All_Query.size(); i++) {
		answer.clear();
		for (int j = 0; j < All_Data.size(); j++) {
			if (LB_cell(All_Query[i].Points, All_Data[j].Points) > epsilon) continue;
			lbcell++;
			double bdis = DFD_LBrow(All_Query[i].Points, All_Data[j].Points);
			if (bdis < epsilon) answer.push_back(j);
		}
		/*cout << i << ": ";
		for (int j = 0; j < answer.size(); j++) {
			cout << answer[j] << " ";
		}
		cout << endl;*/
		anssum = anssum + answer.size();
	}
	fout << "lbcell: " << lbcell << endl;
	fout << "baseline answer: " << anssum << endl;
}

void similarity_search_BLGroup() {
	init_time();
	grouping::Init();
	out_time("grouping Init: ");
	int anssum = 0, lbcell = 0, Grouplb = 0, Groupub = 0;
	for (int i = 0; i < All_Query.size(); i++) {
		answer.clear();
		for (int j = 0; j < All_Data.size(); j++) {
			if (LB_cell(All_Query[i].Points, All_Data[j].Points) > epsilon) continue;
			lbcell++;
			grouping::GroupTra NewQ;
			grouping::DivideTrajectory(All_Query[i], i, NewQ);

			if (grouping::GlbDFD_LBrow(NewQ, grouping::AllGTra[j]) > epsilon) continue;
			Grouplb++;

			if (grouping::GubDFD_LBrow(NewQ, grouping::AllGTra[j]) < epsilon) answer.push_back(j);
			else {
				Groupub++;
				double bdis = DFD_LBrow(All_Query[i].Points, All_Data[j].Points);
				if (bdis < epsilon) answer.push_back(j);
			}
		}
		/*cout << i << ": ";
		for (int j = 0; j < answer.size(); j++) {
		cout << answer[j] << " ";
		}
		cout << endl;*/
		anssum = anssum + answer.size();
	}
	fout << "lbcell: " << lbcell << endl;
	fout << "Grouplb: " << Grouplb << endl;
	fout << "Groupub: " << Groupub << endl;
	fout << "BLGroup answer: " << anssum << endl;
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
		Range_KDsearch(Start_KD, Start_set, All_Query[i].Points[0].latitude, All_Query[i].Points[0].longitude, root_Start);
		Range_KDsearch(End_KD, End_set, All_Query[i].Points.back().latitude, All_Query[i].Points.back().longitude, root_End);	
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
				if (!Range_KDsearch_forband(All_KD[j], All_Query[i].Points[p].latitude, All_Query[i].Points[p].longitude, All_KDroot[j])) { flag = true; break; }
			}
			if (flag) { update_lbub(j, lbj, ubj, dfdj); continue; }
			flag = false;
			for (int p = 0; p < All_Data[j].Points.size(); p++) {
				if (!Range_KDsearch_forband(All_KD[i], All_Data[j].Points[p].latitude, All_Data[j].Points[p].longitude, All_KDroot[i])) { flag = true; break; }
			}
			if (flag) { update_lbub(j, lbj, ubj, dfdj); continue; }
			lbband++;

			//UB_greedy
			double ub_temp = DFD_greedy(All_Query[i].Points, All_Data[j].Points);
			ubj = min(ubj, ub_temp);
			if (ub_temp < epsilon) { update_lbub(j, lbj, ubj, dfdj); answer.push_back(j); }
			else {
				ubgreedy++;
				dfdj= EP_DFD(All_Query[i].Points, All_Data[j].Points);
			//	dfdj = DFD_LBrow(All_Query[i].Points, All_Data[j].Points);
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

//void similarity_search_rtree() {
//
//	string filename = "G:\\work\\DFD_convoy\\experimence_results\\similarity_search\\Geolife\\RTree_" + to_string(TestNumber) + "_" + to_string(Rtree::Capacity) + ".txt";
//	init_time();
//	Rtree::Build(filename.c_str(), All_Data.size());
//	out_time("Rtree building: ");
//	init_time();
//	Rtree::ReadFromDisk(filename.c_str());
//
//	int anssum = 0, lbcell = 0, ubrtree = 0, lbrtree = 0;
//	for (int i = 0; i < All_Query.size(); i++) {
//		Rtree::Answer.clear();
//		Rtree::MBR MBR;
//		Rtree::CreateMBR(MBR, i);
//		Rtree::ExtendMBR(MBR, epsilon);
//		Rtree::Entry Q;
//		Rtree::CreateEntry(Q, i, -1, MBR);
//		Rtree::Candidate.clear();
//		Rtree::RangeQueryMemory(Rtree::Tree[Rtree::root], Q);
//		ubrtree = ubrtree + Rtree::Answer.size() + Rtree::Candidate.size();
//		lbrtree = lbrtree + Rtree::Candidate.size();
//		for (int j = 0; j < Rtree::Candidate.size(); j++) {
//			if (LB_cell(All_Query[i].Points, All_Data[Rtree::Candidate[j]].Points) > epsilon) continue;
//			lbcell++;
//			//double bdis = double_DFD(All_Query[i].Points, All_Data[Rtree::Candidate[j]].Points);
//			double bdis = DFD_LBrow(All_Query[i].Points, All_Data[Rtree::Candidate[j]].Points);
//			if (bdis < epsilon) Rtree::Answer.push_back(Rtree::Candidate[j]);
//		}
//		anssum = anssum + Rtree::Answer.size();
//	}
//	fout << "ubrtree: " << ubrtree << endl;
//	fout << "lbrtree: " << lbrtree << endl;
//	fout << "lbcell: " << lbcell << endl;
//	fout << "Rtree answer: " << anssum << endl;
//}

//void similarity_search_mtreeBL() {
//
//	string filename = "G:\\work\\DFD_convoy\\experimence_results\\similarity_search\\Geolife\\MTreeBulkLoad_" + to_string(TestNumber) + "_" + to_string(MtreeBulkLoad::Capacity) + ".txt";
//	/*init_time();
//	MtreeBulkLoad::Build(filename.c_str(), All_Data.size());
//	out_time("MtreeBulkLoad building: ");
//	init_time();*/
//	MtreeBulkLoad::ReadFromDisk(filename.c_str());
//
//	/*int res = 0;
//	for (int i = 0; i < Mtree::Tree.size()-1; i++) {
//	res=res+Mtree::Tree[i].entry_num;
//	}
//	cout << res / Mtree::Tree.size() << endl;
//	system("pause");*/
//
//	int anssum = 0, lbcell = 0, ubmtree = 0, lbmtree = 0;
//	MtreeBulkLoad::DisNum = 0;
//	for (int i = 0; i < All_Query.size(); i++) {
//		MtreeBulkLoad::Answer.clear();
//		MtreeBulkLoad::Candidate.clear();
//		MtreeBulkLoad::Entry Q;
//		MtreeBulkLoad::CreateEntry(Q, i, -1, 0, epsilon);
//		MtreeBulkLoad::RangeQueryMemory(MtreeBulkLoad::Tree[MtreeBulkLoad::root], Q, 0);
//		ubmtree = ubmtree + MtreeBulkLoad::Candidate.size() + MtreeBulkLoad::Answer.size();
//		lbmtree = lbmtree + MtreeBulkLoad::Candidate.size();
//		for (int j = 0; j < MtreeBulkLoad::Candidate.size(); j++) {
//			if (LB_cell(All_Query[i].Points, All_Data[MtreeBulkLoad::Candidate[j]].Points) > epsilon) continue;
//			lbcell++;
//			//double bdis = double_DFD(All_Query[i].Points, All_Data[Rtree::Candidate[j]].Points);
//			double bdis = DFD_LBrow(All_Query[i].Points, All_Data[MtreeBulkLoad::Candidate[j]].Points);
//			if (bdis < epsilon) MtreeBulkLoad::Answer.push_back(MtreeBulkLoad::Candidate[j]);
//		}
//		anssum = anssum + MtreeBulkLoad::Answer.size();
//	}
//	fout << "MtreeBulkLoad::DisNum: " << MtreeBulkLoad::DisNum << endl;
//	fout << "ubmtree: " << ubmtree << endl;
//	fout << "lbmtree: " << lbmtree << endl;
//	fout << "lbcell: " << lbcell << endl;
//	fout << "MtreeBulkLoad answer: " << anssum << endl;
//}

//void similarity_search_mtree() {
//
//	string filename = "G:\\work\\DFD_convoy\\experimence_results\\similarity_search\\Geolife\\MTree_" + to_string(TestNumber) + "_" + to_string(Mtree::Capacity) + "(balanced).txt";
//	init_time();
//	Mtree::Build(filename.c_str(), All_Data.size());
//	out_time("Mtree building: ");
//	init_time();
//	Mtree::ReadFromDisk(filename.c_str());
//
//	/*int res = 0;
//	for (int i = 0; i < Mtree::Tree.size()-1; i++) {
//		res=res+Mtree::Tree[i].entry_num;
//	}
//	cout << res / Mtree::Tree.size() << endl;
//	system("pause");*/
//
//	int anssum = 0, lbcell = 0, ubmtree = 0, lbmtree = 0;
//	Mtree::DisNum = 0;
//	for (int i = 0; i < All_Query.size(); i++) {
//		Mtree::Answer.clear();
//		Mtree::Candidate.clear();
//		Mtree::Entry Q;
//		Mtree::CreateEntry(Q, i, NULL, 0, epsilon);
//		Mtree::RangeQueryMemory(Mtree::Tree[Mtree::root], Q, 0);
//		ubmtree = ubmtree + Mtree::Candidate.size() + Mtree::Answer.size();
//		lbmtree = lbmtree + Mtree::Candidate.size();
//		for (int j = 0; j < Mtree::Candidate.size(); j++) {
//			if (LB_cell(All_Query[i].Points, All_Data[Mtree::Candidate[j]].Points) > epsilon) continue;
//			lbcell++;
//			//double bdis = double_DFD(All_Query[i].Points, All_Data[Rtree::Candidate[j]].Points);
//			double bdis = DFD_LBrow(All_Query[i].Points, All_Data[Mtree::Candidate[j]].Points);
//			if (bdis < epsilon) Mtree::Answer.push_back(Mtree::Candidate[j]);
//		}
//		anssum = anssum + Mtree::Answer.size();
//	}
//	fout << "Mtree::DisNum: " << Mtree::DisNum << endl;
//	fout << "ubmtree: " << ubmtree << endl;
//	fout << "lbmtree: " << lbmtree << endl;
//	fout << "lbcell: " << lbcell << endl;
//	fout << "Mtree answer: " << anssum << endl;
//}

//void similarity_search_Index_EPplus() {
//	int anssum = 0;
//	int lbcell = 0, lbcorner = 0, lbband = 0;
//	memset(RepMatrix, false, sizeof(RepMatrix));
//	set<int> Start_set, End_set, Intersect_set, All_set, res_set;
//	for (int i = 0; i < All_Query.size(); i++) {
//		answer.clear();
//		//LB_cell
//		Start_set.clear(); End_set.clear(); Intersect_set.clear();
//		Range_KDsearch(Start_KD, Start_set, All_Query[i].Points[0].latitude, All_Query[i].Points[0].longitude, root_Start);
//		Range_KDsearch(End_KD, End_set, All_Query[i].Points.back().latitude, All_Query[i].Points.back().longitude, root_End);
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
//				if (!Range_KDsearch_forband(All_KD[j], All_Query[i].Points[p].latitude, All_Query[i].Points[p].longitude, All_KDroot[j])) { flag = true; break; }
//			}
//			if (flag) continue;
//			flag = false;
//			for (int p = 0; p < All_Data[j].Points.size(); p++) {
//				All_set.clear();
//				if (!Range_KDsearch_forband(All_KD[i], All_Data[j].Points[p].latitude, All_Data[j].Points[p].longitude, All_KDroot[i])) { flag = true; break; }
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
//		Range_KDsearch(Start_KD, Start_set, All_Query[i].Points[0].latitude, All_Query[i].Points[0].longitude, root_Start);
//		Range_KDsearch(End_KD, End_set, All_Query[i].Points.back().latitude, All_Query[i].Points.back().longitude, root_End);
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
//				Range_KDsearch(All_KD[j], &All_set, All_Query[i].Points[p].latitude, All_Query[i].Points[p].longitude, All_KDroot[j]);
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
//		Range_KDsearch(Start_KD, Start_set, All_Query[i].Points[0].latitude, All_Query[i].Points[0].longitude, root_Start);
//		Range_KDsearch(End_KD, End_set, All_Query[i].Points.back().latitude, All_Query[i].Points.back().longitude, root_End);
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
//				if (!Range_KDsearch_forband(All_KD[j], All_Query[i].Points[p].latitude, All_Query[i].Points[p].longitude, All_KDroot[j])) { flag = true; break; }
//			}
//			if (flag) continue;
//			flag = false;
//			for (int p = 0; p < All_Data[j].Points.size(); p++) {
//				All_set.clear();
//				if (!Range_KDsearch_forband(All_KD[i], All_Data[j].Points[p].latitude, All_Data[j].Points[p].longitude, All_KDroot[i])) { flag = true; break; }
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
//		Range_KDsearch(Start_KD, Start_set, All_Query[i].Points[0].latitude, All_Query[i].Points[0].longitude, root_Start);
//		Range_KDsearch(End_KD, End_set, All_Query[i].Points.back().latitude, All_Query[i].Points.back().longitude, root_End);
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
//				if (!Range_KDsearch_forband(All_KD[j], All_Query[i].Points[p].latitude, All_Query[i].Points[p].longitude, All_KDroot[j])) { flag = true; break; }
//			}
//			if (flag) continue;
//			flag = false;
//			for (int p = 0; p < All_Data[j].Points.size(); p++) {
//				All_set.clear();
//				if (!Range_KDsearch_forband(All_KD[i], All_Data[j].Points[p].latitude, All_Data[j].Points[p].longitude, All_KDroot[i])) { flag = true; break; }
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
	//	Range_KDsearch(Start_KD, Start_set, All_Query[i].Points[0].latitude, All_Query[i].Points[0].longitude, root_Start);
	//	Range_KDsearch(End_KD, End_set, All_Query[i].Points.back().latitude, All_Query[i].Points.back().longitude, root_End);
	//	set_intersection(Start_set.begin(), Start_set.end(), End_set.begin(), End_set.end(), inserter(Intersect_set, Intersect_set.end()));

	//	set<int>::iterator it;
	//	for (it = Intersect_set.begin(); it != Intersect_set.end(); it++) {
	//		int j = *it;
	//		if (i == j) continue;
	//		double dfd = double_DFD(All_Query[i].Points, All_Data[j].Points);
	//		fout << fixed;
	//		fout << i << ' ' << j << ' ' << ' ' << setprecision(10) <<dfd << endl;
	//	}
	//}

	// all pairs
	/*for (int i = 0; i < All_Query.size(); i++) {
		for (int j = i + 1; j < All_Query.size(); j++) {
			double dfd= double_DFD(All_Query[i].Points, All_Data[j].Points);
			fout << fixed;
			fout << i << ' ' << j << ' ' << ' ' << setprecision(10) << dfd << endl;
		}
	}*/

	int Knear = 50;
	int Tj[100];
	double TjD[100];
	for (int i = 0; i < All_Query.size(); i++) {
		for (int j = 0; j < Knear; j++) {
			Tj[j] = -1;
			TjD[j] = INFINITY;
		}

		for (int j = 0; j < All_Query.size(); j++) {
			double DisP1 = double_dist(All_Query[j].Points[0], All_Query[i].Points[0]);
			double DisP2 = double_dist(All_Query[j].Points[All_Query[j].Points.size()-1], All_Query[i].Points[All_Query[i].Points.size()-1]);
			for (int k = 0; k < Knear; k++) {
				if (max(DisP1, DisP2) < TjD[k]) {
					Tj[k] = j;
					TjD[k] = max(DisP1, DisP2);
					break;
				}
			}
		}

		for (int j = 0; j < Knear; j++) {
			double dfd= double_DFD(All_Query[i].Points, All_Data[Tj[j]].Points);
			fout << fixed;
			fout << i << ' ' << Tj[j] << ' ' << ' ' << setprecision(10) << dfd << endl;
		}
	}

}
#endif
