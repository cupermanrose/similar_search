#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <windows.h>

#include <BLtoDist.h>
#include <Init.h>
#include <KD_tree.h>
#include <mytime.h>
#include <dfd.h>
#include <similarity.h>
#include <M_tree.h>
#include <R_tree.h>
#include <M_tree_bulkload.h>

using namespace std;

int main(int argc, char **argv) {

		init_time();
		//epsilon = 0.001;//0.1km
		epsilon = 0.01;//1km
		init(argv[1]);
		
		fout << "query&data files= 000.." << query_number << endl;
		fout << "TestNumber= " << All_Query.size() << endl;
		fout << "epsilon= " << epsilon << endl;
		//init_KD();
		out_time("Init ");

		/*{
			init_time();

			similarity_search();

			out_time("brute force: ");
		}*/

		/*{
			init_time();

			similarity_search_baseline();

			out_time("baseline: ");
		}*/

		{
			init_time();

			similarity_search_BLGroup();

			out_time("BLGroup: ");
		}

		{
			init_time();
			similarity_search_mtree();
			out_time("mtree: ");
		}

		{
			init_time();
			similarity_search_mtreeBL();
			out_time("mtreeBL: ");
		}

		/*{
			init_time();
			similarity_search_rtree();
			out_time("rtree: ");
		}*/


		{
			init_time();
			init_KD();
			out_time("init KD: ");

			init_time();

			similarity_search_triangle();

			out_time("triangle: ");
		}
		

	/*{
		epsilon = 0.05;

		init_time();
		cnt = 0;
		for (int i = 0; i < All_Query.size(); i++) {
			int last = cnt;
			for (int j = 0; j < All_Data.size(); j++) {
				double ff = LB_cell(&All_Query[i].Points, &All_Data[j].Points);
				if (ff < epsilon) cnt++;
			}
		}
		cout << cnt << endl;
		out_time("LBcell");

		init_time();
		cnt = 0;
		set<int> Start_set, End_set, Intersect_set;
		for (int i = 0; i < All_Query.size(); i++) {
			answer.clear();
			Start_set.clear(); End_set.clear(); Intersect_set.clear();
			Range_KDsearch(Start_KD, Start_set, All_Query[i].Points[0].latitude, All_Query[i].Points[0].longitude, root_Start);
			Range_KDsearch(End_KD, End_set, All_Query[i].Points.back().latitude, All_Query[i].Points.back().longitude, root_End);

			set_intersection(Start_set.begin(), Start_set.end(), End_set.begin(), End_set.end(), inserter(Intersect_set, Intersect_set.end()));
			cnt = cnt + Intersect_set.size();
		}
		cout << cnt << endl;
		out_time("KDindex");

		system("pause");
	}*/	

	cout << "finish" << endl;
	system("pause");
}