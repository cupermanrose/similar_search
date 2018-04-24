#ifndef _INIT_H
#define _INIT_H

using namespace std;

const int MaxLen = 100000;
const int data_number = 181;
const int query_number = 181;
const int TestNumber = 10000+100; // datasize: 10000 querysize:100
const int DataSize = 1000;
const bool addpoints = false;
fstream fout("G:\\work\\DFD_convoy\\experimence_results\\similarity_search\\Geolife\\test.txt", ios::out);

struct Trajectory {
	vector<Point> Points;
	int number;
	string filename;
};

struct IntTrajectory {
	vector<IntPoint> Points;
};

vector<IntTrajectory> Int_Data;
vector<IntTrajectory> Int_Query;

vector<Trajectory> All_Data;
vector<Trajectory> All_Query;
vector<Trajectory> TempSet;
vector<string> filepath;
//vector<pair<int, int>> Vstack;
//const int MatrixLen = 5000;
//bool RepMatrix[MatrixLen][MatrixLen];
//int Searchedi[MatrixLen * MatrixLen], Searchedj[MatrixLen * MatrixLen];
long long cnt, tot, ave;
bool dfd_flag;
int root_Start, root_End;
int add_len;
double Length_ave;
double lbj, ubj;

int segLen; // update the last element that inside or outside
int seg_num[105]; //temp array
vector<vector<pair<int, double>>> ExactDFD;
bool DisRecord[MaxLen]; // store distances of Q_i (DFD.h, grouping.h)

string tostring(int i) {
	string temps = to_string(i);
	while (temps.length() < 3) temps = '0' + temps;
	return (temps);
} // add pre zeros
void FindAll(const char* index, vector<string>* filepath) {
	std::string sindex = index; 
	WIN32_FIND_DATA ffd;
	HANDLE hFind = INVALID_HANDLE_VALUE;
	hFind = FindFirstFile((sindex+"*.plt").c_str(), &ffd);
	while (true) {
		(*filepath).push_back(sindex + ffd.cFileName);
		if (FindNextFile(hFind, &ffd) == 0) break;
	}
}

void init(string preffix) {
	
	Point tempP;

	TempSet.clear();
	for (int i = 0; i < query_number; i++) {	// first 10 people as query data
		string index = preffix + tostring(i) + "\\Trajectory\\"; // folder's named rule
		filepath.clear();
		FindAll(index.c_str(), &filepath); // find all files in this folder

		for (int j = 0; j < filepath.size(); j++) { // traverse all files in this folder
			if (TempSet.size() == TestNumber) break;
			fstream fin(filepath[j].c_str(), ios::in);
			if (!fin.is_open()) { fin.close(); continue; } // ignore invalid files
			Trajectory tempTra; tempTra.Points.clear(); //generate a trajectory
			tempTra.number = i; tempTra.filename = filepath[j];
			char temp[256];
			for (int k = 0; k < 6; k++) fin.getline(temp, 256);//ignore first 6 lines
			while (!fin.eof()) {
				double latitude, longitude; char ch;
				fin >> latitude >> ch >> longitude; fin.getline(temp, 256);// read b&l, and ignore the rest of a line
				toRad(tempP, latitude, longitude);
				tempTra.Points.push_back(tempP);
			}
			fin.close();

			if (addpoints) {
				Trajectory tempTra_add; tempTra_add.Points.clear(); //generate a new trajectory
				tempTra_add.number = i; tempTra_add.filename = filepath[j];
				tempTra_add.Points.push_back(tempTra.Points[0]);

				int addnumber = add_len / tempTra.Points.size();

				for (int k = 1; k < tempTra.Points.size(); k++) {
					if (addnumber > 1) {
						double platitude = tempTra.Points[k - 1].latitude;
						double plongitude = tempTra.Points[k - 1].longitude;
						double latitude_offset = (tempTra.Points[k].latitude - platitude) / addnumber;
						double longitude_offset = (tempTra.Points[k].longitude - plongitude) / addnumber;
						for (int kk = 1; kk < addnumber; kk++) {
							platitude = platitude + latitude_offset;
							plongitude = plongitude + longitude_offset;
							toRad(tempP, platitude, plongitude);
							tempTra_add.Points.push_back(tempP);
						}
					}
					tempTra_add.Points.push_back(tempTra.Points[k]);
				}
				TempSet.push_back(tempTra_add);
			}
			else TempSet.push_back(tempTra);
		}
	}

	////	read all query first and store them in All_Query
	//All_Query.clear();

	//for (int i = 0; i < query_number; i++) {	// first 10 people as query data
	//	string index = preffix + tostring(i) + "\\Trajectory\\"; // folder's named rule
	//	filepath.clear();
	//	FindAll(index.c_str(), &filepath); // find all files in this folder

	//	for (int j = 0; j < filepath.size(); j++) { // traverse all files in this folder
	//		if (All_Query.size() == TestNumber) break;
	//		fstream fin(filepath[j].c_str(), ios::in);
	//		if (!fin.is_open()) { fin.close(); continue; } // ignore invalid files
	//		Trajectory tempTra; tempTra.Points.clear(); //generate a trajectory
	//		tempTra.number = i; tempTra.filename = filepath[j];
	//		char temp[256];
	//		for (int k = 0; k < 6; k++) fin.getline(temp, 256);//ignore first 6 lines
	//		while (!fin.eof()) {
	//			double latitude, longitude; char ch;
	//			fin >> latitude >> ch >> longitude; fin.getline(temp, 256);// read b&l, and ignore the rest of a line
	//			toRad(tempP, latitude, longitude);
	//			tempTra.Points.push_back(tempP);
	//		}
	//		fin.close();

	//		if (addpoints) {
	//			Trajectory tempTra_add; tempTra_add.Points.clear(); //generate a new trajectory
	//			tempTra_add.number = i; tempTra_add.filename = filepath[j];
	//			tempTra_add.Points.push_back(tempTra.Points[0]);

	//			int addnumber = add_len / tempTra.Points.size();
	//		
	//			for (int k = 1; k < tempTra.Points.size(); k++) {
	//				if (addnumber > 1) {
	//					double platitude = tempTra.Points[k - 1].latitude;
	//					double plongitude = tempTra.Points[k - 1].longitude;
	//					double latitude_offset = (tempTra.Points[k].latitude - platitude) / addnumber;
	//					double longitude_offset = (tempTra.Points[k].longitude - plongitude) / addnumber;
	//					for (int kk = 1; kk < addnumber; kk++) {
	//						platitude = platitude + latitude_offset;
	//						plongitude = plongitude + longitude_offset;
	//						toRad(tempP, platitude, plongitude);
	//						tempTra_add.Points.push_back(tempP);
	//					}
	//				}
	//				tempTra_add.Points.push_back(tempTra.Points[k]);
	//			}
	//			All_Query.push_back(tempTra_add);
	//		}
	//		else All_Query.push_back(tempTra);
	//	}
	//}

	//// then read all data and store them in All_Data

	//All_Data.clear();

	//for (int i =0; i < data_number; i++) {	// total 182 folder, the rest of people's trajectories
	//	string index = preffix + tostring(i) + "\\Trajectory\\"; // fold's named rule
	//	filepath.clear();
	//	FindAll(index.c_str(), &filepath); // find all files in this folder

	//	for (int j = 0; j < filepath.size(); j++) { // traverse all files in this folder
	//		if (All_Data.size() == TestNumber) break;
	//		fstream fin(filepath[j].c_str(), ios::in);
	//		if (!fin.is_open()) { fin.close(); continue; } // ignore invalid files
	//		Trajectory tempTra; tempTra.Points.clear(); //generate a trajectory
	//		tempTra.number = i; tempTra.filename = filepath[j];
	//		char temp[256];
	//		for (int k = 0; k < 6; k++) fin.getline(temp, 256);//ignore first 6 lines
	//		while (!fin.eof()) {
	//			double latitude, longitude; char ch;
	//			fin >> latitude >> ch >> longitude; fin.getline(temp, 256);// read b&l, and ignore the rest of a line

	//			toRad(tempP, latitude, longitude);
	//			tempTra.Points.push_back(tempP);
	//		}
	//		fin.close();

	//		if (addpoints) {
	//			Trajectory tempTra_add; tempTra_add.Points.clear(); //generate a new trajectory
	//			tempTra_add.number = i; tempTra_add.filename = filepath[j];
	//			tempTra_add.Points.push_back(tempTra.Points[0]);

	//			int addnumber = add_len / tempTra.Points.size();

	//			for (int k = 1; k < tempTra.Points.size(); k++) {
	//				if (addnumber > 1) {
	//					double platitude = tempTra.Points[k - 1].latitude;
	//					double plongitude = tempTra.Points[k - 1].longitude;
	//					double latitude_offset = (tempTra.Points[k].latitude - platitude) / addnumber;
	//					double longitude_offset = (tempTra.Points[k].longitude - plongitude) / addnumber;
	//					for (int kk = 1; kk < addnumber; kk++) {
	//						platitude = platitude + latitude_offset;
	//						plongitude = plongitude + longitude_offset;
	//						toRad(tempP, platitude, plongitude);
	//						tempTra_add.Points.push_back(tempP);
	//					}
	//				}
	//				tempTra_add.Points.push_back(tempTra.Points[k]);
	//			}
	//			All_Data.push_back(tempTra_add);
	//		}
	//		else All_Data.push_back(tempTra);
	//	}
	//}

	// random choose 100 query in All_Query
	All_Data.clear();
	All_Query.clear();
	for (int i = 0; i < TestNumber; i++) {
		if ((i % 100 == 0) && (All_Query.size() < 100)) {
			All_Query.push_back(TempSet[i]);
		}
		else {
			if (All_Data.size() == DataSize) continue;
			All_Data.push_back(TempSet[i]);
		}
	}

	cout << All_Data.size() << endl;
	cout << All_Query.size() << endl;
	//{//read exact dfd value
	//	fstream fin("G:\\work\\DFD_convoy\\experimence_results\\similarity_search\\Geolife\\exact_DFD.txt", ios::in);
	//	//fstream fin("G:\\work\\DFD_convoy\\experimence_results\\similarity_search\\Geolife\\Exact_1000_50.txt", ios::in);
	//	ExactDFD.clear(); int prex = -1;
	//	vector<pair<int, double>> tempvec; tempvec.clear();
	//	for (int i = 0; i < 19000; i++) {
	//	//for (int i = 0; i < All_Data.size(); i++) {
	//		ExactDFD.push_back(tempvec);
	//	}
	//	while (!fin.eof()) {
	//		int x, y;
	//		double dis;
	//		fin >> x >> y >> dis;
	//		pair<int, double> temppair = make_pair(y, dis);
	//		ExactDFD[x].push_back(temppair);
	//	}
	//	fin.close();
	//}
	//

	//Length_ave = 0.0;
	//for (int i = 0; i < All_Data.size(); i++) {
	//	Length_ave = Length_ave + All_Data[i].Points.size();
	//}
	//Length_ave = Length_ave / All_Data.size();
}

#endif
