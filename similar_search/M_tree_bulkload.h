#ifndef _MTREEBULKLOAD_H
#define _MTREEBULKLOAD_H

#include <iostream> 
#include <algorithm>
#include <fstream>
#include <vector>
#include <time.h> 
#include <BLtoDist.h>
#include <grouping.h>

using namespace std;

namespace MtreeBulkLoad {
	const int Capacity = 20;
	const int LBCapacity = Capacity / 2;

	struct Node;

	struct Entry {
		int Tid, son_node; // Tid: feature value of the object, son_node: offset from the start address
		double dis_p, radius; // dis_p: distance with parent object, radius: cover radius
	};

	struct Node {
		Entry entries[Capacity];
		int entry_num, ID, height, CenterT; // entry number in the node; CenterT is the center; ID: the tree ID of the node
		bool leaf, root; // leaf and root label
	};

	vector<Node> Tree;
	vector<int> Candidate, Answer;
	int root, DisNum;


	double GetDistance(Entry& A, Entry& B) { // distance
		double dfd = double_DFD(All_Data[A.Tid].Points, All_Data[B.Tid].Points);
		return dfd;
	}

	double GetDisWithCenter(int CenterT, Entry& B) { // distance
		double dfd = double_DFD(All_Data[CenterT].Points, All_Data[B.Tid].Points);
		return dfd;
	}

	void AddAllEntry(Node& CurNode) { // add to the answer directly
		if (!CurNode.leaf) {
			for (int i = 0; i < CurNode.entry_num; i++) {
				AddAllEntry(Tree[CurNode.entries[i].son_node]);
			}
		}
		else {
			for (int i = 0; i < CurNode.entry_num; i++) {
				//Candidate.push_back(CurNode.entries[i].Tid);
				Answer.push_back(CurNode.entries[i].Tid);
			}
		}
	}
	void RangeQueryMemory(Node& CurNode, Entry& Q, double DisPQ) {
		if (!CurNode.leaf) {
			for (int i = 0; i < CurNode.entry_num; i++) {
				if (abs(DisPQ - CurNode.entries[i].dis_p) <= (Q.radius + CurNode.entries[i].radius)) {
					DisNum++;
					double DisRQ = GetDistance(Q, CurNode.entries[i]);
					if (DisRQ + CurNode.entries[i].radius <= Q.radius) { // Q include CurNode.entries[i]
						AddAllEntry(Tree[CurNode.entries[i].son_node]);
					}
					else {
						if (DisRQ <= Q.radius + CurNode.entries[i].radius) {
							RangeQueryMemory(Tree[CurNode.entries[i].son_node], Q, DisRQ);
						}
					}
				}
			}
		}
		else {
			for (int i = 0; i < CurNode.entry_num; i++) {
				if (abs(DisPQ - CurNode.entries[i].dis_p) <= Q.radius) {
					//double DisRQ = GetDistance(Q, CurNode.entries[i]);
					Candidate.push_back(CurNode.entries[i].Tid);
				}
			}
		}
	}
	
	void RangeQueryLoose(Node& CurNode, Entry& Q, double LBPQ, double UBPQ) {
		if (!CurNode.leaf) {
			for (int i = 0; i < CurNode.entry_num; i++) {
				if (((LBPQ - CurNode.entries[i].dis_p) <= (Q.radius + CurNode.entries[i].radius)) &&
					((CurNode.entries[i].dis_p - UBPQ) <= (Q.radius + CurNode.entries[i].radius))) {
					DisNum++;
					double LBRQ = grouping::GlbDFD(grouping::AllGTra[Q.Tid], grouping::AllGTra[CurNode.entries[i].Tid]);
					double UBRQ = grouping::GubDFD(grouping::AllGTra[Q.Tid], grouping::AllGTra[CurNode.entries[i].Tid]);

					if (UBRQ + CurNode.entries[i].radius <= Q.radius) { // Q include CurNode.entries[i]
						AddAllEntry(Tree[CurNode.entries[i].son_node]);
					}
					else {
						if (LBRQ <= Q.radius + CurNode.entries[i].radius) {
							RangeQueryLoose(Tree[CurNode.entries[i].son_node], Q, LBRQ, UBRQ);
						}
					}
				}
			}
		}
		else {
			for (int i = 0; i < CurNode.entry_num; i++) {
				if (((LBPQ - CurNode.entries[i].dis_p) <= Q.radius) &&
					((CurNode.entries[i].dis_p - UBPQ) <= Q.radius)) {
					//double DisRQ = GetDistance(Q, CurNode.entries[i]);
					Candidate.push_back(CurNode.entries[i].Tid);
				}
			}
		}
	}
	
	void CreateEntry(Entry& NewEntry, int Tid, int son_node, double dis_p, double radius) {
		NewEntry.Tid = Tid;
		NewEntry.son_node = son_node;
		NewEntry.dis_p = dis_p;
		NewEntry.radius = radius;
		return;
	}
	int CreateNode(int entry_num, int H, bool root, bool leaf) {
		Node NewNode;
		NewNode.entry_num = entry_num;
		NewNode.height = H;
		NewNode.root = root;
		NewNode.leaf = leaf;
		NewNode.ID = Tree.size();
		NewNode.CenterT = -1;
		Tree.push_back(NewNode);
		return NewNode.ID;
	}
	void Random_Kobject(vector<int>& RanSet, int Tsize, int Rsize) {
		/* initialize random seed: */
		srand(time(NULL));
		bool* flag = new bool[Tsize];
		for (int i = 0; i < Tsize; i++) flag[i] = false;
		RanSet.clear();
		while (RanSet.size() != Rsize) {
			int ran = rand() % Tsize;
			if (!flag[ran]) {
				RanSet.push_back(ran); flag[ran] = true;
			}
		}
		delete[] flag;
		return;
	}
	void Assign_AllObject(vector<Entry>& S, vector<int>& fk, int* AssignG, int K) {
		
		Random_Kobject(fk, S.size(), K); // build sampling set

		for (int i = 0; i < S.size(); i++) { // assign to the nearest sampling
			int MinJ = -1; double MinD = INFINITY;
			for (int j = 0; j < fk.size(); j++) {
				if (i == fk[j]) {
					MinJ = fk[j]; break;
				}
				double dist = GetDistance(S[i], S[fk[j]]);
				if (dist < MinD) {
					MinJ = fk[j]; MinD = dist;
				}
			}
			AssignG[i] = MinJ;
		}

		for (int i = 0; i < fk.size(); i++) { // is fk[i]'s objects< LBCapacity remove fk[i]
			int cnt = 0;
			for (int j = 0; j < S.size(); j++) { // count 
				if (AssignG[j] == fk[i]) cnt++;
			}
			if (cnt < LBCapacity) {
				for (int ii = 0; ii < S.size(); ii++) { // add to another sampling
					if (AssignG[ii] == fk[i]) {
						int MinJ = -1; double MinD = INFINITY;
						for (int j = 0; j < fk.size(); j++) {
							if ((j == i) || (fk[j] == -1)) continue;
							double dist = GetDistance(S[ii], S[fk[j]]);
							if (dist < MinD) {
								MinJ = fk[j]; MinD = dist;
							}
						}
						AssignG[ii] = MinJ;
					}
				}
				fk[i] = -1; // remove later
			}
		}

		for (vector<int>::iterator It = fk.begin(); It != fk.end(); ) {
			if (*It == -1) {
				It = fk.erase(It);
			}
			else It++;
		}

		if (fk.size() == 1) Assign_AllObject(S, fk, AssignG, K);
		return;
	}
	void SpiltHighTree(vector<Entry>& F, int Now, int H) {
		if (Tree[Now].height == H + 1) {
			for (int i = 0; i < Tree[Now].entry_num; i++) {
				F.push_back(Tree[Now].entries[i]);
			}
			return;
		}
		for (int i = 0; i < Tree[Now].entry_num; i++) {
			SpiltHighTree(F, Tree[Now].entries[i].son_node, H);
		}
		return;
	}

	int BulkLoading(vector<Entry>& S, int H, bool root, bool leaf) { // S-> object-set; H-> height(hmin); root/leaf of the whole tree;
	
		if (S.size() <= Capacity) {
			int NewNode=CreateNode(S.size(), H + 1, root, leaf);
			for (int i = 0; i < S.size(); i++) Tree[NewNode].entries[i] = S[i];
			return NewNode;
		}
		vector<Entry> SuperTree; SuperTree.clear(); // one entry corresponds one node

		vector<int> fk; // fk: sampling set
		int* AssignG = new int[S.size()];	// AssignG: assign result
		Assign_AllObject(S, fk, AssignG, Capacity);
		 
		vector<Entry> Fi;
		for (int i = 0; i < fk.size(); i++) {
			Fi.clear();
			for (int j = 0; j < S.size(); j++) {
				if (AssignG[j] == fk[i]) Fi.push_back(S[j]);
			}
			S[fk[i]].son_node = BulkLoading(Fi, H, false, leaf); // sub-tree cannot be root
			Tree[S[fk[i]].son_node].CenterT = S[fk[i]].Tid;
			if (Tree[S[fk[i]].son_node].entry_num < LBCapacity) { // underfilled, split
				for (int j = 0; j < Tree[S[fk[i]].son_node].entry_num; j++) {
					SuperTree.push_back(Tree[S[fk[i]].son_node].entries[j]);
				}
			}
			else {
				SuperTree.push_back(S[fk[i]]);
			}
		}
		Fi.clear(); vector<Entry>(Fi).swap(Fi); // free memory
		delete[] AssignG;

		int Hmin = INT_MAX;
		for (int i = 0; i < SuperTree.size(); i++) {
			Hmin = min(Hmin, Tree[SuperTree[i].son_node].height);
		}

		vector<Entry> TempF;
		vector<Entry> F; F.clear();
		for (vector<Entry>::iterator It = SuperTree.begin(); It != SuperTree.end(); It++) {
			 TempF.clear();
			 if (Tree[It->son_node].height > Hmin) {
				 SpiltHighTree(TempF, It->son_node, Hmin);
				 for (int i = 0; i < TempF.size(); i++) {
					 F.push_back(TempF[i]);
				 }
			 }
			 else F.push_back(*It);
		}
		TempF.clear(); vector<Entry>(TempF).swap(TempF); // free memory
		SuperTree.clear(); vector<Entry>(SuperTree).swap(SuperTree); // free memory
		return BulkLoading(F, Hmin, root, false); // supertree cannot be leaf;
	}
	void UpdateRadius(int Now) {
		if (Tree[Now].leaf) {
			for (int i = 0; i < Tree[Now].entry_num; i++) {
				Tree[Now].entries[i].radius = 0;
				Tree[Now].entries[i].dis_p = GetDisWithCenter(Tree[Now].CenterT, Tree[Now].entries[i]);
			}
			return;
		}

		for (int i = 0; i < Tree[Now].entry_num; i++) {
			double MaxR = 0; int Son = Tree[Now].entries[i].son_node;
			UpdateRadius(Son);
			for (int j = 0; j < Tree[Son].entry_num; j++) {
				MaxR = max(MaxR, Tree[Son].entries[j].radius + Tree[Son].entries[j].dis_p);
			}
			Tree[Now].entries[i].radius = MaxR;
			if (!Tree[Now].root) Tree[Now].entries[i].dis_p = GetDisWithCenter(Tree[Now].CenterT, Tree[Now].entries[i]);
		}
		return;
	}
	void WriteToDisk(const char* filename) {
		fstream MTfile(filename, ios::out | ios::binary);
		MTfile.write((char *)&root, sizeof(int));
		for (int i = 0; i < Tree.size(); i++) {
			char* NewNode = reinterpret_cast<char *>(&Tree[i]);
			MTfile.write(NewNode, sizeof(Node));
		}
		MTfile.close();
		return;
	}

	void ReadFromDisk(const char* filename) {
		Tree.clear();
		fstream MTfile(filename, ios::in | ios::binary);
		MTfile.read((char*)&root, sizeof(int));
		int i = 0;
		while (!MTfile.eof()) {
			Node* NewNode = new Node;
			MTfile.read((char *)NewNode, sizeof(Node));
			Tree.push_back(*NewNode);
			i++;
		}
		cout << "read size: " << i << endl;
		MTfile.close();
		return;
	}
	void Build(const char* filename, int Datasize) { // a building example //char* filename
		Tree.clear();
		vector<Entry> S; S.clear();
		for (int i = 0; i < Datasize; i++) {
			Entry NewEntry;
			CreateEntry(NewEntry, i, -1, 0, 0);
			S.push_back(NewEntry);
		}
		root = BulkLoading(S, -1, true, true);
		S.clear(); vector<Entry>(S).swap(S); // free memory
		UpdateRadius(root);
		WriteToDisk(filename);
		cout << "MtreeBulkLoad root: " << root << endl;
		cout << "MtreeBulkLoad size: " << Tree.size() << endl;
	}
}


#endif