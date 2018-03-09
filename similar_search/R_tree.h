#ifndef _RTREE_H
#define _RTREE_H

#include <iostream> 
#include <algorithm>
#include <fstream>
#include <vector>
#include <BLtoDist.h>
#include <init.h>

using namespace std;

namespace Rtree {
	
	const int dimension = 2;
	const int Capacity = 20;

	struct MBR {
		double U[dimension], L[dimension];
	};

	struct Entry {
		int TID, SonNode;
		MBR MBR;
	};

	struct Node {
		Entry entries[Capacity];
		int ID, PID, entry_num;
		bool leaf, root;
	};

	vector<Node> Tree;
	vector<Entry> NN, N1, N2;
	vector<int> Candidate;
	int root, cnt;
	
	void ClearMBR(MBR& MBR) {
		for (int i = 0; i < dimension; i++) {
			MBR.L[i] = INFINITY;
			MBR.U[i] = -INFINITY;
		}
		return;
	}
	
	void CreateMBR(MBR& MBR, int TID) {
		ClearMBR(MBR);
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < All_Data[TID].Points.size(); j++) {
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
	
	void ExtendMBR(MBR& MBR, double d) {
		for (int i = 0; i < dimension; i++) {
			MBR.L[i] = MBR.L[i] - d;
			MBR.U[i] = MBR.U[i] + d;
		}
		return;
	}
	
	bool EmbodyMBR(MBR& M1, MBR& M2) { // if M2 in M1
		for (int i = 0; i < dimension; i++) {
			if (M2.L[i] < M1.L[i] - eps) return false;
			if (M2.U[i] > M1.U[i] + eps) return false;
		}
		return true;
	}
	
	void CreateEntry(Entry& NewEntry, int TID, int SonNode, MBR MBR) {
		NewEntry.TID = TID;
		NewEntry.SonNode = SonNode;
		NewEntry.MBR = MBR;
		return;
	}
	
	void CreateNode(int PID, int entry_num, bool leaf, bool root) {
		Node* NewNode = new Node;
		NewNode->ID = Tree.size();
		NewNode->PID = PID;
		NewNode->entry_num = entry_num;
		NewNode->leaf = leaf;
		NewNode->root = root;
		Tree.push_back(*NewNode);
		delete NewNode;
		return;
	}
	
	void UpdateMBR(int N) {
		if (Tree[N].root) return; // root doesn't have parent
		int PID = Tree[N].PID;
		for (int i = 0; i < Tree[PID].entry_num; i++) {
			if (Tree[N].ID == Tree[PID].entries[i].SonNode) {
				ClearMBR(Tree[PID].entries[i].MBR);
				for (int j = 0; j < Tree[N].entry_num; j++) {
					for (int k = 0; k < dimension; k++) {
						Tree[PID].entries[i].MBR.L[k] = min(Tree[PID].entries[i].MBR.L[k], Tree[N].entries[j].MBR.L[k]);
						Tree[PID].entries[i].MBR.U[k] = max(Tree[PID].entries[i].MBR.U[k], Tree[N].entries[j].MBR.U[k]);
					}
				}
				return;
			}
		}
		return;
	}
	
	void GetMerge(MBR& MergeMBR, MBR& M1, MBR& M2) {
		for (int i = 0; i < dimension; i++) {
			MergeMBR.L[i] = min(M1.L[i], M2.L[i]);
			MergeMBR.U[i] = max(M1.U[i], M2.U[i]);
		}
		return;
	}
	
	double GetArea(MBR& MBR) {
		double Area= 1;
		for (int i = 0; i < dimension; i++) {
			Area = Area*(MBR.L[i] - MBR.U[i]);
		}
		return abs(Area);
	}
	
	bool GetIntersection(MBR& M1, MBR& M2) {
		double eps = 1.0E-10;
		for (int i = 0; i < dimension; i++) {
			if (M1.L[i] > M2.U[i] + eps) return false;
			if (M2.L[i] > M1.U[i] + eps) return false;
		}
		return true;
	}
	
	pair<int, int> Promote(int N) {
		pair<int, int> temppair;
		double MaxArea = -INFINITY;
		for (int i = 0; i < NN.size(); i++) {
			for (int j = i + 1; j < NN.size(); j++) {
				MBR MergeMBR;
				GetMerge(MergeMBR, NN[i].MBR, NN[j].MBR);
				double Area = GetArea(MergeMBR);
				if (Area > MaxArea) {
					MaxArea = Area; 
					temppair.first = i; temppair.second = j;
				}
			}
		}
		return temppair;
	}
	
	void Partition(int Op1, int Op2) {
		N1.clear(); N2.clear();
		N1.push_back(NN[Op1]); MBR M1 = N1[0].MBR;
		N2.push_back(NN[Op2]); MBR M2 = N2[0].MBR;
		int restEntry = NN.size() - 2;
		for (int i = 0; i < NN.size(); i++) {
			if ((i == Op1) || (i == Op2)) continue;
			if (restEntry <= (Capacity >> 1 - N1.size())) {
				N1.push_back(NN[i]);
			}
			else {
				if (restEntry <= (Capacity >> 1 - N2.size())) {
					N2.push_back(NN[i]);
				}
				else {
					MBR MM1; 
					GetMerge(MM1, M1, NN[i].MBR);
					double A1 = GetArea(MM1) - GetArea(M1);
					MBR MM2;
					GetMerge(MM2, M2, NN[i].MBR);
					double A2 = GetArea(MM2) - GetArea(M2);
					if (A1 < A2) {
						N1.push_back(NN[i]);
						M1 = MM1;
					}
					else {
						N2.push_back(NN[i]);
						M2 = MM2;
					}
				}
			}

			restEntry--;
		}
		return;
	}
	
	void Split(int N, Entry& CurEntry) {
		NN.clear();
		for (int i = 0; i < Tree[N].entry_num; i++) {
			NN.push_back(Tree[N].entries[i]);
		}
		NN.push_back(CurEntry);
		pair<int, int> Op = Promote(N);
		Partition(Op.first, Op.second);

		// update CurNode, store N1
		Tree[N].entry_num = N1.size();
		for (int i = 0; i < N1.size(); i++) {
			Tree[N].entries[i] = N1[i];
		}
		UpdateMBR(N);// update MBR!!!

		// Create a new Node to store N2
		CreateNode(Tree[N].PID, N2.size(), Tree[N].leaf, false); // NewNode cant be root
		int NewNode = Tree.back().ID;
		for (int i = 0; i < N2.size(); i++) {
			Tree[NewNode].entries[i] = N2[i];
			if (!Tree[NewNode].leaf) {
				Tree[Tree[NewNode].entries[i].SonNode].PID = Tree[NewNode].ID;// update son_node's parent_node!!!!
			}
		}

		// Create a new Entry2 to be the NewNode's parent Entry
		// update Entry2's MBR

		MBR NewMBR;
		ClearMBR(NewMBR);
		for (int j = 0; j < Tree[NewNode].entry_num; j++) {
			for (int k = 0; k < dimension; k++) {
				NewMBR.L[k] = min(NewMBR.L[k], Tree[NewNode].entries[j].MBR.L[k]);
				NewMBR.U[k] = max(NewMBR.U[k], Tree[NewNode].entries[j].MBR.U[k]);
			}
		}
		Entry Entry2;
		CreateEntry(Entry2, -1, Tree[NewNode].ID, NewMBR); // Entry in non-leaf node has no TID
		
		if (Tree[N].root) {
			CreateNode(-1, 2, false, true);
			root = Tree.size() - 1; Tree[N].root = false; // update root status
			Entry Entry1; // Create a new Entry1 in root
			MBR TempMBR; // Temp value
			ClearMBR(TempMBR);
			CreateEntry(Entry1, -1, Tree[N].ID, TempMBR);
			Tree[root].entries[0] = Entry1; Tree[N].PID = Tree[root].ID; UpdateMBR(N);
			Tree[root].entries[1] = Entry2; Tree[NewNode].PID = Tree[root].ID; //UpdateMBR(NewNode);
		}
		else {
			int PID = Tree[N].PID;
			if (Tree[PID].entry_num < Capacity) {
				Tree[PID].entries[Tree[PID].entry_num++] = Entry2; //UpdateMBR(NewNode);
			}
			else Split(PID, Entry2);
		}

		return;
	}
	
	void Insert(int N, Entry& CurEntry) {
		if (!Tree[N].leaf) {
			int MinEntry = -1; double MinArea = INFINITY;
			for (int i = 0; i < Tree[N].entry_num; i++) {
				MBR MergeMBR;
				GetMerge(MergeMBR, Tree[N].entries[i].MBR, CurEntry.MBR);
				double IncArea = GetArea(MergeMBR) - GetArea(Tree[N].entries[i].MBR);
				if (IncArea < MinArea) {
					MinArea = IncArea;
					MinEntry = i;
				}
			}

			Insert(Tree[N].entries[MinEntry].SonNode, CurEntry);
			UpdateMBR(N);
		}
		else {
			if (Tree[N].entry_num < Capacity) {
				Tree[N].entries[Tree[N].entry_num++] = CurEntry;
				UpdateMBR(N);
			}
			else Split(N, CurEntry);
		}
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
	
	void AddAllEntry(Node& CurNode) {
		if (!CurNode.leaf) {
			for (int i = 0; i < CurNode.entry_num; i++) {
				AddAllEntry(Tree[CurNode.entries[i].SonNode]);
			}
		}
		else {
			for (int i = 0; i < CurNode.entry_num; i++) {
				Candidate.push_back(CurNode.entries[i].TID);
			}
		}
	}

	void RangeQueryMemory(Node& CurNode, Entry& Q) {
		if (!CurNode.leaf) {
			for (int i = 0; i < CurNode.entry_num; i++) {
				if (EmbodyMBR(Q.MBR,CurNode.entries[i].MBR)) { // Q include CurNode.entries[i]
					AddAllEntry(Tree[CurNode.entries[i].SonNode]);
				}
				else {
					if (GetIntersection(CurNode.entries[i].MBR, Q.MBR)) {
						RangeQueryMemory(Tree[CurNode.entries[i].SonNode], Q);
					}
				}
			}
		}
		else {
			for (int i = 0; i < CurNode.entry_num; i++) {
				if (GetIntersection(CurNode.entries[i].MBR, Q.MBR)) {
					Candidate.push_back(CurNode.entries[i].TID);
				}
			}
		}
		return;
	}
	
	void Build(const char* filename, int Datasize) { // a building example //char* filename
		Tree.clear();
		MBR EmptyMBR;
		for (int i = 0; i < dimension; i++) EmptyMBR.L[i] = EmptyMBR.U[i] = 0;
		CreateNode(-1, 0, true, true); // create root;
		root = Tree.back().ID;
		for (int i = 0; i < Datasize; i++) {
			MBR MBR;
			CreateMBR(MBR, i);
			Entry NewEntry;
			CreateEntry(NewEntry, i, -1, MBR);
			Insert(root, NewEntry);
		}

		WriteToDisk(filename);
		cout << "Rtree root: " << root << endl;
		cout << "Rtree size: " << Tree.size() << endl;
	}
}
#endif
