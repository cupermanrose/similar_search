#ifndef _MTREE_H
#define _MTREE_H

#include <iostream> 
#include <algorithm>
#include <fstream>
#include <vector>
#include <BLtoDist.h>

using namespace std;

namespace Mtree {
	const int Capacity = 20;
	
	struct Node;

	struct Entry {
		int Tid,son_node; // Tid: feature value of the object, son_node: offset from the start address
		double dis_p, radius; // dis_p: distance with parent object, radius: cover radius
	};

	struct Node {
		Entry entries[Capacity]; 
		int entry_num, center, ID, PID; // entry number in the node, center is the center Entry, ID: the tree ID of the node, PID: parent ID
		bool leaf, root; // leaf and root label
	};

	vector<Node> Tree;
	vector<Entry> NN, N1, N2;
	vector<int> Nin, RQans;
	int root;

	double GetDistance(Entry& A, Entry& B) { // distance
		double dfd = double_DFD(All_Query[A.Tid].Points, All_Data[B.Tid].Points);
		return dfd;
	}

	void AddAllEntry(Node& CurNode) {
		if (!CurNode.leaf) {
			for (int i = 0; i < CurNode.entry_num; i++) {
				AddAllEntry(Tree[CurNode.entries[i].son_node]);
			}
		}
		else {
			for (int i = 0; i < CurNode.entry_num; i++) {
				RQans.push_back(CurNode.entries[i].Tid);
			}
		}
	}

	void RangeQueryMemory(Node& CurNode, Entry& Q, double DisPQ) {
		/*double DisPQ;
		if (!CurNode.root) {
			DisPQ = GetDistance(Q, CurNode.entries[CurNode.center]);
		}
		else DisPQ = 0;*/

		if (!CurNode.leaf) {
			for (int i = 0; i < CurNode.entry_num; i++) {
				if (abs(DisPQ - CurNode.entries[i].dis_p) <= (Q.radius + CurNode.entries[i].radius)) {
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
					double DisRQ = GetDistance(Q, CurNode.entries[i]);
					if (DisRQ <= Q.radius) RQans.push_back(CurNode.entries[i].Tid);
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

	void CreateNode(int entry_num, int PID, int center, bool leaf, bool root) {
		Node* NewNode = new Node;
		NewNode->entry_num = entry_num;
		NewNode->PID = PID;
		NewNode->center = center;
		NewNode->leaf = leaf;
		NewNode->root = root;
		NewNode->ID = Tree.size();
		Tree.push_back(*NewNode);
		delete NewNode;
		return;
	}

	pair<int, int> Promote(int N) {
		pair<int, int> temppair;
		if (Tree[N].root) { // root doesn't have center, select the first two entries;
			temppair.first = 0; temppair.second = 1;
			return temppair;
		}
		double MaxDis = -1;
		temppair.first = Tree[N].center;
		// Find the farthest Entry with center be the Op2
		for (int i = 0; i < NN.size(); i++) {
			if (i == Tree[N].center) continue;
			if (NN[i].dis_p >= MaxDis) {
				MaxDis = NN[i].dis_p;
				temppair.second = i;
			}
		}
		return temppair;
	}

	pair<int, int> Partition(int Op1, int Op2) {
		pair<int, int> temppair; // store new center entries of two node 
		N1.clear(); N2.clear();
		for (int i = 0; i < NN.size(); i++) {
			double d1 = GetDistance(NN[i], NN[Op1]);
			double d2 = GetDistance(NN[i], NN[Op2]);
			if (((d1 < d2) || (i == Op1)) && (i != Op2)) { // NN[i] is closer to Op1, op2 cannot be in the N1
				NN[i].dis_p = d1; // update dis_p
				N1.push_back(NN[i]);
				if (i == Op1) temppair.first = N1.size() - 1;
			}
			else {
				NN[i].dis_p = d2;
				N2.push_back(NN[i]);
				if (i == Op2) temppair.second = N2.size() - 1;
			}
		}
		return temppair;
	}

	void UpdateRadius(int N) {
		if (Tree[N].root) return; // root need not to update radius
		int PID = Tree[N].PID;
		for (int i = 0; i < Tree[PID].entry_num; i++) {
			if (Tree[N].ID == Tree[PID].entries[i].son_node) {
				Tree[PID].entries[i].radius = 0;
				if (Tree[N].leaf) {
					for (int j = 0; j < Tree[N].entry_num; j++) {
						Tree[PID].entries[i].radius = max(Tree[PID].entries[i].radius, Tree[N].entries[j].dis_p);
					}
				}
				else {
					for (int j = 0; j < Tree[N].entry_num; j++) {
						Tree[PID].entries[i].radius = max(Tree[PID].entries[i].radius, Tree[N].entries[j].dis_p + Tree[N].entries[j].radius);
					}
				}
				return;
			}
		}
	}

	void Split(int N, Entry& CurEntry) {
		NN.clear(); 
		// let NN= CurNode.entries + CurEntry; 
		for (int i = 0; i < Tree[N].entry_num; i++) {
			NN.push_back(Tree[N].entries[i]);
		}
		NN.push_back(CurEntry); // CurEntry is the last element
		
		// let CurNode's center entry always be the Op.first, and promote Op.second;
		pair<int,int> Op = Promote(N);
		Op = Partition(Op.first, Op.second); // Op is the new center of two nodes
		
		// update CurNode, store N1
		Tree[N].entry_num = N1.size();
		Tree[N].center = Op.first;
		for (int i = 0; i < N1.size(); i++) {
			Tree[N].entries[i] = N1[i];
		}
		UpdateRadius(N);// update radius!!!

		// Create a new Node to store N2
		CreateNode(N2.size(), Tree[N].PID, Op.second, Tree[N].leaf, false); // NewNode cant be root
		int NewNode = Tree.back().ID;
		for (int i = 0; i < N2.size(); i++) {
			Tree[NewNode].entries[i] = N2[i];
			if (!Tree[NewNode].leaf) {
				Tree[Tree[NewNode].entries[i].son_node].PID = Tree[NewNode].ID;// update son_node's parent_node!!!!
			}
		}
		
		// Create a new Entry2 to be the NewNode's parent Entry
		Entry Entry2;
		CreateEntry(Entry2, N2[Op.second].Tid, Tree[NewNode].ID, 0, 0);
		// update Entry2's radius
		for (int i = 0; i < Tree[NewNode].entry_num; i++) {
			if (Tree[NewNode].leaf) {
				Entry2.radius = max(Entry2.radius, Tree[NewNode].entries[i].dis_p);
			}
			else Entry2.radius = max(Entry2.radius, Tree[NewNode].entries[i].dis_p + Tree[NewNode].entries[i].radius);
		}

		if (Tree[N].root) {
			CreateNode(2, -1, -1, false, true);
			root = Tree.size() - 1; Tree[N].root = false; // update root status
			Entry Entry1; // Create a new Entry1 && Entries in root don't have dis_p
			CreateEntry(Entry1, N1[Op.first].Tid, Tree[N].ID, 0, 0);
			Tree[root].entries[0] = Entry1; Tree[N].PID = Tree[root].ID; UpdateRadius(N);
			Tree[root].entries[1] = Entry2; Tree[NewNode].PID = Tree[root].ID; UpdateRadius(NewNode);
		}
		else {
			int PID = Tree[N].PID;
			if (!Tree[PID].root) Entry2.dis_p = GetDistance(Tree[PID].entries[Tree[PID].center], Entry2);
			if (Tree[PID].entry_num < Capacity) {
				Tree[PID].entries[Tree[PID].entry_num++] = Entry2; UpdateRadius(NewNode);
			}
			else Split(PID, Entry2);
		}

		return;
	}

	void Insert(int N, Entry& CurEntry) {
		if (!Tree[N].leaf) {
			double* Dis = new double[Capacity]; // get distances with all Or and On
			for (int i = 0; i < Tree[N].entry_num; i++) {
				Dis[i] = GetDistance(Tree[N].entries[i], CurEntry);
			}
			// all element Or in Nin: d(Or,On)<=Or.radius
			Nin.clear();
			for (int i = 0; i < Tree[N].entry_num; i++) {
				int OrID = Tree[N].entries[i].Tid;
				if (Dis[i] < Tree[N].entries[i].radius) Nin.push_back(i);
			}

			int MinEntry = -1; double MinDis = INFINITE;
			if (!Nin.empty()) { // d(Or, On) is minimum
				for (int i = 0; i < Nin.size(); i++) {
					if (MinDis >= Dis[Nin[i]]) {
						MinEntry = Nin[i];
						MinDis = Dis[Nin[i]];
					}
				}
			}
			else { // d(Or,On)-r(Or) is minimum
				for (int i = 0; i < Tree[N].entry_num; i++) {
					if (MinDis >= Dis[i]) {
						MinEntry = i;
						MinDis = Dis[i];
					}
				}
			}

			CurEntry.dis_p = Dis[MinEntry]; // update dis_p
			delete[] Dis;
			Insert(Tree[N].entries[MinEntry].son_node, CurEntry);
			UpdateRadius(N);
		}
		else {
			if (Tree[N].entry_num < Capacity) {
				Tree[N].entries[Tree[N].entry_num++] = CurEntry;
				UpdateRadius(N);
			}
			else Split(N, CurEntry);
		}
	}

	void WriteToDisk(const char* filename) {
		fstream MTfile(filename, ios::out|ios::binary);
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
		CreateNode(0, -1, -1, true, true); // create root;
		root = Tree.back().ID;
		for (int i = 0; i < Datasize; i++) {
			Entry NewEntry;
			CreateEntry(NewEntry, i, -1, 0, 0);
			Insert(root, NewEntry);
		}
		WriteToDisk(filename);
		cout << "Mtree root: " << root << endl;
		cout << "Mtree size: " << Tree.size() << endl;
	}

}


#endif