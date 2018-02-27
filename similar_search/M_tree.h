#ifndef _MTREE_H
#define _MTREE_H

#include <iostream>     
#include <algorithm>
#include <DFD.h>
#include <Init.h>
#include <vector>

using namespace std;

// better struct
// each node store a main object.Tid
// *****use the best first search to accelerate DFD*****

struct Object_mtree {
	int Tid;
	double dis_p, radius;
	int son_node;
};

struct node_mtree {
	vector<Object_mtree> objects;
	bool leaf, root;
	int ID, parent_object, parent_node;	// ID in M_tree
};

const int Capacity_mtree = 20;
int nodenum_mtree;	//number of mtree nodes
int root_mtree;
vector<node_mtree> M_tree;
vector<Object_mtree> N1, N2, NN;
vector<int> results, N_in;
Object_mtree O_n;
node_mtree N_temp;

double Get_distance(int a, int b) {
	return double_DFD(&(All_Data[a].Points), &(All_Data[b].Points));
}

double Getdisp_mtree(int Oid, int Np, int Op) {
	if ((Np == -1) || (Op == -1)) return 0;
	return Get_distance(Oid, M_tree[Np].objects[Op].Tid);
}

void CreateRoot_mtree() {
	M_tree.clear();	//	clear the mtree
	node_mtree tempnode;
	tempnode.objects.clear();
	tempnode.parent_node = -1;
	tempnode.parent_object = -1;
	tempnode.leaf = true;
	tempnode.root = true;
	tempnode.ID = 0;
	nodenum_mtree = 0; root_mtree = 0;
	M_tree.push_back(tempnode);
	return;
}

void CreateNewNode_mtree(node_mtree* N, int ID, bool leaf, bool root, int parent_object, int parent_node) {
	(*N).objects.clear();
	(*N).parent_node = parent_node;
	(*N).parent_object = parent_object;
	(*N).leaf = leaf;
	(*N).root = root;
	(*N).ID = ID;
	return;
}

void RangeQuery_mtree(int N, int Q_Tid, double search_r) {
	int Np = M_tree[N].parent_node;
	int Op = M_tree[N].parent_object;
	double dpq;
	if (!M_tree[N].root) {
		dpq= Get_distance(M_tree[Np].objects[Op].Tid, Q_Tid);
	} 
	else dpq = 0;

	if (!M_tree[N].leaf) {
		for (int i = 0; i < M_tree[N].objects.size(); i++) {
			double dpr = M_tree[N].objects[i].dis_p;
			if ((abs(dpq - dpr) < search_r + M_tree[N].objects[i].radius) || (M_tree[N].root)) {	// root has no parent, so can not satisfy triangle inequality
				double dqr = Get_distance(M_tree[N].objects[i].Tid, Q_Tid);
				if (dqr < search_r + M_tree[N].objects[i].radius) {
					int son_id = M_tree[N].objects[i].son_node;
					RangeQuery_mtree(son_id, Q_Tid, search_r);
				}
			}
		}
	}
	else {
		for (int i = 0; i < M_tree[N].objects.size(); i++) { 
			double dpr = M_tree[N].objects[i].dis_p;
			if ((abs(dpq - dpr) < search_r) || (M_tree[N].root)) { // root has no parent, so can not satisfy triangle inequality
				double dqr = Get_distance(M_tree[N].objects[i].Tid, Q_Tid);
				if (dqr < search_r) results.push_back(M_tree[N].objects[i].Tid);
			}
		}
	}
	return;
}

pair<int, int> Promote_mtree(int N) { // M_LB_DIST
	pair<int, int> temppair;

	if (M_tree[N].root) { // root no parent, select the first two objects
		temppair.first = M_tree[N].objects[0].Tid;
		temppair.second = M_tree[N].objects[1].Tid;
		return temppair;
	}

	double Maxd = 0; int Maxi = 0;
	for (int i = 0; i < NN.size(); i++) {
		if (NN[i].dis_p >= Maxd) {
			Maxd = NN[i].dis_p; Maxi = i;
		}
	}

	temppair.first = M_tree[M_tree[N].parent_node].objects[M_tree[N].parent_object].Tid;
	temppair.second = NN[Maxi].Tid;
	return temppair;
}

void Partition_mtree(int N, int O1, int O2) { // O1 is parent_object.Tid, O2 is promoted Tid
	N1.clear(); N2.clear();
	for (int i = 0; i < NN.size(); i++) {
		int Oid = NN[i].Tid;
		double d1 = Get_distance(O1, Oid);
		double d2 = Get_distance(O2, Oid);
		// then update son_node's parent_objects (if son_node is exist)
		if (d1 < d2) {	// choose nearest object, N1.size <= capacity
			NN[i].dis_p = d1;
			N1.push_back(NN[i]);
			if (NN[i].son_node>=0) M_tree[NN[i].son_node].parent_object = N1.size() - 1;
		}
		else {
			NN[i].dis_p = d2;
			N2.push_back(NN[i]);
			if (NN[i].son_node>=0) M_tree[NN[i].son_node].parent_object = N2.size() - 1;
		}
	}	
}

void ClearNode_mtree(node_mtree* N) {
	(*N).objects.clear();
	(*N).parent_node = NULL;
	(*N).parent_object = NULL;
	(*N).leaf = false;
	(*N).root = false;
	(*N).ID = 0;
	return;
}

Object_mtree CreateNewObject_mtree(int Tid, double radius, double dis_p, int son_node) {
	Object_mtree O;
	O.Tid = Tid;
	O.radius = radius;
	O.dis_p = dis_p;
	O.son_node = son_node;
	return O;
}

void UpdateR_mtree(int N,int O) {//M_tree[N].objects[O] need to be updated. 
	M_tree[N].objects[O].radius = 0;
	int N_son = M_tree[N].objects[O].son_node;
	int pid = M_tree[N].objects[O].Tid;
	if (M_tree[N_son].leaf) {
		for (int i = 0; i < M_tree[N_son].objects.size(); i++) {
			int Oid = M_tree[N_son].objects[i].Tid;
			//M_tree[N_son].objects[i].dis_p = Get_distance(Oid, pid);
			M_tree[N].objects[O].radius = max(M_tree[N].objects[O].radius, M_tree[N_son].objects[i].dis_p);
		}
	}
	else {
		for (int i = 0; i < M_tree[N_son].objects.size(); i++) {
			int Oid = M_tree[N_son].objects[i].Tid;
			//M_tree[N_son].objects[i].dis_p = Get_distance(Oid, pid);
			M_tree[N].objects[O].radius = max(M_tree[N].objects[O].radius, M_tree[N_son].objects[i].dis_p + M_tree[N_son].objects[i].radius);
		}
	}
	return;
}

void split_mtree(int N, Object_mtree O_n) {
	NN.clear(); NN.push_back(O_n);	// add new entry O_n to the N's objects
	for (int i = 0; i < M_tree[N].objects.size(); i++) {
		NN.push_back(M_tree[N].objects[i]);
	}

	pair<int, int> O = Promote_mtree(N);	// Promote two entries(Trajectory ID) from NN
	Partition_mtree(N, O.first, O.second);	// Partition NN into 2 groups, N1 and N2
	CreateNewNode_mtree(&N_temp, -1, M_tree[N].leaf, M_tree[N].root, -1, M_tree[N].parent_node);	// Create a new node N_temp, N and N_temp are brothers in mtree
	M_tree[N].objects = N1;
	N_temp.objects = N2;	// save N1 into N.objects, save N2 into N_temp.objects
	nodenum_mtree++; N_temp.ID = nodenum_mtree; M_tree.push_back(N_temp);   int N_t = nodenum_mtree;	// add N_temp into mtree as a new node, N_t is the pointer
	//if not leaf, update objects' son_node in N2 and N1
	if (!M_tree[N_t].leaf) {
		for (int i = 0; i < M_tree[N_t].objects.size(); i++) {
			int son_id = M_tree[N_t].objects[i].son_node;
			M_tree[son_id].parent_node = N_t;
		}
	}
	if (!M_tree[N].leaf) {
		for (int i = 0; i < M_tree[N].objects.size(); i++) {
			int son_id = M_tree[N].objects[i].son_node;
			M_tree[son_id].parent_node = N;
		}
	}
	Object_mtree O_temp = CreateNewObject_mtree(O.second, 0, 0, N_t);	// Create a new object O_temp, N_t is the subtree
	
	// add O_temp to the parent node
	if (M_tree[N].root) {	// split root of mtree
		M_tree[N].root = false; M_tree[N_t].root = false;	// modify .root 
		// create a new root node
		node_mtree newroot;
		CreateNewNode_mtree(&newroot, -1, false, true, -1, -1);
		nodenum_mtree++; root_mtree = nodenum_mtree; newroot.ID = nodenum_mtree;
		M_tree.push_back(newroot); // add newroot

		// (*N).parent_object is NULL, so create a new object
		Object_mtree O_1 = CreateNewObject_mtree(O.first, 0, 0, M_tree[N].ID);
		// add O_1 and O_temp into new root, update parent_object and parent_node
		M_tree[root_mtree].objects.push_back(O_1);
		M_tree[N].parent_object = M_tree[root_mtree].objects.size() - 1;
		M_tree[N].parent_node = root_mtree;
		M_tree[root_mtree].objects.push_back(O_temp);
		M_tree[N_t].parent_object = M_tree[root_mtree].objects.size() - 1;
		M_tree[N_t].parent_node = root_mtree;
	}
	else {
		int Np = M_tree[N].parent_node;
		if (M_tree[Np].objects.size() < Capacity_mtree) {
			O_temp.dis_p = Getdisp_mtree(O_temp.Tid, Np, M_tree[N].parent_object);
			M_tree[Np].objects.push_back(O_temp);
			M_tree[N_t].parent_object = M_tree[Np].objects.size() - 1;
			M_tree[N_t].parent_node = Np;
		}
		else split_mtree(Np, O_temp);
	}

	// N.parent_object as O1, update (O1,O2).raidus and dis_p in (N1,N2);
	int Np = M_tree[N].parent_node;
	UpdateR_mtree(Np, M_tree[N].parent_object);
	UpdateR_mtree(Np, M_tree[N_t].parent_object);
	return;
}

void GetNin_mtree(int N, Object_mtree* O_n) {
	N_in.clear();
	for (int i = 0; i < M_tree[N].objects.size(); i++) { // Get N_in :dfd(Tra_id,Tid_r)< radius_r
		int Tid = M_tree[N].objects[i].Tid;
		if (Get_distance(Tid, (*O_n).Tid) < M_tree[N].objects[i].radius) N_in.push_back(i);
	}
	return;
}

int GetMin_fromNin_mtree(int N, Object_mtree* O_n) {
	double Mind = INFINITY; int Mini = -1;
	for (int i = 0; i < N_in.size(); i++) {
		int Oid = M_tree[N].objects[N_in[i]].Tid;
		double dis_temp = Get_distance(Oid, (*O_n).Tid);
		if (dis_temp < Mind) {
			Mind = dis_temp; Mini = N_in[i];
		}
	}
	return Mini;
}

int GetMin_fromN_mtree(int N, Object_mtree* O_n) {
	double Mind = INFINITY; int Mini = -1;
	for (int i = 0; i < M_tree[N].objects.size(); i++) {
		int Oid = M_tree[N].objects[i].Tid;
		double dis_temp = Get_distance(Oid, (*O_n).Tid) - M_tree[N].objects[i].radius;
		if (dis_temp < Mind) {
			Mind = dis_temp; Mini = i;
		}
	}
	// update Oid.radius
	M_tree[N].objects[Mini].radius = Get_distance(M_tree[N].objects[Mini].Tid, (*O_n).Tid);
	return Mini;
}

void insert_mtree(int N, Object_mtree* O_n) {
	if (!M_tree[N].leaf) { // Mtree_id is not a leaf node store routing objects	
		GetNin_mtree(N, O_n);
		int O_r = -1;
		if (!N_in.empty()) {
			O_r = GetMin_fromNin_mtree(N, O_n);// choose node Oid, min(d(Tra_id,Tid_r))
		}
		else {
			O_r = GetMin_fromN_mtree(N, O_n);  // N_in is empty, choose node Oid, min(d(Tra_id,Tid_r)-radius_r)
		}
		insert_mtree(M_tree[N].objects[O_r].son_node, O_n);
	}
	else { // node_is is a leaf node, store entries
		if (M_tree[N].objects.size() < Capacity_mtree) { // not full
			(*O_n).dis_p = Getdisp_mtree((*O_n).Tid, M_tree[N].parent_node, M_tree[N].parent_object);//update O_n.dis_p
			M_tree[N].objects.push_back(*O_n);
			if (!M_tree[N].root) {
				int Np = M_tree[N].parent_node;
				UpdateR_mtree(Np, M_tree[N].parent_object);
			}
		}
		else { // if this node larger than capacity, then spilt it 
			split_mtree(N,*O_n);
		}
	}
	return;
}

void build_mtree() {
	CreateRoot_mtree();	// create the root node
	for (int i = 0; i < All_Data.size(); i++) {
	//for (int i = 0; i < 30; i++) {
		O_n.Tid = i; O_n.dis_p = 0; O_n.radius = 0; O_n.son_node = -1;
		cout << i << endl;
		//if (i == 10) system("pause");
		insert_mtree(root_mtree,&O_n);
	}
	cout << "M_tree build finish!" << endl;
	return;
}
#endif