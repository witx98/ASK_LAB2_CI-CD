#include <iostream>

using namespace std;

//H - wysokoœæ siatki, W - szerokoœæ siatki, nW - liczba elementów na szerokoœci, 
//nN - liczba wierzcho³ków, nE - liczba elementów (12 kwadratów)
//nH - liczba elementów wysokoœci

struct GlobalData
{
	double H, W;
	unsigned int nH, nW;
	double deltaX;
	double deltaY;
	int nN;
	int nE;
	double k;
	double cp;
	double Ro;
	double alfa;
	double t_alfa;
	double t0;
	double time;
	double timeStep;
	double numberOfIterations;
	GlobalData(double h, double w, unsigned int nh, unsigned int nw, double k, double cp, double Ro, double alfa, double t_alfa, double t0, double time, double timeStep) {
		this->H = h;
		this->W = w;
		this->nH = nh;
		this->nW = nw;

		this->k = k;
		this->cp = cp;
		this->Ro = Ro;
		this->alfa = alfa;
		this->t_alfa = t_alfa;
		this->t0 = t0;
		this->time = time;
		this->timeStep = timeStep;

		this->deltaX = this->W / (this->nW - 1);
		this->deltaY = this->H / (this->nH - 1);
		this->nN = this->nH * this->nW;  //liczba wêz³ów
		this->nE = (this->nH - 1) * (this->nW - 1); //liczba elementów
		this->numberOfIterations = time / timeStep;
	}
	GlobalData() {
		this->H = 0;
		this->W = 0;
		this->nH = 0;
		this->nW = 0;
		this->deltaX = 0;
		this->deltaY = 0;
		this->nN = 0;  //liczba wêz³ów
		this->nE = 0; //liczba elementów

		this->k = 0;
		this->cp = 0;
		this->Ro = 0;
		this->alfa = 0;
		this->t_alfa = 0;
		this->t0 = 0;
		this->time = 0;
		this->timeStep = 0;
		this->numberOfIterations = 0;
	}
};

struct Node
{ //wspó³rzêdne wierzcho³ka
	int nodeNumber;
	double x;
	double y;
	double t0;
	bool Bwc = 0;
};

struct Element
{
	Node ID[4];
	double** Hl;
	double** Cl;
	double** Hbc;
	double** HLBc;
	double* Pl;
};

class FEM_GRID {
public:
	GlobalData data;
	Node* nodesTable;
	Element* elementsTable;

public:
	FEM_GRID(GlobalData data) {
		this->data = data;
		nodesTable = new Node[this->data.nN];
		elementsTable = new Element[this->data.nE];
	}
private:
	void fillNodeTable() {
		int nodeIterator = 0;
		//wype³nianie tablicy node'ow
		for (int i = 0; i < data.nW; i++)
		{
			for (int j = 0; j < data.nH; j++)
			{ //przypisanie wspó³rzêdnych do wierzcho³ków 
				nodesTable[nodeIterator].x = data.deltaX * i;
				nodesTable[nodeIterator].y = data.deltaY * j;
				nodesTable[nodeIterator].nodeNumber = nodeIterator;
				nodesTable[nodeIterator].t0 = data.t0;

				if (nodesTable[nodeIterator].x == 0 || nodesTable[nodeIterator].y == 0 || nodesTable[nodeIterator].x == data.W || nodesTable[nodeIterator].y == data.H) {
					nodesTable[nodeIterator].Bwc = 1;
				}
				nodeIterator++;
			}

		}
	}
	void fillElementsTable() {
		int nodeCounter = 0;
		for (size_t i = 0; i < data.nE; i++)
		{

			elementsTable[i].ID[0] = nodesTable[nodeCounter];
			elementsTable[i].ID[1] = nodesTable[nodeCounter + data.nH];
			elementsTable[i].ID[2] = nodesTable[nodeCounter + data.nH + 1];
			elementsTable[i].ID[3] = nodesTable[nodeCounter + 1];

			nodeCounter++;
			if ((nodeCounter + 1) % data.nH == 0) {
				nodeCounter++;
			}
		}
	}
public:
	void createMesh() {
		fillNodeTable();
		fillElementsTable();
	}
	void showElement(int elementNumber) {
		cout << "Element " << elementNumber << " sklada sie z:" << endl;

		for (size_t i = 0; i < 4; i++)
		{
			cout << "ID" << i << ": " << elementsTable[elementNumber].ID[i].nodeNumber << endl;
			cout << "\tx: " << elementsTable[elementNumber].ID[i].x << endl;
			cout << "\ty: " << elementsTable[elementNumber].ID[i].y << endl;
			cout << "\tWarunek brzegowy: " << elementsTable[elementNumber].ID[i].Bwc << endl;
			cout << "\tTemperatura pocz¹tkowa: " << elementsTable[elementNumber].ID[i].t0<< endl;
		}
	}
};

struct SOE {
	double** Hglobal;
	double** Cglobal;
	double** HLBcglobal;
	double* Pglobal;
	double** HFinal;
	double* PFinal;
};