#include <iostream>
#include<fstream>
#include <string>
#include <math.h>
#include <float.h>
#include "Structures.h"
#include "Elem4.h"
#include "Functions.h"
using namespace std;
 
    
double** H_fun(double factor, Elem4 * e4, double x[],double y[]) {
    //Macierze [ipc*ipc]
    //detJ

    //Macierze [size][ipc*ipc]
    //jacobian
    //reverseJacobian
    //dNdX
    //dNdY

    //Macierze [ipc*ipc][size]
    //transposed_dNdX
    //transposed_dNdY

    //Macierze [ipc*ipc][size][size]
    //dNdX_il
    //dNdY_il
    //multiplications

    //Macierze [size][size]
    //H

    double** jacobian = fill_Jacobian_Table(e4->dNdEta, e4->dNdKsi,e4->ipc, e4->size, x, y);
    double* detJ = fill_DetJ_Table(jacobian, e4->ipc);
    double** reverseJacobian = fill_Reverse_Jacobian_Table(jacobian, detJ, e4->size, e4->ipc);
    double** dNdX = fill_dNdX(reverseJacobian, e4->dNdEta, e4->dNdKsi, e4->size, e4->ipc);
    double** dNdY = fill_dNdY(reverseJacobian, e4->dNdEta, e4->dNdKsi, e4->size, e4->ipc);
    double** transposed_dNdX = transpose_matrix(e4->size, e4->ipc, dNdX);
    double** transposed_dNdY = transpose_matrix(e4->size, e4->ipc, dNdY);
    double*** dNdX_il = iloczyn(transposed_dNdX, dNdX, e4->size, e4->ipc);
    double*** dNdY_il = iloczyn(transposed_dNdY, dNdY, e4->size, e4->ipc);
    double*** multiplications = count_multiplications(factor, dNdX_il,  dNdY_il,  detJ, e4->size, e4->ipc);
    double** H = count_H(multiplications, e4->size, e4->ipc, e4->weights);

    //e4->show_Ksi_Eta_Tables();
    //show_Jacobian_Table(jacobian, e4->ipc, e4->size);
    //show_DetJ_Table(detJ, e4->ipc);
    //show_Reverse_Jacobian_Table(reverseJacobian, e4->size, e4->ipc);
    //show_dNdX(dNdX, e4->size, e4->ipc);
    //show_dNdY(dNdY, e4->size, e4->ipc);
    //show_transposed_dNdX(transposed_dNdX, e4->size, e4->ipc);
    //show_transposed_dNdY(transposed_dNdY, e4->size, e4->ipc);
    //show_iloczynX(dNdX_il, e4->size, e4->ipc);
    //show_iloczynY(dNdY_il, e4->size, e4->ipc);
    //show_multiplications(multiplications, e4->size, e4->ipc);
    return H;
}


    
double ** count_C(double** N, double** transposed_N, double *DetJ, int size, int ipc, double cp, double Ro, double* weights) {
    double** C = new double* [size];
    for (size_t i = 0; i < size; i++)
    {
        C[i] = new double[size];
        for (size_t j = 0; j < size; j++)
        {
            C[i][j] = 0;
        }
    }
    double* Weight = new double[ipc * ipc];
    int w = 0;
    for (int i = 0; i < ipc; i++) {
        for (int j = 0; j < ipc; j++) {
            Weight[w] = weights[i] * weights[j];
            w++;
        }
    }
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            for (size_t k = 0; k < ipc*ipc; k++)
            {
                C[i][j] += cp*Ro*N[i][k] * transposed_N[k][j]* DetJ[k]*Weight[k];
            }
        }
    }
    return C;
}
double** C_fun(Elem4* e4, double *x, double *y, double cp, double Ro) {
    //e4->show_N_Table();
    double** Transposed_N = transpose_matrix(e4->size, e4->ipc, e4->N);
    /*cout << endl << "Transposed N" << endl;
    for (int i = 0; i < e4->ipc * e4->ipc; i++)
    {
        for (int j = 0; j < e4->size; j++)
        {
            cout << "[" << i << "][" << j << "] = " << Transposed_N[i][j] << " ";
        }
        cout << endl;
    }*/
    //Macierze [size][ipc*ipc]
    
    double** jacobian = fill_Jacobian_Table(e4->dNdEta, e4->dNdKsi, e4->ipc, e4->size, x, y);
    double* detJ = fill_DetJ_Table(jacobian, e4->ipc);
    double** C = count_C(e4->N, Transposed_N, detJ, e4->size, e4->ipc, cp, Ro, e4->weights);
    return C;
}
    
double** Hbc_fun(Elem4* e4, Element singleElement, double alfa) {
    double** Hbc = new double* [e4->size];
    for (size_t i = 0; i < e4->size; i++)
    {
        Hbc[i] = new double[e4->size];
        for (size_t j = 0; j < e4->size; j++)
        {
            Hbc[i][j] = 0;
        }
    }
    double*** partialHbc = new double** [e4->size];
    for (size_t i = 0; i < e4->size; i++)
    {
        partialHbc[i] = new double* [e4->size];
        for (size_t j = 0; j < e4->size; j++)
        {
            partialHbc[i][j] = new double[e4->size];
            for (size_t k = 0; k < e4->size; k++)
            {
                partialHbc[i][j][k] = 0;
            }
        }
    }
    int integrationPointIndex = 0;
    for (size_t i = 0; i < e4->size; i++)
    {
        int secondIndex = i + 1;
        if (secondIndex == 4)
        {
            secondIndex = 0;
        }
        if (singleElement.ID[i].Bwc == 1 && singleElement.ID[secondIndex].Bwc == 1) {
            for (size_t j = 0; j < e4->ipc; j++)
            {
                for (size_t k = 0; k < e4->size; k++)
                {
                    for (size_t l = 0; l < e4->size; l++)
                    {
                        partialHbc[i][k][l] += (e4->NBc[k][integrationPointIndex] * e4->NBc[l][integrationPointIndex]) * e4->weights[j] * alfa * distance(singleElement.ID[i], singleElement.ID[secondIndex]) / 2;
                    }
                }
                integrationPointIndex++;
            }
        }
        else
        {
            integrationPointIndex += e4->ipc;
        }
    }

    for (size_t i = 0; i < e4->size; i++)
    {
        for (size_t j = 0; j < e4->size; j++)
        {
            for (size_t k = 0; k < e4->size; k++)
            {
                Hbc[i][j] += partialHbc[k][i][j];
            }
        }
    }
    return Hbc;
}
double** HLBc_fun(double** Hl, double** Hbc, int size) {
    double** HLBc = new double* [size];
    for (size_t i = 0; i < size; i++)
    {
        HLBc[i] = new double[size];
        for (size_t j = 0; j < size; j++)
        {
            HLBc[i][j] = Hl[i][j] + Hbc[i][j];
        }   
    }
    return HLBc;
}
double* P_fun(Elem4* e4, Element singleElement, double alfa, double t_alfa) {
    double* P = new double[e4->size];
    for (size_t i = 0; i < e4->size; i++)
    {
        P[i] = 0;
    }
    double** P_partial = new double*[e4->size];
    for (size_t i = 0; i < e4->size; i++)
    {
        P_partial[i] = new double[e4->size];
        for (size_t j = 0; j < e4->size; j++)
        {
            P_partial[i][j] = 0;
        }
    }

    int integrationPointIndex = 0;
    for (size_t i = 0; i < e4->size; i++)
    {
        int secondIndex = i + 1;
        if (secondIndex == 4)
        {
            secondIndex = 0;
        }
        if (singleElement.ID[i].Bwc == 1 && singleElement.ID[secondIndex].Bwc == 1) {
            for (size_t j = 0; j < e4->ipc; j++)
            {
                for (size_t k = 0; k < e4->size; k++)
                {
                    P_partial[i][k] += -1*(e4->NBc[k][integrationPointIndex] * e4->weights[j] * alfa * t_alfa * (distance(singleElement.ID[i], singleElement.ID[secondIndex]) / 2));
                }
            integrationPointIndex++;
            }
        }
        else
        {
            integrationPointIndex += e4->ipc;
        }
    }
    for (size_t i = 0; i < e4->size; i++)
    {
        for (size_t j = 0; j < e4->size; j++)
        {
            P[i] += P_partial[j][i];
        }
    }
    return P;
}
void count_local_matrix(FEM_GRID* mesh, Elem4* e4) {
    double* xx = new double[e4->size];
    double* yy = new double[e4->size];
    for (size_t i = 0; i < mesh->data.nE; i++)
    {
        //cout << "Nr" << i << endl;
        for (size_t j = 0; j < e4->size; j++)
        {
            xx[j] = mesh->elementsTable[i].ID[j].x;
            yy[j] = mesh->elementsTable[i].ID[j].y;
        }
        mesh->elementsTable[i].Hl = H_fun(mesh->data.k, e4, xx, yy);
        mesh->elementsTable[i].Cl = C_fun(e4, xx, yy, mesh->data.cp, mesh->data.Ro);
        mesh->elementsTable[i].Hbc = Hbc_fun(e4, mesh->elementsTable[i], mesh->data.alfa);
        mesh->elementsTable[i].HLBc = HLBc_fun(mesh->elementsTable[i].Hl, mesh->elementsTable[i].Hbc,e4->size);
        mesh->elementsTable[i].Pl = P_fun(e4, mesh->elementsTable[i], mesh->data.alfa, mesh->data.t_alfa);
        //show_H(mesh->elementsTable[i].Hl, e4->size);
        //show_H(mesh->elementsTable[i].Hbc, e4->size);
        //show_P(mesh->elementsTable[i].Pl, e4->size);
    }     
}
void agregate_global_matrix(FEM_GRID* mesh, SOE* sol, int size) {
    sol->Hglobal = new double* [mesh->data.nN];
    sol->Cglobal = new double* [mesh->data.nN];
    sol->HLBcglobal = new double* [mesh->data.nN];
    sol->Pglobal = new double[mesh->data.nN];
    sol->PFinal = new double[mesh->data.nN];

    sol->HFinal = new double *[mesh->data.nN];
    for (size_t i = 0; i < mesh->data.nN; i++)
    {
        sol->Hglobal[i] = new double[mesh->data.nN];
        sol->Cglobal[i] = new double[mesh->data.nN];
        sol->HLBcglobal[i] = new double[mesh->data.nN];
        sol->Pglobal[i] = 0;
        sol->PFinal[i] = 0;
        sol->HFinal[i] = new double[mesh->data.nN];

        for (size_t j = 0; j < mesh->data.nN; j++)
        {
            sol->Hglobal[i][j] = 0.0;
            sol->Cglobal[i][j] = 0.0;
            sol->HLBcglobal[i][j] = 0.0;
            sol->HFinal[i][j] = 0.0;
        }
    }

    for (size_t i = 0; i < mesh->data.nE; i++)
    {
        Element temp = mesh->elementsTable[i];
        for (size_t j = 0; j < size; j++)
        {
            sol->Pglobal[temp.ID[j].nodeNumber] += mesh->elementsTable[i].Pl[j];
            for (size_t k = 0; k < size; k++)
            {
                sol->Hglobal[temp.ID[k].nodeNumber][temp.ID[j].nodeNumber] += mesh->elementsTable[i].Hl[j][k];
                sol->Cglobal[temp.ID[k].nodeNumber][temp.ID[j].nodeNumber] += mesh->elementsTable[i].Cl[j][k];
                sol->HLBcglobal[temp.ID[k].nodeNumber][temp.ID[j].nodeNumber] += mesh->elementsTable[i].HLBc[j][k];
            }
        }
    }
}

void FEM(FEM_GRID* mesh, SOE* sol, string outputFileName) {
    cout << "-------------------------UKLAD ROWNAN-------------------------\n";
    double* t1 = new double[mesh->data.nN];
    for (size_t i = 0; i < mesh->data.nN; i++)
    {
        t1[i] = 0;
    }
    ofstream outputFile(outputFileName);
    outputFile << "Time[s]" << ";" << "MinTemp" << ";" << "MaxTemp" << endl;
    for (size_t i = 0; i < mesh->data.numberOfIterations; i++)
    {
        cout << "\n\n-------------------------Iteration " << i << " -------------------------\n";
        for (size_t j = 0; j < mesh->data.nN; j++)
        {
            t1[j] = 0;
            sol->PFinal[j] = sol->Pglobal[j];
            for (size_t k = 0; k < mesh->data.nN; k++)
            {
                sol->HFinal[j][k] = sol->HLBcglobal[j][k] + sol->Cglobal[j][k] / mesh->data.timeStep;
                sol->PFinal[j] += -(sol->Cglobal[j][k] / mesh->data.timeStep) * mesh->nodesTable[k].t0;
            }
            sol->PFinal[j] *= -1;
        }
        //show_Global_Matrix(sol->HFinal, mesh->data.nN);
        //cout << endl;
        //show_P_Global_Matrix(sol->PFinal, mesh->data.nN);

        double** ukladRownan = combineTables(sol->HFinal, sol->PFinal, mesh->data.nN);
        Gauss(ukladRownan, t1, mesh->data.nN, mesh->data.nN + 1);
        //cout << "\n\nTeperatury:\n";
        for (size_t i = 0; i < mesh->data.nN; i++)
        {
            //cout << t1[i] << " ";
            mesh->nodesTable[i].t0 = t1[i];
        }
        cout << "\n\nMaksymalna i minimalna temperatura w iteracji " << i << ": \n";
        double min = findMin(t1, mesh->data.nN);
        double max = findMax(t1, mesh->data.nN);
        cout << "Max: " << max << endl;
        cout << "Min: " << min << endl;
        outputFile << (i + 1) * mesh->data.timeStep << ";" << min << ";" << max << endl;
    }
    outputFile.close();
} 

int main()
{
    //////////WYBOR TEST CASE'U/////////
    //string fileName = "TestCase1.txt";
	string fileName = "TestCase2.txt";
	string linia;
	fstream plik;
	double H, W;
	int nH, nW, fileCounter = 0;
    double k, cp, Ro, alfa, t_alfa, t0, time, timeStep;
	plik.open(fileName, ios::in);		
	if (plik.good() == true)			
	{
		while (!plik.eof())
		{
			getline(plik, linia, ',');		 
			if (fileCounter == 0) {
				double i = atof(linia.c_str());
				H = i;
			}
			if (fileCounter == 1) {
				double i = atof(linia.c_str());
				W = i;
			}
			if (fileCounter == 2) {
				int i = atoi(linia.c_str());
				nH = i;
			}
			if (fileCounter == 3) {
				int i = atoi(linia.c_str());
				nW = i;
			}
            if (fileCounter == 4) {
                double i = atof(linia.c_str());
                k = i;
            }
            if (fileCounter == 5) {
                double i = atof(linia.c_str());
                cp = i;
            }
            if (fileCounter == 6) {
                double i = atof(linia.c_str());
                Ro = i;
            }
            if (fileCounter == 7) {
                double i = atof(linia.c_str());
                alfa = i;
            }
            if (fileCounter == 8) {
                double i = atof(linia.c_str());
                t_alfa = i;
            }
            if (fileCounter == 9) {
                double i = atof(linia.c_str());
                t0 = i;
            }
            if (fileCounter == 10) {
                double i = atof(linia.c_str());
                time = i;
            }
            if (fileCounter == 11) {
                double i = atof(linia.c_str());
                timeStep = i;
            }
			fileCounter++;
		}
		plik.close();
	}
    
	GlobalData data(H, W, nH, nW, k, cp, Ro, alfa, t_alfa, t0, time, timeStep);
	FEM_GRID *FG = new FEM_GRID(data);
	FG->createMesh();
    Elem4* test = new Elem4(4, 3);

    //Liczenbie macierzy lokalnych
    count_local_matrix(FG, test);
 
    //Agregacja macierzy lokalnych
    SOE* solution = new SOE();
    agregate_global_matrix(FG, solution, test->size);
    //show_Global_Matrix(solution->Hglobal, FG->data.nN);
    //show_Global_Matrix(solution->Cglobal, FG->data.nN);
    //show_Global_Matrix(solution->HLBcglobal, FG->data.nN);
    //show_P_Global_Matrix(solution->Pglobal, FG->data.nN);
    cout << FG->data.nN << endl;
    //Rozwiazanie calego ukladu
    FEM(FG, solution, "temperatures.csv");
    
	return 0;
}
