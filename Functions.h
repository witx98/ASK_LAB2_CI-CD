#include <iostream>

using namespace std;

double** fill_Jacobian_Table(double** dNdEta, double** dNdKsi, int ipc, int size, double x[], double y[]) {
    double** Jacobian = new double*[size];
    for (size_t i = 0; i < size; i++)
    {
        Jacobian[i] = new double[ipc * ipc];
        for (size_t j = 0; j < ipc*ipc; j++)
        {
            Jacobian[i][j] = 0;
        }
    }
    for (int i = 0; i < ipc * ipc; i++)
    {
        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                if (j == 0)
                {
                    Jacobian[j][i] += (x[k] * dNdEta[k][i]);
                }
                else if (j == 1)
                {
                    Jacobian[j][i] += (x[k] * dNdKsi[k][i]);
                }
                else if (j == 2)
                {
                    Jacobian[j][i] += (y[k] * dNdEta[k][i]);
                }
                else if (j == 3)
                {
                    Jacobian[j][i] += (y[k] * dNdKsi[k][i]);
                }
            }
        }
    }
    return Jacobian;
}
void show_Jacobian_Table(double** jacobian, int ipc, int size) {
    cout << "\n" << "Jacobian" << endl;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < ipc * ipc; j++)
        {
            cout << "[" << i << "][" << j << "] = " << jacobian[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}
double* fill_DetJ_Table(double** jacobian, int ipc) {
    double* DetJ = new double[ipc * ipc];
    for (int i = 0; i < ipc * ipc; i++)
    {
        DetJ[i] = jacobian[0][i] * jacobian[3][i] - jacobian[1][i] * jacobian[2][i];
    }
    return DetJ;
}
void show_DetJ_Table(double* detJ, int ipc) {
    cout << "\n" << "DetJ" << endl;
    for (int i = 0; i < ipc * ipc; i++)
    {
        cout << "[" << i << "]=" << detJ[i] << " ";
    }
    cout << "\n";
}
double** fill_Reverse_Jacobian_Table(double** jacobian, double* detJ, int size, int ipc) {
   
    double** ReverseJacobian = new double* [size];
    for (size_t i = 0; i < size; i++)
    {
        ReverseJacobian[i] = new double[ipc * ipc];
        for (size_t j = 0; j < ipc * ipc; j++)
        {
            ReverseJacobian[i][j] = 0;
        }
    }

    int iter = 0;
    for (int i = 0; i < ipc * ipc; i++)
    {
        for (int j = (size)-1; j >= 0; j--)
        {
            if (iter % (size - 1) == 0) {
                ReverseJacobian[iter][i] = jacobian[j][i] / detJ[i];
            }
            else {
                ReverseJacobian[iter][i] = -(jacobian[j][i] / detJ[i]);
                if (ReverseJacobian[iter][i] == -0.0) {
                    ReverseJacobian[iter][i] = 0.0;
                }
            }
            iter++;
        }
        iter = 0;
    }
    return ReverseJacobian;
}
void show_Reverse_Jacobian_Table(double** reverseJacobian, int size, int ipc) {
    cout << "\n" << "Reverse Jacobian" << endl;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < ipc * ipc; j++)
        {
            cout << "[" << i << "][" << j << "] = " << reverseJacobian[i][j] << " ";
        }
        cout << endl;
    }
    cout << "\n";
}
double** fill_dNdX(double** reverseJacobian, double** dNdEta, double** dNdKsi, int size, int ipc) {
    double** dNdX = new double* [size];
    for (size_t i = 0; i < size; i++)
    {
        dNdX[i] = new double[ipc * ipc];
        for (size_t j = 0; j < ipc * ipc; j++)
        {
            dNdX[i][j] = 0;
        }
    }
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < ipc * ipc; j++)
        {
            dNdX[i][j] = reverseJacobian[0][j] * dNdKsi[i][j] + reverseJacobian[1][j] * dNdEta[i][j];
        }
    }
    return dNdX;
}

double** fill_dNdY(double** reverseJacobian, double** dNdEta, double** dNdKsi, int size, int ipc) {
    double** dNdY = new double* [size];
    for (size_t i = 0; i < size; i++)
    {
        dNdY[i] = new double[ipc * ipc];
        for (size_t j = 0; j < ipc * ipc; j++)
        {
            dNdY[i][j] = 0;
        }
    }
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < ipc * ipc; j++)
        {
            dNdY[i][j] = reverseJacobian[2][j] * dNdKsi[i][j] + reverseJacobian[3][j] * dNdEta[i][j];
        }
    }
    return dNdY;
}

void show_dNdX(double** dNdX, int size, int ipc) {
    cout << endl << "dNdX" << endl;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < ipc * ipc; j++)
        {
            cout << "[" << i << "][" << j << "] = " << dNdX[i][j] << " ";
        }
        cout << endl;
    }
}
void show_dNdY(double** dNdY, int size, int ipc) {
    cout << endl << "dNdY" << endl;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < ipc * ipc; j++)
        {
            cout << "[" << i << "][" << j << "] = " << dNdY[i][j] << " ";
        }
        cout << endl;
    }
}

double** transpose_matrix(int size, int ipc, double** Table) {
    double** Transposed_Table = new double* [ipc * ipc];
    for (size_t i = 0; i < ipc * ipc; i++)
    {
        Transposed_Table[i] = new double[size];
        for (size_t j = 0; j < size; j++)
        {
            Transposed_Table[i][j] = 0;
        }
    }
    for (size_t i = 0; i < ipc * ipc; i++)
    {
        for (size_t j = 0; j < size; j++)
        {
            Transposed_Table[i][j] = Table[j][i];
        }
    }
    return Transposed_Table;

}
void show_transposed_dNdX(double** transposed_dNdX, int size, int ipc) {
    cout << endl << "Transposed dNdX" << endl;
    for (int i = 0; i < ipc * ipc; i++)
    {
        for (int j = 0; j < size; j++)
        {
            cout << "[" << i << "][" << j << "] = " << transposed_dNdX[i][j] << " ";
        }
        cout << endl;
    }
}
void show_transposed_dNdY(double** transposed_dNdY, int size, int ipc) {
    cout << endl << "Transposed dNdY" << endl;
    for (int i = 0; i < ipc * ipc; i++)
    {
        for (int j = 0; j < size; j++)
        {
            cout << "[" << i << "][" << j << "] = " << transposed_dNdY[i][j] << " ";
        }
        cout << endl;
    }
}

double*** iloczyn(double** transposed_dNdI, double** dNdI, int size, int ipc) {
    double*** dNdI_il = new double** [ipc* ipc];
    for (size_t i = 0; i < ipc*ipc; i++)
    {
        dNdI_il[i] = new double* [size];
        for (size_t j = 0; j < size; j++)
        {
            dNdI_il[i][j] = new double[size];
            for (size_t k = 0; k < size; k++)
            {
                dNdI_il[i][j][k] = 0;
            }
        }
    }
    for (int i = 0; i < ipc * ipc; i++)
    {
        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                dNdI_il[i][j][k] = dNdI[j][i] * transposed_dNdI[i][k];
            }
        }
    }
    return dNdI_il;
}
void show_iloczynX(double*** dNdX_il, int size, int ipc) {
    cout << endl << "//////Iloczyny X//////" << endl;
    for (size_t i = 0; i < ipc * ipc; i++)
    {
        cout << endl << "Iloczyn X" << i << ": " << endl;
        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                cout << "[" << j << "]" << "[" << k << "] = " << dNdX_il[i][j][k] << " ";
            }
            cout << endl;
        }
    }
}
void show_iloczynY(double*** dNdY_il, int size, int ipc) {
    cout << endl << "//////Iloczyny Y//////" << endl;
    for (size_t i = 0; i < ipc * ipc; i++)
    {
        cout << endl << "Iloczyn Y" << i << ": " << endl;
        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                cout << "[" << j << "]" << "[" << k << "] = " << dNdY_il[i][j][k] << " ";
            }
            cout << endl;
        }
    }

}
double*** count_multiplications(double Factor, double*** dNdX_il, double*** dNdY_il, double* DetJ, int size, int ipc) {
    double*** Multiplications = new double** [ipc * ipc];
    for (size_t i = 0; i < ipc * ipc; i++)
    {
        Multiplications[i] = new double* [size];
        for (size_t j = 0; j < size; j++)
        {
            Multiplications[i][j] = new double[size];
            for (size_t k = 0; k < size; k++)
            {
                Multiplications[i][j][k] = 0;
            }
        }
    }
    for (int i = 0; i < ipc * ipc; i++)
    {
        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                Multiplications[i][j][k] = (dNdX_il[i][j][k] + dNdY_il[i][j][k]) * Factor * DetJ[i];
            }
        }
    }
    return Multiplications;
}
void show_multiplications(double*** multiplications, int size, int ipc) {
    for (size_t i = 0; i < ipc * ipc; i++)
    {
        cout << endl << "Macierz H" << i << ": " << endl;
        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                cout << "[" << j << "]" << "[" << k << "] = " << multiplications[i][j][k] << " ";
            }
            cout << endl;
        }
    }
}
double** count_H(double*** Multiplications, int size, int ipc, double* weights) {
    double** H = new double*[size];
    for (size_t i = 0; i < size; i++)
    {
        H[i] = new double[size];
        for (size_t j = 0; j < size; j++)
        {
            H[i][j] = 0;
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
            for (size_t k = 0; k < ipc * ipc; k++)
            {
                H[i][j] += Multiplications[k][i][j] * Weight[k];
            }
        }
    }
    return H;
}

double distance(Node p1, Node p2) {

    return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

int Gauss(double** uklad_rownan, double* Wyniki, int liczba_wierszy, int dlugosc_wiersza) {
    for (size_t i = 0; i < liczba_wierszy; i++) //sprowadzenie ukladu to postaci trojkatnej 
    {
        for (size_t j = i + 1; j < liczba_wierszy; j++)
        {
            double iloraz = uklad_rownan[j][i] / uklad_rownan[i][i];	//Obliczamy iloraz wartosci jednego wiersza przez drugi 
            uklad_rownan[j][i] = 0;		//zerujemy tworzac "schodki"
            for (size_t k = i + 1; k < dlugosc_wiersza; k++)
            {
                uklad_rownan[j][k] = uklad_rownan[j][k] - uklad_rownan[i][k] * iloraz; //Obliczamy wartosci kolejnych wspolczynnikow w wierszu
            }
        }
    }
    if (uklad_rownan[liczba_wierszy - 1][liczba_wierszy - 1] == 0 && uklad_rownan[liczba_wierszy - 1][liczba_wierszy] != 0) {
        return -1;	//Funkcja zakonczy sie z wynikiem -1 jesli uklad jest sprzeczny
    }
    if (uklad_rownan[liczba_wierszy - 1][liczba_wierszy - 1] == 0 && uklad_rownan[liczba_wierszy - 1][liczba_wierszy] == 0) {
        return 0;	//Funkcja zakonczy sie z wynikiem 0 jesli uklad jest nieoznaczony 
    }
    for (int w = liczba_wierszy - 1; w > -1; w--) {		//Obliczenie wartosci kolejnych X
        Wyniki[w] = uklad_rownan[w][liczba_wierszy];
        for (size_t i = w + 1; i < liczba_wierszy; i++)
        {
            double iloczyn = (uklad_rownan[w][i] * Wyniki[i]);
            double wynik = Wyniki[w];
            int roznica = wynik - iloczyn;
            /*if (w == 0) {
                Wyniki[w] = roznica;
            }
            else {*/
            Wyniki[w] = wynik - iloczyn;
            //}
        }
        Wyniki[w] = Wyniki[w] / uklad_rownan[w][w];
    }
    return 1;
}
double** combineTables(double** matrixH, double* matrixP, int size) {
    double** combinedTable = new double* [size];
    for (size_t i = 0; i < size; i++)
    {
        combinedTable[i] = new double[size + 1];
        for (size_t j = 0; j < size + 1; j++)
        {
            combinedTable[i][j] = 0;
        }
    }
    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = 0; j < size; j++)
        {
            combinedTable[i][j] = matrixH[i][j];
        }
    }
    for (size_t i = 0; i < size; i++)
    {
        combinedTable[i][size] = matrixP[i];
    }
    return combinedTable;
}
void save_global_Combined_Table(string fileName, double** combinedTable, int size) {
    ofstream outputFile(fileName);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size + 1; j++)
        {
            outputFile << combinedTable[i][j] << ";";
        }
        outputFile << endl;
    }
    outputFile.close();
}
void save_globalH(string fileName, double** H, int size) {
    ofstream outputFile(fileName);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            outputFile << H[i][j] << ";";
        }
        outputFile << endl;
    }
    outputFile.close();
}
void show_H(double** H, int size) {
    cout << endl;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            cout << "H[" << i << "][" << j << "] = " << H[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
void show_P(double* P, int size) {
    cout << endl << "P lokalne" << endl;
    for (int i = 0; i < size; i++)
    {
        cout << "P[" << i << "] = " << P[i] << " ";
    }
    cout << endl << endl;
}
void show_Global_Matrix(double** matrix, int size) {
    cout << endl;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
void show_P_Global_Matrix(double* matrix, int size) {
    for (int j = 0; j < size; j++)
    {
        cout << matrix[j] << " ";
    }
    cout << endl;
} 
double findMax(double *Table, int size) {
    double max = Table[0];
    for (size_t i = 0; i < size; i++)
    {
        if (max < Table[i]) {
            max = Table[i];
        }
    }
    return max;
}
double findMin(double* Table, int size) {
    double min = Table[0];
    for (size_t i = 0; i < size; i++)
    {
        if (min > Table[i]) {
            min = Table[i];
        }
    }
    return min;
}