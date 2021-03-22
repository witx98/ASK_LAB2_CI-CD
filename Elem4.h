#include <iostream>

using namespace std;

class Elem4 {
public:
    double* ksi;
    double* eta;
    double* weights;
    int size;
    int ipc;
    double** dNdEta;
    double** dNdKsi;
    double** N;

    double* BcKsi;
    double* BcEta;
    double** NBc;
public:
    Elem4(int size, int ipc) {
        this->size = size;
        this->ipc = ipc;
        if (this->ipc == 2) {
            ksi = new double[ipc * ipc]{ -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0), -1.0 / sqrt(3.0) };
            eta = new double[ipc * ipc]{ -1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
            weights = new double[ipc] {1, 1};

            BcKsi = new double[ipc * size]{ -1.0 / sqrt(3), 1.0 / sqrt(3), 1, 1, 1.0 / sqrt(3), -1.0 / sqrt(3), -1, -1 };
            BcEta = new double[ipc * size]{ -1, -1, -1.0 / sqrt(3), 1.0 / sqrt(3), 1, 1, 1.0 / sqrt(3), -1.0 / sqrt(3) };
        }
        else if (this->ipc == 3) {
            ksi = new double[ipc * ipc]{ -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0) };
            eta = new double[ipc * ipc]{ -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), 0, 0, 0, sqrt(3.0 / 5.0), sqrt(3.0 / 5.0), sqrt(3.0 / 5.0) };
            weights = new double[ipc] {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

            BcKsi = new double[ipc * size]{ -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0), 1, 1, 1, sqrt(3.0 / 5.0), 0, -sqrt(3.0 / 5.0), -1, -1, -1 };
            BcEta = new double[ipc * size]{ -1, -1, -1, -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0), 1, 1, 1, sqrt(3.0 / 5.0), 0, -sqrt(3.0 / 5.0) };
        }
        else if (this->ipc == 4) {
            double valF = sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
            double valS = sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
            ksi = new double[ipc * ipc]{ -valF, -valS, valS, valF, -valF, -valS, valS, valF, -valF, -valS, valS, valF, -valF, -valS, valS, valF };
            eta = new double[ipc * ipc]{ -valF, -valF, -valF, -valF, -valS, -valS, -valS, -valS, valS, valS, valS, valS, valF, valF, valF, valF };

            double weight1 = (18.0 - sqrt(30)) / 36.0;
            double weight2 = (18.0 + sqrt(30)) / 36.0;
            weights = new double[ipc] {weight1, weight2, weight2, weight1};

            BcKsi = new double[ipc * 4]{ -valF, -valS, valS, valF, 1, 1, 1, 1, valF, valS, -valS, -valF, -1, -1, -1, -1 };
            BcEta = new double[ipc * 4]{ -1, -1, -1, -1, -valF, -valS, valS, valF, 1, 1, 1, 1, valF, valS, -valS, -valF };
        }
        this->dNdEta = new double* [this->size];
        this->dNdKsi = new double* [this->size];
        this->N = new double* [this->size];

        for (size_t i = 0; i < this->size; i++)
        {
            this->dNdEta[i] = new double[this->ipc * this->ipc];
            this->dNdKsi[i] = new double[this->ipc * this->ipc];
            this->N[i] = new double[this->ipc * this->ipc];

            for (size_t j = 0; j < this->ipc * this->ipc; j++)
            {
                this->dNdEta[i][j] = 0;
                this->dNdKsi[i][j] = 0;
                this->N[i][j] = 0;


            }
        }
        this->NBc = new double* [this->size];
        for (size_t i = 0; i < this->size; i++)
        {
            this->NBc[i] = new double[this->ipc * this->size];

            for (size_t j = 0; j < this->ipc * this->size; j++)
            {
                this->NBc[i][j] = 0;
            }
        }
        fill_Ksi_Eta_Tables();
        fill_N_Table();
        fill_NBc_Table();
    }
    void fill_Ksi_Eta_Tables() {
        for (int i = 0; i < this->ipc * this->ipc; i++)
        {
            for (int j = 0; j < size; j++)
            {
                if (j == 0)
                {
                    this->dNdKsi[j][i] = (-0.25) * (1.0 - this->eta[i]);
                    this->dNdEta[j][i] = (-0.25) * (1.0 - this->ksi[i]);
                }
                else if (j == 1)
                {
                    this->dNdKsi[j][i] = (-0.25) * (1.0 + this->eta[i]);
                    this->dNdEta[j][i] = (0.25) * (1.0 - this->ksi[i]);
                }
                else if (j == 2)
                {
                    this->dNdKsi[j][i] = (0.25) * (1.0 + this->eta[i]);
                    this->dNdEta[j][i] = (0.25) * (1.0 + this->ksi[i]);
                }
                else if (j == 3)
                {
                    this->dNdKsi[j][i] = (0.25) * (1.0 - this->eta[i]);
                    this->dNdEta[j][i] = (-0.25) * (1.0 + this->ksi[i]);
                }
            }
        }
    }
    void show_Ksi_Eta_Tables() {
        cout << "\n" << "Eta Table" << endl;
        for (int i = 0; i < this->size; i++)
        {
            for (int j = 0; j < this->ipc * this->ipc; j++)
            {
                cout << "[" << i << "][" << j << "] = " << this->dNdEta[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "\n";
        cout << "\n" << "Ksi Table" << endl;
        for (int i = 0; i < this->size; i++)
        {
            for (int j = 0; j < this->ipc * this->ipc; j++)
            {
                cout << "[" << i << "][" << j << "] = " << this->dNdKsi[i][j] << " ";
            }
            cout << "\n";
        }
    }
    void fill_N_Table() {
        for (int i = 0; i < this->ipc * this->ipc; i++)
        {
            for (int j = 0; j < this->size; j++)
            {
                if (j == 0)
                {
                    this->N[j][i] = (0.25) * ((1.0 - this->ksi[i]) * (1.0 - this->eta[i]));
                }
                else if (j == 1)
                {
                    this->N[j][i] = (0.25) * ((1.0 + this->ksi[i]) * (1.0 - this->eta[i]));
                }
                else if (j == 2)
                {
                    this->N[j][i] = (0.25) * ((1.0 + this->ksi[i]) * (1.0 + this->eta[i]));
                }
                else if (j == 3)
                {
                    this->N[j][i] = (0.25) * ((1.0 - this->ksi[i]) * (1.0 + this->eta[i]));
                }
            }
        }
    }
    void show_N_Table() {
        cout << "\n" << "Table N" << endl;
        for (int i = 0; i < this->size; i++)
        {
            for (int j = 0; j < this->ipc * this->ipc; j++)
            {
                cout << "[" << i << "][" << j << "] = " << this->N[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "\n";
    }
    void fill_NBc_Table() {
        for (int i = 0; i < this->size * this->ipc; i++)
        {
            for (int j = 0; j < this->size; j++)
            {
                if (j == 0)
                {
                    this->NBc[j][i] = (0.25) * ((1.0 - this->BcKsi[i]) * (1.0 - this->BcEta[i]));
                }
                else if (j == 1)
                {
                    this->NBc[j][i] = (0.25) * ((1.0 + this->BcKsi[i]) * (1.0 - this->BcEta[i]));
                }
                else if (j == 2)
                {
                    this->NBc[j][i] = (0.25) * ((1.0 + this->BcKsi[i]) * (1.0 + this->BcEta[i]));
                }
                else if (j == 3)
                {
                    this->NBc[j][i] = (0.25) * ((1.0 - this->BcKsi[i]) * (1.0 + this->BcEta[i]));
                }
            }
        }
    }
};