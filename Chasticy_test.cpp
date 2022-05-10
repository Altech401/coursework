#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream> 
#define PI 3.14159265

const int N = 1000, R0 = 20, Rc = 30, H = 20, Z = 13;
const double mc2 = 0.51099, Emin = 0.02, Wmin = 0.000000000001, rd1[3] = { 0, 0, 10}, rd2[3] = { 0, 0, 20 }, rd3[3] = { 30, 0, 20 }, rd4[3] = { -30, 0, 20 }, r_e = 2.8 * 0.0000000000001;
using namespace std;
const double sig_comp[19] = { 0.1845, 0.181, 0.176, 0.170,  0.165,  0.160,  0.154,   0.146,  0.135, 0.123, 0.11,   0.097,  0.0905,  0.0835,  0.0735,  0.0645,  0.0555,    0.0466,    0.0374 },
               sig_ph[19] = { 16.61,  5.29,  1.943, 0.575 , 0.2405, 0.1222, 0.06115, 0.0259, 0.01,  0.002, 0.0015, 0.0003, 0.00015, 0.00008, 0.00004, 0.00002, 0.000007,  0.000002,  0.000000005 },
             sig_pair[19] = { 0,      0,     0,     0,      0,      0,      0,       0,      0,     0,     0,      0,      0,       0,       0,       0,       0.0000835, 0.0004185, 0.00129 },
                sigma[19] = { 45.3452,  14.7717,  5.7213,  2.0115,  1.09485 , 0.76194 , 0.580905,  0.46413,  0.3915,  0.3375,  0.30105,  0.26271,  0.244755,  0.225666,  0.198558,  0.174204,  0.150094,  0.126955,  0.104463 };
bool Fl = false;
int N0=N;

      

/*
struct Chst
{
    float N0, cos_theta[N], E[N], x[N], y[N], z[N], Omega[N][3], W[N]; int K[N];
};
*/
//Структуры
/*struct Chst
{ 
    int NN=30;
    double* cos_theta = new double[NN];
    double* E = new double[NN];
    double* x = new double[NN];
    double* y = new double[NN];
    double* z = new double[NN];
    double* W = new double[NN];
    double* Omega1 = new double[NN];
    double* Omega2 = new double[NN];
    double* Omega3 = new double[NN];
    int* K= new int[30];
    bool FL1 = true;
};
*/

struct Chst
{
    double cos_theta[30], E[30], x[30], y[30], z[30], Omega1[30], Omega2[30], Omega3[30],W[30];
    int K[30];
    bool FL1 = true;
};

struct ETA
{
    double Eta1[19];
    double Eta2[19];
    double Eta3[19];
    double Eta4[19];
};
/*
struct FAIL
{
    double* x = new double[N];
    double* y = new double[N];
    double* z = new double[N];
    bool FL1 = true;
};

*/

typedef struct ETA Struct2;
typedef struct Chst Struct1;

//Struct1 A;

//FAIL B[N];


Struct1 energy_group(Struct1 s,int i)
{
    
        if ((0.01 <= s.E[i]) and (s.E[i] <= 0.015))
        {
            s.K[i] = 0;
        }
        if ((0.015 < s.E[i]) and (s.E[i] <= 0.02))
        {
            s.K[i] = 1;
        }
        if ((0.02 < s.E[i]) and (s.E[i] <= 0.03))
        {
            s.K[i] = 2;
        }
        if ((0.03 < s.E[i]) and (s.E[i] <= 0.04))
        {
            s.K[i] = 3;
        }
        if ((0.04 < s.E[i]) and (s.E[i] <= 0.05))
        {
            s.K[i] = 4;
        }
        if ((0.05 < s.E[i]) and (s.E[i] <= 0.06))
        {
            s.K[i] = 5;
        }
        if ((0.06 < s.E[i]) and (s.E[i] <= 0.08))
        {
            s.K[i] = 6;
        }
        if ((0.08 < s.E[i]) and (s.E[i] <= 0.1))
        {
            s.K[i] = 7;
        }
        if ((0.1 < s.E[i]) and (s.E[i] <= 0.15))
        {
            s.K[i] = 8;
        }
        if ((0.15 < s.E[i]) and (s.E[i] <= 0.2))
        {
            s.K[i] = 9;
        }
        if ((0.2 < s.E[i]) and (s.E[i] <= 0.3))
        {
            s.K[i] = 10;
        }
        if ((0.3 < s.E[i]) and (s.E[i] <= 0.4))
        {
            s.K[i] = 11;
        }
        if ((0.4 < s.E[i]) and (s.E[i] <= 0.5))
        {
            s.K[i] = 12;
        }
        if ((0.5 < s.E[i]) and (s.E[i] <= 0.6))
        {
            s.K[i] = 13;
        }
        if ((0.6 < s.E[i]) and (s.E[i] <= 0.8))
        {
            s.K[i] = 14;
        }
        if ((0.8 < s.E[i]) and (s.E[i] <= 1.0))
        {
            s.K[i] = 15;
        }
        if ((1.0 < s.E[i]) and (s.E[i] <= 1.5))
        {
            s.K[i] = 16;
        }
        if ((1.5 < s.E[i]) and (s.E[i] <= 2.0))
        {
            s.K[i] = 17;
        }
        if ((2.0 < s.E[i]) and (s.E[i] <= 3.0))
        {
            s.K[i] = 18;
        }
    
    return s;
}

Struct1 StartParameters(Struct1 s, int i)
{
    s.E[i] = 2.5;
    s.W[i] = 1;
    //Разыгрывание начальных координат

    s.x[i] = ((double)rand() / (double)RAND_MAX) * 20;
    s.y[i] = ((double)rand() / (double)RAND_MAX) * 40 - 20;
    s.z[i] = 0;
    if (sqrt(pow(s.x[i], 2) + pow(s.y[i], 2)) > R0)
    {
        s = StartParameters(s, i);
    }
    
    return s;
}  

Struct1 FirstItCoord(Struct1 s, int i)
{
    double L,psi,theta;
    int K;
    K = s.K[0];
    s.E[i] = 2.5;
    s.W[i] = 1;
    L = (-log((double)rand() / (double)RAND_MAX)) / (sigma[s.K[0]]);
    psi = 2 * PI * ((double)rand() / (double)RAND_MAX);
    theta = PI / 2 * ((double)rand() / (double)RAND_MAX);
    s.x[i] = s.x[i-1] + L * cos(theta) * cos(psi);
    s.y[i] = s.y[i-1] + L * cos(theta) * sin(psi);
    s.z[i] = s.z[i-1] + L * sin(theta);
    s.Omega1[i] = (s.x[i] - s.x[i - 1]) / L;
    s.Omega2[i] = (s.y[i] - s.y[i - 1]) / L;
    s.Omega3[i] = (s.z[i] - s.z[i - 1]) / L;
    
    return s;
}  

/*
void OutInFile0(Struct1 s)
{
    ofstream XYZ;
    XYZ.open("XYZ0.xyz");

    for (int i = 0; i < *s.N0; i++)
    {
        XYZ << s.x[i] << "   " << s.y[i] << "   " << s.z[i] << endl;
        
    }
    XYZ.close();

}
void OutInFile1(Struct1 s)
{
    ofstream XYZ;
    XYZ.open("XYZ1.xyz");

    for (int i = 0; i < *s.N0; i++)
    {
        XYZ << s.x[i] << "   " << s.y[i] << "   " << s.z[i] << endl;

    }
    XYZ.close();

}
*/

Struct1 Proverka(Struct1 s, int i)
{
    /*int j = 0, N1;
    N1 = N0;
    */
    if ((sqrt(pow(s.x[i], 2) + pow(s.y[i], 2)) > Rc) or (s.z[i] > H) or (s.z[i] < 0) or (s.E[i] < Emin) or (s.W[i]<Wmin))
    {
        s.FL1 = false;
        //s.NN = i;
    }

    /*
    *s.N0 = j;
    
    if (Fl == true)
    {
        cout << "Статистический вес или энергия меньше минимального значения у : " << N1 - *s.N0 <<" частиц"<< endl;
        cout << "Осталось: " << *s.N0 << endl;
        cout << " " << endl;
        Fl = false;
    }
    else
    {
        cout << "Ушло за пределы решения задачи: " << N1 - *s.N0 << endl;
        cout << "Осталось: " << *s.N0 << endl;
        cout << " " << endl;
    }
    */
    return s;
    
}

Struct1 TypeOfIt(Struct1 s, int i)
{
    int n_comp = 0, n_ph=0, n_pair=0, N1,K;
    N1 = N0;
    K = s.K[1];
    if (sig_comp[s.K[i-1]] * 2.7 / sigma[s.K[i-1]] >= ((double)rand() / (double)RAND_MAX))
    {
        n_comp++;
    }
    else
    {
        if ((sig_comp[s.K[i-1]] + sig_ph[s.K[i-1]]) * 2.7 / sigma[s.K[i-1]] >= ((double)rand() / (double)RAND_MAX))
        {
            n_ph++;
            s.FL1 = false;
            //s.NN = i;
        }
        else
        {
            if ((sig_comp[s.K[i-1]] + sig_ph[s.K[i-1]] + sig_pair[s.K[i-1]]) * 2.7 / sigma[s.K[i-1]] >= ((double)rand() / (double)RAND_MAX))
            {
                n_pair++;
                s.FL1 = false;
                //s.NN = i;
            }
            else
            {
                n_pair++;
                s.FL1 = false;
                //s.NN = i;
            }
        }
    }
    
    
    return s;
    
}

Struct1 energyAfterIt(Struct1 s, int i)
{
    double alpha1, x, p, G1, G2;
    bool f;
    
    alpha1 = s.E[i-1] / mc2;
    f = true;
    while (f == true)
    {
        G1 = (double)rand() / (double)RAND_MAX;
        G2 = (double)rand() / (double)RAND_MAX;
        x = alpha1 * (1 + 2 * alpha1 * G1) / (1 + 2 * alpha1);
        p = x / alpha1 + alpha1 / x + (1 / alpha1 - 1 / x) * (2 + 1 / alpha1 - 1 / x);
        if (G2 * (1 + 2 * alpha1 + 1 / (1 + 2 * alpha1)) < p)
        {
            s.E[i] = x * mc2;
            s.cos_theta[i] = 1 - 1 / x + 1 / alpha1;
            f = false;
        }
    }

    
return s;

}

Struct1 StaticalWeihgt(Struct1 s, int i)
{
    s.W[i] = s.W[i-1] * sig_comp[s.K[i]]*2.7 / sigma[s.K[i]];
 
        return s;
}

Struct2 Vklad_v_detector(Struct1 s, Struct2 E, int i)
{
    double alpha, U,sig_c_up, sig_c_down, W_x_sigc, Sigma_Rho, delta_r21, delta_r22, delta_r23, delta_r24,Et1,Et2, Et3, Et4;
    
    alpha = s.E[i] / mc2;
    U = s.cos_theta[i];
    sig_c_up = Z * pow(r_e, 2) / 2 / pow(1+alpha*(1-U),2)*(1+pow(U,2)+pow(alpha*(1-U),2)/(1 + alpha*(1-U)));
    sig_c_down = 2 * PI * Z * pow(r_e, 2) * ((1 + alpha) / pow(alpha, 2) * (2 * (1 + alpha) / (1 + 2 * alpha) - log(1 + 2 * alpha) / alpha) + log(1 + 2 * alpha) / (2 * alpha) - (1 + 3 * alpha) / pow(1 + 2 * alpha, 2));
    W_x_sigc = s.W[i] * sig_c_up / sig_c_down;
    Sigma_Rho = sigma[s.K[i]];
    delta_r21 = pow(s.x[i] - rd1[0], 2) + pow(s.y[i] - rd1[1], 2) + pow(s.z[i] - rd1[2], 2);
    delta_r22 = pow(s.x[i] - rd2[0], 2) + pow(s.y[i] - rd2[1], 2) + pow(s.z[i] - rd2[2], 2);
    delta_r23 = pow(s.x[i] - rd3[0], 2) + pow(s.y[i] - rd3[1], 2) + pow(s.z[i] - rd3[2], 2);
    delta_r24 = pow(s.x[i] - rd4[0], 2) + pow(s.y[i] - rd4[1], 2) + pow(s.z[i] - rd4[2], 2);
    Et1 = (W_x_sigc * exp(-Sigma_Rho * sqrt(delta_r21))) / delta_r21;
    Et2 = (W_x_sigc * exp(-Sigma_Rho * sqrt(delta_r22))) / delta_r22;
    Et3 = (W_x_sigc * exp(-Sigma_Rho * sqrt(delta_r23))) / delta_r23;
    Et4 = (W_x_sigc * exp(-Sigma_Rho * sqrt(delta_r24))) / delta_r24;
    E.Eta1[s.K[i]] += Et1;
    E.Eta2[s.K[i]] += Et2;
    E.Eta3[s.K[i]] += Et3;
    E.Eta4[s.K[i]] += Et4;
    
    return E;
}

Struct1 Generate_new_Interaction(Struct1 s, int i)
{
    double omega01, omega02, omega03, U, psi, L;
    
    psi = 2 * PI * (double)rand() / (double)RAND_MAX;
    U = s.cos_theta[i-1];
    omega01 = s.Omega1[i-1];
    omega02 = s.Omega2[i-1];
    omega03 = s.Omega3[i-1];
    s.Omega3[i] = omega03 * U + sqrt((1 - pow(U, 2)) * (1 - pow(omega03, 2))) * cos(psi);
    s.Omega2[i] = (omega02 * (U - omega03 * s.Omega3[i]) + omega01 * sin(psi) * sqrt((1 - pow(U, 2)) * (1 - pow(omega03, 2)))) / (1 - pow(omega03, 2));
    s.Omega1[i] = (omega01 * (U - omega03 * s.Omega3[i]) - omega02 * sin(psi) * sqrt((1 - pow(U, 2)) * (1 - pow(omega03, 2)))) / (1 - pow(omega03, 2));
    L = -log((double)rand() / (double)RAND_MAX) / sigma[s.K[i-1]];
    s.x[i] = s.x[i-1] + L * s.Omega1[i];
    s.y[i] = s.y[i-1] + L * s.Omega2[i];
    s.z[i] = s.z[i-1] + L * s.Omega3[i];

    return s;
}



//////////////////////////////////////////////////////////////////////////////////

int main()
{
    using namespace std;
    int n = 0;
    ofstream XYZ;
    //XYZ.open("XYZ.xyz");
    setlocale(LC_ALL, "Russian");
    Struct1* A = new Struct1[N];
    Struct2* E = new Struct2;




    srand(time(NULL));

    for (int i = 0; i < 19; i++)
    {
        (*E).Eta1[i] = 0;
        (*E).Eta2[i] = 0;
        (*E).Eta3[i] = 0;
        (*E).Eta4[i] = 0;

    }
    for (int j = 0; j < N; j++)
    {
        n = 0;
        A[j] = StartParameters(A[j],n);
        A[j] = energy_group(A[j],n);
        n = 1;
        A[j] = FirstItCoord(A[j],n);
        A[j] = Proverka(A[j],n);

        if (A[j].FL1 == true)
        {
            
            A[j] = TypeOfIt(A[j],n);
            if (A[j].FL1 == true)
            {
                A[j] = energyAfterIt(A[j], n);
                A[j] = energy_group(A[j], n);
                A[j] = StaticalWeihgt(A[j], n);
                A[j] = Proverka(A[j], n);
                
                if (A[j].FL1 == true)
                {
                    *E = Vklad_v_detector(A[j], *E, n);
                    while (N0 > 0)
                    {
                        n++;
                    
                        A[j] = Generate_new_Interaction(A[j], n);
                        A[j] = energyAfterIt(A[j], n);
                        A[j] = energy_group(A[j], n);
                        A[j] = StaticalWeihgt(A[j], n);
                        A[j] = Proverka(A[j], n);
                        
                        if (A[j].FL1 == true)
                        {
                            
                            A[j] = TypeOfIt(A[j], n);

                            if (A[j].FL1 == true)
                            {
                                A[j] = energyAfterIt(A[j], n);
                                A[j] = energy_group(A[j], n);
                                A[j] = StaticalWeihgt(A[j], n);
                                A[j] = Proverka(A[j], n);
                                if (A[j].FL1 == true)
                                {
                                    *E = Vklad_v_detector(A[j], *E, n);
                                }
                            }
                            else break;
                        }
                        else break;
                    }
                }
            }
        }
    }


    char q[255] = "test";
    char buffer[33];
    for (int i = 0; i < N; i++)
    {
        q[4] = 0;
        _itoa_s(i, buffer, 10);
        strcat_s(q, buffer);
        strcat_s(q, ".txt");
        strcat_s(buffer, ".txt");
        ofstream out(buffer);
        for (int j = 0; j < 30; j++)
        {
            if (A[i].x[j] <= -100 ) break;
            out << A[i].x[j] << "   " << A[i].y[j] << "   " << A[i].z[j] << "   " << endl;
        }
    }


    for (int i = 0; i < 19; i++)
    {
        cout << E->Eta1[i]<< "  ";

    }
    cout << endl;
    for (int i = 0; i < 19; i++)
    {
        cout << E->Eta2[i] << "  ";

    }
    cout << endl;
    for (int i = 0; i < 19; i++)
    {
        cout << E->Eta3[i] << "  ";

    }
    cout << endl;
    for (int i = 0; i < 19; i++)
    {
        cout << E->Eta4[i] << "  ";

    }
    //XYZ.close();
}
