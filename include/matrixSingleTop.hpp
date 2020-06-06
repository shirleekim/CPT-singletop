#ifndef matrixSingleTop_h
#define matrixSingleTop_h
#define matrixSingleTop_cxx

#include "../include/const.hpp"

#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
#include <TString.h>

#ifndef __CINT__
//___________Matrix Element class for SM and SME in single-top production case___________//
class MatrixSingleTop {

     private:

//***********Initiales values*************//

    // P
    TLorentzVector pMother1;
    TLorentzVector b;
    TLorentzVector pt;
    TLorentzVector p3;
    int nature;
    int pMother1_PID;
    int p3_PID;
    int b_PID;
    int pt_PID; // ajoute par sukyung
    double pTelec;
    double etaElec;
    double pTmu;
    double etamu;
    double pTb;
    double eta_b;
    double eta_non_b;
    double pTnon_b;

    public:
//**********Calculated Values************//

    //For SM
    double Mbg;       double calculateMbg();
    double Mbq;       double calculateMbq();
    double Mbqbar;    double calculateMbqbar();

    //For SME
    TLorentzVector Mbg_mu;           TLorentzVector calculateMbg_mu();
    TLorentzVector Mbq_mu;           TLorentzVector calculateMbq_mu();
    TLorentzVector Mbqbar_mu;        TLorentzVector calculateMbqbar_mu();

//*************** Constructor *************//

        // Takes the 4-momenti and creat MatrixSingleTop,

            MatrixSingleTop();
            ~MatrixSingleTop();

            MatrixSingleTop (TLorentzVector pMother1_user, TLorentzVector b_user, TLorentzVector pt_user, TLorentzVector p3_user, int pMother1_PID_user, int p3_PID_user, int b_PID_user, int nature_user,double pTmu_user, double etamu_user, double pTelec_user, double etaElec_user, double pTb_user, double pTnon_b_user, double eta_b_user, double eta_non_b_user);

//****** Recover values previously obtained ******//

    double getMbg();
    double getMbq();
    double getMbqbar();
    TLorentzVector getMbg_mu();
    TLorentzVector getMbq_mu();
    TLorentzVector getMbqbar_mu();
    void checkckm();

//****** Tools to simplify the expression of the matrix elements *******//

    static double normMinkowski2 (TLorentzVector p);
    //static void writeNumMatrix (Tstring s, int n, TLorentzVector m);
};

#endif
#endif
