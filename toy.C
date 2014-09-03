#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include "TRandom3.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"
#include "TF1.h"
#include "TFile.h"
#include "TString.h"

#define PI 3.14159265

using namespace std;

void p4(TLorentzVector v, TString name="") {
    cout << name << endl;
    cout << "(px,py,pz,e): (" << v.Px() << "," << v.Py() << "," << v.Pz() << "," << v.E() << ") " 
        << "(pt,eta,phi,et): " << v.Pt() << "," << v.Eta() << "," << v.Phi() << "," << v.Et() << endl << endl;
}

float deltaR(TLorentzVector v1, TLorentzVector v2) {
    float dPhi = min(fabs(v1.Phi()-v2.Phi()), 2.0*PI-fabs(v1.Phi()-v2.Phi()));
    float dEta = fabs(v1.Eta() - v2.Eta());
    return sqrt(dEta*dEta + dPhi*dPhi);
}

void doStuff() {

    TRandom3 r;
    TF1 f1("f1","sin(x)",0,PI);

    TLorentzVector top, b, W, q1, q2;
    float dRqq, dRq1b, dRq2b;


    TFile f("data.root","recreate");
    // f.SetCompressionLevel(1); 
    TTree tree("Events","Events");
    tree.Branch("Wp4","TLorentzVector",&W);
    tree.Branch("bp4","TLorentzVector",&b);
    tree.Branch("topp4","TLorentzVector",&top);
    tree.Branch("q1p4","TLorentzVector",&q1);
    tree.Branch("q2p4","TLorentzVector",&q2);

    tree.Branch("dRqq", &dRqq, "dRqq/F");
    tree.Branch("dRq1b", &dRq1b, "dRq1b/F");
    tree.Branch("dRq2b", &dRq2b, "dRq2b/F");

    double massW = 80.4;
    double massQuark = 0.01; // 10 MeV let's say XXX (?)
    double massTop = 173; 
    double massB = 0.042; // 4.2 MeV

    for(double energyTop = 400; energyTop < 1200; energyTop += 25) {

        cout << "energyTop: " << energyTop << endl;

        // double energyTop = 700;

        for(int i = 0; i < 20000; i++) {
            // TLorentzVector top;
            double pTop = sqrt(energyTop*energyTop - massTop*massTop);
            double phiTop = r.Uniform(-PI,PI);
            double thetaTop = r.Uniform(0,PI);
            top.SetPxPyPzE( pTop*cos(phiTop)*sin(thetaTop), pTop*sin(phiTop)*sin(thetaTop), pTop*cos(thetaTop), energyTop );
            TVector3 boostFromTop = top.BoostVector();

            // TLorentzVector b;
            double energyB = (massTop*massTop - massW*massW + massB*massB)/(2*massTop);
            double pB = sqrt(energyB*energyB - massB*massB);
            double phiB = r.Uniform(-PI,PI);
            double thetaB = r.Uniform(0,PI);
            b.SetPxPyPzE( pB*cos(phiB)*sin(thetaB), pB*sin(phiB)*sin(thetaB), pB*cos(thetaB), energyB );
            b.Boost(boostFromTop);

            // TLorentzVector W;
            double energyW = (massTop*massTop + massW*massW - massB*massB)/(2*massTop); // W is unboosted right now
            double pW = sqrt(energyW*energyW - massW*massW);
            double phiW = -phiB;
            double thetaW = -thetaB;
            W.SetPxPyPzE( pW*cos(phiW)*sin(thetaW), pW*sin(phiW)*sin(thetaW), pW*cos(thetaW), energyW );

            W.Boost(boostFromTop);

            TVector3 boostFromW = W.BoostVector();

            double energyQuark = W.E() / 2; // let's say, if qq masses same
            double pQuark = sqrt(energyQuark*energyQuark - massQuark*massQuark);
            double phiQuark = r.Uniform(-PI,PI);
            double thetaQuark = f1.GetRandom();

            //back to back qq
            // TLorentzVector q1, q2;
            q1.SetPxPyPzE( pQuark*cos(phiQuark)*sin(thetaQuark), pQuark*sin(phiQuark)*sin(thetaQuark), pQuark*cos(thetaQuark), energyQuark );
            q2.SetPxPyPzE(-pQuark*cos(phiQuark)*sin(thetaQuark),-pQuark*sin(phiQuark)*sin(thetaQuark),-pQuark*cos(thetaQuark), energyQuark );

            q1.Boost(boostFromW);
            q2.Boost(boostFromW);

            dRqq = deltaR(q1,q2);
            dRq1b = deltaR(q1,b);
            dRq2b = deltaR(q2,b);

            // p4(q1, "q1");
            // p4(q2, "q2");
            // p4(W, "W");
            // p4(b, "b");
            // p4(top, "top");


            // std::cout << " deltaR(q1,q2): " << deltaR(q1,q2) << std::endl;
            // std::cout << " deltaR(q1,b): " << deltaR(q1,b)  << std::endl;
            // std::cout << " deltaR(q2,b): " << deltaR(q2,b)  << std::endl;

            tree.Fill();


        }
    }
    tree.Print();
    tree.Write();
}


void toy() {

    // TLorentzVector * q1, q2, b, W, top;
    doStuff();
}
