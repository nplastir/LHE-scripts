#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TLine.h"
#include "TMath.h"
#include "TLorentzVector.h"

// Function to calculate the transverse momentum of a particle
double calculate_pt(double px, double py) {
    return sqrt(px * px + py * py);
}

// Function to calculate the invariant mass of a 4-particle system
double calculate_invariant_mass(double px1, double py1, double pz1, double energy1,
                                double px2, double py2, double pz2, double energy2,
                                double px3, double py3, double pz3, double energy3,
                                double px4, double py4, double pz4, double energy4) {
    double total_energy = energy1 + energy2 + energy3 + energy4;
    double total_px = px1 + px2 + px3 + px4;
    double total_py = py1 + py2 + py3 + py4; 
    double total_pz = pz1 + pz2 + pz3 + pz4; 
    double invariant_mass_sq = (total_energy) * (total_energy) 
                             - (total_px) * (total_px) 
                             - (total_py) * (total_py)
                             - (total_pz) * (total_pz);
    return sqrt(std::max(0.0, invariant_mass_sq)); // Ensure non-negative square root
}

// Function to calculate the HT of a 4-particle system
double calculate_HT(double px1, double py1, double px2, double py2, double px3, double py3, double px4, double py4){
    double pt1 = calculate_pt(px1, py1);
    double pt2 = calculate_pt(px2, py2);
    double pt3 = calculate_pt(px3, py3);
    double pt4 = calculate_pt(px4, py4);
    return pt1+pt2+pt3+pt4;
}

// Function that adds the overflow in the last bin 
void add_overflow(TH1F* hist) {
    int nbins = hist->GetNbinsX();
    double e1 = hist->GetBinError(nbins);
    double e2 = hist->GetBinError(nbins + 1);
    hist->AddBinContent(nbins, hist->GetBinContent(nbins + 1));
    hist->SetBinError(nbins, sqrt(e1 * e1 + e2 * e2));
    hist->SetBinContent(nbins + 1, 0);
    hist->SetBinError(nbins + 1, 0);
}

// Function that reads the ROOT file, creates the histograms and saves them as .png
int ttbb(const std::string& inputFileName, const std::string& outputFileName) {
    // Open the ROOT file
    TFile *file = TFile::Open(inputFileName.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // Get the TTree from the file
    TTree *tree;
    file->GetObject("events", tree); // Replace "events" with the actual name of your TTree

    // Set up the branches
    std::vector<float> *px = 0, *py = 0, *pz = 0, *energy = 0, *mass = 0;
    std::vector<int> *pid = 0, *mother1 = 0, *mother2 = 0;
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("energy", &energy);
    tree->SetBranchAddress("pid", &pid);
    tree->SetBranchAddress("mother1", &mother1);
    tree->SetBranchAddress("mother2", &mother2);

    // Create the histograms 
    // TH1F *name = new TH1F(const_char* name, const_char* title, Int_t nbinsx, Double_t xlow, Double_t xup)
    
    TH1F *hist_invariant_mass = new TH1F("invariant_mass_hist", "Invariant Mass; Mass [GeV]; Events", 50, 0, 2500);
    TH1F *hist_Ht = new TH1F("Ht_hist", "HT; pt [GeV]; Events", 50, 0, 1500);

    // Loop over the events in the TTree
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        for (size_t j = 0; j < px->size(); j++) {
            for (size_t k = 0; k < px->size(); k++) {
                for (size_t l = 0; l < px->size(); l++) {
                    for (size_t m = 0; m < px->size(); m++) {
                        if ((*pid)[j] == 6 && (*pid)[k] == -6 && (*pid)[l] == 5 && (*pid)[m] == -5) {
                            if ((*mother1)[l] < 2 && (*mother1)[m] < 2) {
                                // Calculate the invariant mass of the tt+bb system
                                double invariant_mass = calculate_invariant_mass((*px)[j], (*py)[j], (*pz)[j], (*energy)[j],
                                                                                 (*px)[k], (*py)[k], (*pz)[k], (*energy)[k],
                                                                                 (*px)[l], (*py)[l], (*pz)[l], (*energy)[l],
                                                                                 (*px)[m], (*py)[m], (*pz)[m], (*energy)[m]);
                                // Calculate the HT of the tt+bb system
                                double HT = calculate_HT((*px)[j], (*py)[j],
                                                        (*px)[k], (*py)[k], 
                                                        (*px)[l], (*py)[l], 
                                                        (*px)[m], (*py)[m]);
                                // Fill the histograms
                                hist_invariant_mass->Fill(invariant_mass);
                                hist_Ht->Fill(HT);
                            }
                        }
                    }
                }   
            }
        }
    }

    // Handle overflow bins
    add_overflow(hist_invariant_mass);
    add_overflow(hist_Ht);

    // Save histograms to a new ROOT file
    TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    hist_invariant_mass->Write();
    hist_Ht->Write();
    outputFile->Close();

    TCanvas *canvas = new TCanvas("canvas", "Invariant Mass Distribution", 800, 600);
    hist_invariant_mass->Draw("hist");
    canvas->SetGrid();
    canvas->SaveAs("invariant_mass_4_distribution.png");

    TCanvas *canvas2 = new TCanvas("canvas2", "HT Distribution", 800, 600);
    hist_Ht->Draw("hist");
    canvas2->SetGrid();
    canvas2->SaveAs("HT.png");

    return 0;
}

void run_ttbb(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: root -l -b -q 'ttbb.C(\"input_file.root\", \"output_file.root\")'" << std::endl;
        return;
    }

    std::string inputFileName = argv[1];
    std::string outputFileName = argv[2];

    ttbb(inputFileName, outputFileName);
}

void ttbb(const char* inputFileName, const char* outputFileName) {
    ttbb(std::string(inputFileName), std::string(outputFileName));
}

int main(int argc, char** argv) {
    TApplication theApp("App", &argc, argv);
    run_ttbb(argc, argv);
    theApp.Run();
    return 0;
}
