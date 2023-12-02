#include <TH2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>

void final_plots()
{


TFile *f_dat = TFile::Open("C:/root_v6.24.06/macros/OpenData/Data_Histograms.root","read");
TH1F *h_phi_data = (TH1F*)f_dat->Get("h_phi_mu")->Clone();
TH1F *h_eta_data = (TH1F*)f_dat->Get("h_eta_mu")->Clone();
TH1F *h_pt_data = (TH1F*)f_dat->Get("h_pt_mu")->Clone();
h_phi_data->SetDirectory(0);
h_eta_data->SetDirectory(0);
h_pt_data->SetDirectory(0);
TFile *f_sig = TFile::Open("C:/root_v6.24.06/macros/OpenData/Signal_Histograms.root","read");
TH1F *h_phi_signal = (TH1F*)f_sig->Get("h_phi_mu")->Clone();
TH1F *h_eta_signal = (TH1F*)f_sig->Get("h_eta_mu")->Clone();
TH1F *h_pt_signal = (TH1F*)f_sig->Get("h_pt_mu")->Clone();
h_phi_signal->SetDirectory(0);
h_eta_signal->SetDirectory(0);
h_pt_signal->SetDirectory(0);
TFile *f_bkg = TFile::Open("C:/root_v6.24.06/macros/OpenData/Bkg_Histograms.root","read");
TH1F *h_phi_bkg = (TH1F*)f_bkg->Get("h_phi_mu")->Clone();
TH1F *h_eta_bkg = (TH1F*)f_bkg->Get("h_eta_mu")->Clone();
TH1F *h_pt_bkg = (TH1F*)f_bkg->Get("h_pt_mu")->Clone();
h_phi_bkg->SetDirectory(0);
h_eta_bkg->SetDirectory(0);
h_pt_bkg->SetDirectory(0);

TCanvas *c1 = new TCanvas();
THStack *hs_eta = new THStack("hs_eta","");
h_eta_signal->SetFillColor(kGreen);
h_eta_bkg->SetFillColor(kBlue);
h_eta_data->SetMarkerStyle(20);
hs_eta->Add(h_eta_bkg);
hs_eta->Add(h_eta_signal);
h_eta_data->Draw("e");
hs_eta->Draw("HISTsame");
h_eta_data->Draw("esame");

TCanvas *c2 = new TCanvas();
THStack *hs_phi = new THStack("hs_phi","");
h_phi_signal->SetFillColor(kGreen);
h_phi_bkg->SetFillColor(kBlue);
h_phi_data->SetMarkerStyle(20);
hs_phi->Add(h_phi_bkg);
hs_phi->Add(h_phi_signal);
h_phi_data->Draw("e");
hs_phi->Draw("HISTsame");
h_phi_data->Draw("esame");

TCanvas *c3 = new TCanvas();
THStack *hs_pt = new THStack("hs_pt","");
h_pt_signal->SetFillColor(kGreen);
h_pt_bkg->SetFillColor(kBlue);
h_pt_data->SetMarkerStyle(20);
hs_pt->Add(h_pt_bkg);
hs_pt->Add(h_pt_signal);
h_pt_data->Draw("e");
hs_pt->Draw("HISTsame");
h_pt_data->Draw("esame");
}