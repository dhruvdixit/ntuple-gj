#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <TH1.h>

using namespace std;

#define NTRACK_MAX (1U << 15)
#define NCLUSTER_MAX (1U << 15)
#define NMC_TRUTH_MAX (1U << 15)
#define CLUSTER_NMC_TRUTH_MAX 32
#define NCELL 5
namespace {

  void to_sm_nphi(unsigned int &sm, unsigned int &nphi,
      unsigned int n)
  {
    sm = n < 11520 ? n / 1152 :
      n < 12288 ? 10 + (n - 11520) / 384 :
      n < 16896 ? 12 + (n - 12288) / 768 :
      18 + (n - 16896) / 384;
    nphi = sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;
  }

  void to_sm_ieta_iphi(unsigned int &sm, unsigned int &ieta,
      unsigned int &iphi, unsigned int n)
  {
    unsigned int nphi;

    to_sm_nphi(sm, nphi, n);

    const unsigned int n0 =
      sm < 10 ? sm * 1152 :
      sm < 12 ? 11520 + (sm - 10) * 384 :
      sm < 18 ? 12288 + (sm - 12) * 768 :
      16896 + (sm - 18) * 384;
    const unsigned int n1 = n - n0;

    ieta = 2 * (n1 / (2 * nphi)) + 1 - (n1 % 2);
    iphi = (n1 / 2) % nphi;
  }

  void neta_nphi(unsigned int &neta, unsigned int &nphi,
      const unsigned int sm)
  {
    neta = sm < 12 ? 48 : sm < 18 ? 32 : 48;
    nphi = sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;
  }

  bool inside_edge(unsigned int n, unsigned int d)
  {
    unsigned int sm;
    unsigned int ieta;
    unsigned int iphi;

    to_sm_ieta_iphi(sm, ieta, iphi, n);

    unsigned int neta;
    unsigned int nphi;

    neta_nphi(neta, nphi, sm);

    return (ieta >= d && iphi >= d &&
        ieta < neta - d && iphi < nphi - d);
  }

}

int main(int argc, char *argv[]){

  //Load the ROOT File and TTree
  TFile *root_file = TFile::Open(argv[1]);

  if (root_file == NULL) {
    fprintf(stderr,"%s: %d: Failed to load ROOT file %s \n",__FILE__,__LINE__,argv[1]);
    exit(EXIT_FAILURE);
  }   

  TDirectoryFile *df = dynamic_cast<TDirectoryFile *>
    (root_file->Get("AliAnalysisTaskNTGJ"));

  if (df == NULL) {
    fprintf(stderr,"%s: %d: Failed to load TDirectory %s \n",__FILE__,__LINE__,"AliAnalysisTaskNTGJ");
    exit(EXIT_FAILURE);
  }   

  TTree *_tree_event = dynamic_cast<TTree *>
    (df->Get("_tree_event"));

  if (_tree_event == NULL) {
    fprintf(stderr,"%s: %d: Failed to load TTree %s \n",__FILE__,__LINE__,"_tree_event");
    exit(EXIT_FAILURE);
  }   

  Float_t multiplicity_v0[64];

  _tree_event->SetBranchAddress("multiplicity_v0", multiplicity_v0);

  UInt_t ncluster;
  Float_t cluster_e[NCLUSTER_MAX];
  Float_t cluster_pt[NCLUSTER_MAX];
  Float_t cluster_eta[NCLUSTER_MAX];
  Float_t cluster_phi[NCLUSTER_MAX];
  Float_t cluster_lambda_square[NCLUSTER_MAX][2];
  Float_t cluster_tof[NCLUSTER_MAX];
  Int_t cluster_ncell[NCLUSTER_MAX];
  UShort_t cluster_cell_id_max[NCLUSTER_MAX];
  Float_t cluster_e_max[NCLUSTER_MAX];
  Float_t cluster_e_cross[NCLUSTER_MAX];
  Float_t cluster_s_nphoton[NCLUSTER_MAX][4];
  UInt_t cluster_nmc_truth[NCLUSTER_MAX];
  UShort_t cluster_mc_truth_index[NCLUSTER_MAX][CLUSTER_NMC_TRUTH_MAX];

  _tree_event->SetBranchAddress("ncluster", &ncluster);
  _tree_event->SetBranchAddress("cluster_e", cluster_e);
  _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
  _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
  _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
  _tree_event->SetBranchAddress("cluster_lambda_square",
      cluster_lambda_square);
  _tree_event->SetBranchAddress("cluster_tof", cluster_tof);
  _tree_event->SetBranchAddress("cluster_ncell",
      cluster_ncell);
  _tree_event->SetBranchAddress("cluster_cell_id_max",
      cluster_cell_id_max);
  _tree_event->SetBranchAddress("cluster_e_max",
      cluster_e_max);
  _tree_event->SetBranchAddress("cluster_e_cross",
      cluster_e_cross);
  _tree_event->SetBranchAddress("cluster_nmc_truth",
      cluster_nmc_truth);
  _tree_event->SetBranchAddress("cluster_mc_truth_index",
      cluster_mc_truth_index);

  UInt_t ntrack;
  Float_t track_e[NTRACK_MAX];
  Float_t track_pt[NTRACK_MAX];
  Float_t track_eta[NTRACK_MAX];
  Float_t track_phi[NTRACK_MAX];
  UChar_t track_quality[NTRACK_MAX];
  Float_t track_eta_emcal[NTRACK_MAX];
  Float_t track_phi_emcal[NTRACK_MAX];

  _tree_event->SetBranchAddress("ntrack", &ntrack);
  _tree_event->SetBranchAddress("track_e", track_e);
  _tree_event->SetBranchAddress("track_pt", track_pt);
  _tree_event->SetBranchAddress("track_eta", track_eta);
  _tree_event->SetBranchAddress("track_phi", track_phi);
  _tree_event->SetBranchAddress("track_quality",
      track_quality);
  _tree_event->SetBranchAddress("track_eta_emcal",
      track_eta_emcal);
  _tree_event->SetBranchAddress("track_phi_emcal",
      track_phi_emcal);

  Float_t cell_e[17664];
  Float_t cell_tof[17664];

  _tree_event->SetBranchAddress("cell_e", cell_e);
  _tree_event->SetBranchAddress("cell_tof", cell_tof);

  UInt_t nmc_truth;
  Float_t mc_truth_e[NMC_TRUTH_MAX];
  Short_t mc_truth_first_parent_pdg_code[NMC_TRUTH_MAX];

  _tree_event->SetBranchAddress("nmc_truth", &nmc_truth);
  _tree_event->SetBranchAddress("mc_truth_e", mc_truth_e);
  _tree_event->SetBranchAddress("mc_truth_first_parent_pdg_code",
      mc_truth_first_parent_pdg_code);


  TF1 *prompt_spline = new TF1("prompt_spline","[0]+[1]*TMath::Log(x) + [2]*pow(TMath::Log(x),2) + [3]*pow(TMath::Log(x),3) + [4]*pow(TMath::Log(x),4)",0,50);
  TF1 *non_prompt_spline = new TF1("non_prompt_spline","[0]+[1]*TMath::Log(x) + [2]*pow(TMath::Log(x),2) + [3]*pow(TMath::Log(x),3) + [4]*pow(TMath::Log(x),4)",0,50);

  size_t count_prompt = 0;
  size_t count_nonprompt = 0;

  TH1F *prompt_pT = new TH1F("prompt_pT","prompt pT distribution",400,0,100);
  TH1F *non_prompt_pT = new TH1F("non_prompt_pT","non_prompt pT distribution",400,0,100);
  /* for (Long64_t i = 0; i < 5000; i++) { */
  for (Long64_t i = 0; i < _tree_event->GetEntries(); i++) {
    _tree_event->GetEntry(i);

    if (i % 1000 == 0) {
      fprintf(stderr, "%s:%d: %lld / %lld "
          "prompt = %lu, nonprompt = %lu)\n",
          __FILE__, __LINE__, i,
          _tree_event->GetEntries(), 
          count_prompt, count_nonprompt);
    }

    for (UInt_t j = 0; j < ncluster; j++) {
      if (cluster_e[j] >= 8.0F &&
          inside_edge(cluster_cell_id_max[j],
            (NCELL - 1) / 2) &&
          // EMCAL noise suppression cut (with no particular
          // name by the ALICE experiment)
          cluster_ncell[j] >= 2 &&
          // EMCAL noise suppression cut, which the ALICE
          // collaboration calls the "exoticity cut"
          cluster_e_cross[j] > 0.05 * cluster_e_max[j]) {

        // Loop over all MC (ground) truth particles and
        // decide if the particle is prompt by the highest
        // energy one hitting the cluster
        float mc_truth_e_max = -INFINITY;
        bool prompt = false;

        for (size_t k = 0; k < cluster_nmc_truth[j]; k++) {
          if (mc_truth_e[cluster_mc_truth_index[j][k]] >
              mc_truth_e_max) {
#if 0
            fprintf(stderr, "%s:%d: %d\n",
                __FILE__, __LINE__,
                mc_truth_first_parent_pdg_code
                [cluster_mc_truth_index[j][k]]);
#endif
            // Prompt clusters are from partons, gauge
            // bosons and leptons (electrons). The PDG
            // numbering scheme of Monte Carlo
            // generator particles
            // (http://pdg.lbl.gov/mc_particle_id_contents.html)
            // designates partons, gauge bosons, and
            // leptons to have absolute values < 100
            // and mesons to have absolute values >=
            // 100
            prompt = std::abs(
                mc_truth_first_parent_pdg_code
                [cluster_mc_truth_index[j][k]])
              < 100;
            mc_truth_e_max = mc_truth_e
              [cluster_mc_truth_index[j][k]];

          }//Truth E Check
        }//Cluster MC Truth Index

        // Guard against clusters with empty MC (ground)
        // truth particle association
        if (mc_truth_e_max >= 0) {
          if (prompt) {
            count_prompt++;
            prompt_pT->Fill(cluster_pt[j]);
          }
          else {
            count_nonprompt++;
            non_prompt_pT->Fill(cluster_pt[j]);
            /* std::cout<<"cool"<<std::endl; */
          }
        } 

      }//Cluster Cuts
    }//Clusters
  }//Events
  prompt_pT->Fit(prompt_spline,"S");
  non_prompt_pT->Fit(non_prompt_spline,"S");

  TFile *out_root = new TFile("spline_fits.root","RECREATE");
  prompt_pT->Write(); 
  prompt_spline->Write();
  non_prompt_pT->Write();
  non_prompt_spline->Write();
  out_root->Close();
  return 0;
}
