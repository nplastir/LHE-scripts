import ROOT
import os

def overlay_and_ratio_plot(root_files, output_directory):

    ROOT.gStyle.SetOptTitle(0) #removes the hist title
    # Include axes on all sides 
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    # Number of quantities
    quan_name=[
        "all_px", "all_py", "all_pz", "all_energy", "all_mass", "all_pt", "all_rapidity", "all_pseudorapidity", "all_phi", "all_pdgid", "initial_particles_px", "initial_particles_py", "initial_particles_pz", "initial_particles_energy", "initial_particles_mass", "initial_particles_pt", "initial_particles_rapidity", "initial_particles_pseudorapidity", "initial_particles_phi", "initial_particles_pdgid", "radiation_px", "radiation_py", "radiation_pz", "radiation_energy", "radiation_mass", "radiation_pt", "radiation_rapidity", "radiation_pseudorapidity", "radiation_phi", "radiation_pdgid", "top_px", "top_py", "top_pz", "top_energy", "top_mass", "top_pt", "top_rapidity", "top_phi", "antitop_px", "antitop_py", "antitop_pz", "antitop_energy", "antitop_mass", "antitop_pt", "antitop_rapidity", "antitop_phi", "ttbar_px", "ttbar_py", "ttbar_pz", "ttbar_energy", "ttbar_mass", "ttbar_pt", "ttbar_rapidity", "ttbar_phi", "Wp_px", "Wp_py", "Wp_pz", "Wp_energy", "Wp_mass", "Wp_pt", "Wp_rapidity", "Wp_phi", "Wm_px", "Wm_py", "Wm_pz", "Wm_energy", "Wm_mass", "Wm_pt", "Wm_rapidity", "Wm_phi", "W_px", "W_py", "W_pz", "W_energy", "W_mass", "W_pt", "W_rapidity", "W_phi", "b_all_px", "b_all_py", "b_all_pz", "b_all_energy", "b_all_mass", "b_all_pt", "b_all_rapidity", "b_all_pseudorapidity", "b_all_phi", "bottom_px", "bottom_py", "bottom_pz", "bottom_energy", "bottom_mass", "bottom_pt", "bottom_rapidity", "bottom_pseudorapidity", "bottom_phi", "antibottom_px", "antibottom_py", "antibottom_pz", "antibottom_energy", "antibottom_mass", "antibottom_pt", "antibottom_rapidity", "antibottom_pseudorapidity", "antibottom_phi", "b_from_top_px", "b_from_top_py", "b_from_top_pz", "b_from_top_energy", "b_from_top_mass", "b_from_top_pt", "b_from_top_rapidity", "b_from_top_pseudorapidity", "b_from_top_phi", "b_from_initial_px", "b_from_initial_py", "b_from_initial_pz", "b_from_initial_energy", "b_from_initial_mass", "b_from_initial_pt", "b_from_initial_rapidity", "b_from_initial_pseudorapidity", "b_from_initial_phi", "tt_system_pt", "tt_system_rapidity", "tt_system_HT", "tt_system_invariant_mass", "tt_system_dR", "tt_system_dR_2", "bb_system_pt", "bb_system_rapidity", "bb_system_pseudorapidity", "bb_system_HT", "bb_system_invariant_mass", "bb_system_dR", "prompt_bs_pt", "prompt_bs_rapidity","prompt_bs_pseudorapidity", "prompt_bs_HT", "prompt_bs_invariant_mass", "prompt_bs_dR", "bs_from_top_pt", "bs_from_top_rapidity", "bs_from_top_pseudorapidity", "bs_from_top_HT", "bs_from_top_invariant_mass", "bs_from_top_dR", "top_reco_mass", "antitop_reco_mass", "ttbar_reco_mass", "Wp_reco_mass", "Wm_reco_mass", "W_reco_mass", "quarks_system_pt", "quarks_system_rapidity", "quarks_system_pseudorapidity", "quarks_system_HT", "quarks_system_invariant_mass", "quarks_system_dR", "ud_pairs_system_pt", "ud_pairs_system_rapidity", "ud_pairs_system_pseudorapidity", "ud_pairs_system_HT", "ud_pairs_system_invariant_mass", "ud_pairs_system_dR", "sc_pairs_system_pt", "sc_pairs_system_rapidity", "sc_pairs_system_pseudorapidity", "sc_pairs_system_HT", "sc_pairs_system_invariant_mass", "sc_pairs_system_dR", "cd_pairs_system_pt", "cd_pairs_system_rapidity", "cd_pairs_system_pseudorapidity", "cd_pairs_system_HT", "cd_pairs_system_invariant_mass", "cd_pairs_system_dR", "su_pairs_system_pt", "su_pairs_system_rapidity", "su_pairs_system_pseudorapidity", "su_pairs_system_HT", "su_pairs_system_invariant_mass", "su_pairs_system_dR", "quarks_px", "quarks_py", "quarks_pz", "quarks_energy", "quarks_mass", "quarks_pt", "quarks_rapidity", "quarks_pseudorapidity", "quarks_phi", "quarks_reco_mass", "d_quarks_px", "d_quarks_py", "d_quarks_pz", "d_quarks_energy", "d_quarks_mass", "d_quarks_pt", "d_quarks_rapidity", "d_quarks_pseudorapidity", "d_quarks_phi", "d_quarks_reco_mass", "u_quarks_px", "u_quarks_py", "u_quarks_pz", "u_quarks_energy", "u_quarks_mass", "u_quarks_pt", "u_quarks_rapidity", "u_quarks_pseudorapidity", "u_quarks_phi", "u_quarks_reco_mass", "c_quarks_px", "c_quarks_py", "c_quarks_pz", "c_quarks_energy", "c_quarks_mass", "c_quarks_pt", "c_quarks_rapidity", "c_quarks_pseudorapidity", "c_quarks_phi", "c_quarks_reco_mass", "s_quarks_px", "s_quarks_py", "s_quarks_pz", "s_quarks_energy", "s_quarks_mass", "s_quarks_pt", "s_quarks_rapidity", "s_quarks_pseudorapidity", "s_quarks_phi", "s_quarks_reco_mass", "d_px", "d_py", "d_pz", "d_energy", "d_mass", "d_pt", "d_rapidity", "d_pseudorapidity", "d_phi", "d_reco_mass", "dbar_px", "dbar_py", "dbar_pz", "dbar_energy", "dbar_mass", "dbar_pt", "dbar_rapidity", "dbar_pseudorapidity", "dbar_phi", "dbar_reco_mass", "s_px", "s_py", "s_pz", "s_energy", "s_mass", "s_pt", "s_rapidity", "s_pseudorapidity", "s_phi", "s_reco_mass", "sbar_px", "sbar_py", "sbar_pz", "sbar_energy", "sbar_mass", "sbar_pt", "sbar_rapidity", "sbar_pseudorapidity", "sbar_phi", "sbar_reco_mass", "u_px", "u_py", "u_pz", "u_energy", "u_mass", "u_pt", "u_rapidity", "u_pseudorapidity", "u_phi", "u_reco_mass", "ubar_px", "ubar_py", "ubar_pz", "ubar_energy", "ubar_mass", "ubar_pt", "ubar_rapidity", "ubar_pseudorapidity", "ubar_phi", "ubar_reco_mass", "c_px", "c_py", "c_pz", "c_energy", "c_mass", "c_pt", "c_rapidity", "c_pseudorapidity", "c_phi", "c_reco_mass", "cbar_px", "cbar_py", "cbar_pz", "cbar_energy", "cbar_mass", "cbar_pt", "cbar_rapidity", "cbar_pseudorapidity", "cbar_phi", "cbar_reco_mass", 
        "leptons_px", "leptons_py", "leptons_pz", "leptons_energy", "leptons_mass", "leptons_pt", "leptons_rapidity", "leptons_pseudorapidity", "leptons_phi", "charged_leptons_px", "charged_leptons_py", "charged_leptons_pz", "charged_leptons_energy", "charged_leptons_mass", "charged_leptons_pt", "charged_leptons_rapidity", "charged_leptons_pseudorapidity", "charged_leptons_phi", "neutrinos_px", "neutrinos_py", "neutrinos_pz", "neutrinos_energy", "neutrinos_mass", "neutrinos_pt", "neutrinos_rapidity", "neutrinos_pseudorapidity", "neutrinos_phi", "electron_pairs_px", "electron_pairs_py", "electron_pairs_pz", "electron_pairs_energy", "electron_pairs_mass", "electron_pairs_pt", "electron_pairs_rapidity", "electron_pairs_pseudorapidity", "electron_pairs_phi", "electrons_px", "electrons_py", "electrons_pz", "electrons_energy", "electrons_mass", "electrons_pt", "electrons_rapidity", "electrons_pseudorapidity", "electrons_phi", "positrons_px", "positrons_py", "positrons_pz", "positrons_energy", "positrons_mass", "positrons_pt", "positrons_rapidity", "positrons_pseudorapidity", "positrons_phi", "muon_pairs_px", "muon_pairs_py", "muon_pairs_pz", "muon_pairs_energy", "muon_pairs_mass", "muon_pairs_pt", "muon_pairs_rapidity", "muon_pairs_pseudorapidity", "muon_pairs_phi", "muons_px", "muons_py", "muons_pz", "muons_energy", "muons_mass", "muons_pt", "muons_rapidity", "muons_pseudorapidity", "muons_phi", "antimuons_px", "antimuons_py", "antimuons_pz", "antimuons_energy", "antimuons_mass", "antimuons_pt", "antimuons_rapidity", "antimuons_pseudorapidity", "antimuons_phi", "tau_pairs_px", "tau_pairs_py", "tau_pairs_pz", "tau_pairs_energy", "tau_pairs_mass", "tau_pairs_pt", "tau_pairs_rapidity", "tau_pairs_pseudorapidity", "tau_pairs_phi", "taus_px", "taus_py", "taus_pz", "taus_energy", "taus_mass", "taus_pt", "taus_rapidity", "taus_pseudorapidity", "taus_phi", "antitaus_px", "antitaus_py", "antitaus_pz", "antitaus_energy", "antitaus_mass", "antitaus_pt", "antitaus_rapidity", "antitaus_pseudorapidity", "antitaus_phi", "electron_neutrinos_px", "electron_neutrinos_py", "electron_neutrinos_pz", "electron_neutrinos_energy", "electron_neutrinos_mass", "electron_neutrinos_pt", "electron_neutrinos_rapidity", "electron_neutrinos_pseudorapidity", "electron_neutrinos_phi", "muon_neutrinos_px", "muon_neutrinos_py", "muon_neutrinos_pz", "muon_neutrinos_energy", "muon_neutrinos_mass", "muon_neutrinos_pt", "muon_neutrinos_rapidity", "muon_neutrinos_pseudorapidity", "muon_neutrinos_phi", "tau_neutrinos_px", "tau_neutrinos_py", "tau_neutrinos_pz", "tau_neutrinos_energy", "tau_neutrinos_mass", "tau_neutrinos_pt", "tau_neutrinos_rapidity", "tau_neutrinos_pseudorapidity", "tau_neutrinos_phi",
    ]

    quan_title=[
        "p_{x} [GeV]", "p_{y} [GeV]", "p_{z} [GeV]", "E [GeV]", "m [GeV]", "p_{T} [GeV]", "y", "#eta", "#phi [rad]", "PDG ID", "p_{x} [GeV]", "p_{y} [GeV]", "p_{z} [GeV]", "E [GeV]", "m [GeV]", "p_{T} [GeV]", "y", "#eta", "#phi [rad]", "PDG ID", "p_{x}(l_{j}^{ extra}) [GeV]", "p_{y}(l_{j}^{ extra}) [GeV]", "p_{z}(l_{j}^{ extra}) [GeV]", "E(l_{j}^{ extra}) [GeV]", "m(l_{j}^{ extra}) [GeV]", "p_{T}(l_{j}^{ extra}) [GeV]", "y(l_{j}^{ extra})", "#eta(l_{j}^{ extra})", "#phi(l_{j}^{ extra}) [rad]", "PDG ID", "p_{x}(t) [GeV]", "p_{y}(t) [GeV]", "p_{z}(t) [GeV]", "E(t) [GeV]", "m(t) [GeV]", "p_{T}(t) [GeV]", "y(t)", "#phi(t) [rad]", "p_{x}(#bar{t}) [GeV]", "p_{y}(#bar{t}) [GeV]", "p_{z}(#bar{t}) [GeV]", "E(#bar{t}) [GeV]", "m(#bar{t}) [GeV]", "p_{T}(#bar{t}) [GeV]", "y(#bar{t})", "#phi(#bar{t}) [rad]", "p_{x}(t,#bar{t}) [GeV]", "p_{y}(t,#bar{t}) [GeV]", " p_{z}(t,#bar{t}) [GeV]", "E(t,#bar{t}) [GeV]", "m(t,#bar{t}) [GeV]", "p_{T}(t,#bar{t}) [GeV]", "y(t,#bar{t})", "#phi(t,#bar{t}) [rad]", "p_{x}(W^{+}) [GeV]", "p_{y}(W^{+}) [GeV]", "p_{z}(W^{+}) [GeV]", "E(W^{+}) [GeV]", "m(W^{+}) [GeV]", "p_{T}(W^{+}) [GeV]", "y(W^{+})", "#phi(W^{+}) [rad]", "p_{x}(W^{-}) [GeV]", "p_{y}(W^{-}) [GeV]", "p_{z}(W^{-}) [GeV]", "E(W^{-}) [GeV]", "m(W^{-}) [GeV]", "p_{T}(W^{-}) [GeV]", "y(W^{-})", "#phi(W^{-}) [rad]", "p_{x}(W^{#pm}) [GeV]", "p_{y}(W^{#pm}) [GeV]", "p_{z}(W^{#pm}) [GeV]", "E(W^{#pm}) [GeV]", "m(W^{#pm}) [GeV]", "p_{T}(W^{#pm}) [GeV]", "y(W^{#pm})", "#phi(W^{#pm}) [rad]", "p_{x}(b,#bar{b}^{ all}) [GeV]", "p_{y}(b,#bar{b}^{ all}) [GeV]", "p_{z}(b,#bar{b}^{ all}) [GeV]", "E(b,#bar{b}^{ all}) [GeV]", "m(b,#bar{b}^{ all}) [GeV]", "p_{T}(b,#bar{b}^{ all}) [GeV]", "y(b,#bar{b}^{ all})", "#eta(b,#bar{b}^{ all})", "#phi(b,#bar{b}^{ all}) [rad]", "p_{x}(b^{ all}) [GeV]", "p_{y}(b^{ all}) [GeV]", "p_{z}(b^{ all}) [GeV]", "E(b^{ all}) [GeV]", "m(b^{ all}) [GeV]", "p_{T}(b^{ all}) [GeV]", "y(b^{ all})", "#eta(b^{ all})", "#phi(b^{ all}) [rad]", "p_{x}(#bar{b}^{ all}) [GeV]", "p_{y}(#bar{b}^{ all}) [GeV]", "p_{z}(#bar{b}^{ all}) [GeV]", "E(#bar{b}^{ all}) [GeV]", "m(#bar{b}^{ all}) [GeV]", "p_{T}(#bar{b}^{ all}) [GeV]", "y(#bar{b}^{ all})", "#eta(#bar{b}^{ all})", "#phi(#bar{b}^{ all}) [rad]",  "p_{x}(b,#bar{b}) [GeV]", "p_{y}(b,#bar{b}) [GeV]", "p_{z}(b,#bar{b}) [GeV]", "E(b,#bar{b}) [GeV]", "m(b,#bar{b}) [GeV]", "p_{T}(b,#bar{b}) [GeV]", "y(b,#bar{b})", "#eta(b,#bar{b})", "#phi(b,#bar{b}) [rad]", "p_{x}(b,#bar{b}^{ prompt}) [GeV]", "p_{y}(b,#bar{b}^{ prompt}) [GeV]", "p_{z}(b,#bar{b}^{ prompt}) [GeV]", "E(b,#bar{b}^{ prompt}) [GeV]", "m(b,#bar{b}^{ prompt}) [GeV]", "p_{T}(b,#bar{b}^{ prompt}) [GeV]", "y(b,#bar{b}^{ prompt})", "#eta(b,#bar{b}^{ prompt})", "#phi(b,#bar{b}^{ prompt}) [rad]", "p_{T}(t#bar{t}) [GeV]", "y(t#bar{t})", "H_{T}(t#bar{t}) [GeV]", "m(t#bar{t})", "#DeltaR(t#bar{t})", "#DeltaR2(t#bar{t})", "p_{T}(b#bar{b}^{ all}) [GeV]", "y(b#bar{b}^{ all})","#eta(b#bar{b}^{ all})", "H_{T}(b#bar{b}^{ all}) [GeV]", "m(b#bar{b}^{ all})", "#DeltaR(b#bar{b}^{ all})", "p_{T}(b#bar{b}^{ prompt}) [GeV]", "y(b#bar{b}^{ prompt})", "#eta(b#bar{b}^{ prompt})", "H_{T}(b#bar{b}^{ prompt}) [GeV]", "m(b#bar{b}^{ prompt})", "#DeltaR(b#bar{b}^{ prompt})", "p_{T}(b#bar{b}) [GeV]","y(b#bar{b})", "#eta(b#bar{b})", "H_{T}(b#bar{b}) [GeV]", "m(b#bar{b})", "#DeltaR(b#bar{b})", "m^{reco}(t) [GeV]", "m^{reco}(#bar{t}) [GeV]", "m^{reco}(t,#bar{t}) [GeV]", "m^{reco}(W^{+}) [GeV]", "m^{reco}(W^{-}) [GeV]", "m^{reco}(W^{#pm}) [GeV]", "p_{T}(q'#bar{q}) [GeV]", "y(q'#bar{q})", "#eta(q'#bar{q})", "H_{T}(q'#bar{q}) [GeV]", "m^{reco}(q'#bar{q}) [GeV]", "#DeltaR(q'#bar{q})", "p_{T}(ud) [GeV]", "y(ud)", "#eta(ud)", "H_{T}(ud) [GeV]", "m^{reco}(ud) [GeV]", "#DeltaR(ud)", "p_{T}(cs) [GeV]", "y(cs)", "#eta(cs)", "H_{T}(cs) [GeV]", "m^{reco}(cs) [GeV]", "#DeltaR(cs)", "p_{T}(cd) [GeV]", "y(cd)", "#eta(cd)", "H_{T}(cd) [GeV]", "m^{reco}(cd) [GeV]", "#DeltaR(cd)", "p_{T}(us) [GeV]", "y(us)", "#eta(us)", "H_{T}(us) [GeV]", "m^{reco}(us)", "#DeltaR(us)", "p_{x}(q_{lf}) [GeV]", "p_{y}(q_{lf}) [GeV]", "p_{z}(q_{lf}) [GeV]", "E(q_{lf}) [GeV]", "m(q_{lf}) [GeV]", "p_{T}(q_{lf}) [GeV]", "y(q_{lf})", "#eta(q_{lf})", "#phi(q_{lf}) [rad]", "m^{reco}(q_{lf}) [GeV]", "p_{x}(d, #bar{d}) [GeV]", "p_{y}(d, #bar{d}) [GeV]", "p_{z}(d, #bar{d}) [GeV]", "E(d, #bar{d}) [GeV]", "m(d, #bar{d}) [GeV]", "p_{T}(d, #bar{d}) [GeV]", "y(d, #bar{d})", "#eta(d, #bar{d})", "#phi(d, #bar{d}) [rad]", "m^{reco}(d, #bar{d}) [GeV]", "p_{x}(u, #bar{u}) [GeV]", "p_{y}(u, #bar{u}) [GeV]", "p_{z}(u, #bar{u}) [GeV]", "E(u, #bar{u}) [GeV]", "m(u, #bar{u}) [GeV]", "p_{T}(u, #bar{u}) [GeV]", "y(u, #bar{u})", "#eta(u, #bar{u})", "#phi(u, #bar{u}) [rad]", "m^{reco}(u, #bar{u}) [GeV]", "p_{x}(c, #bar{c}) [GeV]", "p_{y}(c, #bar{c}) [GeV]", "p_{z}(c, #bar{c}) [GeV]", "E(c, #bar{c}) [GeV]", "m(c, #bar{c}) [GeV]", "p_{T}(c, #bar{c}) [GeV]", "y(c, #bar{c})", "#eta(c, #bar{c})", "#phi(c, #bar{c}) [rad]", "m^{reco}(c, #bar{c}) [GeV]", "p_{x}(s, #bar{s}) [GeV]", "p_{y}(s, #bar{s}) [GeV]", "p_{z}(s, #bar{s}) [GeV]", "E(s, #bar{s}) [GeV]", "m(s, #bar{s}) [GeV]", "p_{T}(s, #bar{s}) [GeV]", "y(s, #bar{s})", "#eta(s, #bar{s})", "#phi(s, #bar{s}) [rad]", "m^{reco}(s, #bar{s}) [GeV]", "p_{x}(d) [GeV]", "p_{y}(d) [GeV]", "p_{z}(d) [GeV]", "E(d) [GeV]", "m(d) [GeV]", "p_{T}(d) [GeV]", "y(d)", "#eta(d)", "#phi(d) [rad]", "m^{reco}(d) [GeV]", "p_{x}(#bar{d}) [GeV]", "p_{y}(#bar{d}) [GeV]", "p_{z}(#bar{d}) [GeV]", "E(#bar{d}) [GeV]", "m(#bar{d}) [GeV]", "p_{T}(#bar{d}) [GeV]", "y(#bar{d})", "#eta(#bar{d})", "#phi(#bar{d}) [rad]", "m^{reco}(#bar{d}) [GeV]", "p_{x}(s) [GeV]", "p_{y}(s) [GeV]", "p_{z}(s) [GeV]", "E(s) [GeV]", "m(s) [GeV]", "p_{T}(s) [GeV]", "y(s)", "#eta(s)", "#phi(s) [rad]", "m^{reco}(s) [GeV]", "p_{x}(#bar{s}) [GeV]", "p_{y}(#bar{s}) [GeV]", "p_{z}(#bar{s}) [GeV]", "E(#bar{s}) [GeV]", "m(#bar{s}) [GeV]", "p_{T}(#bar{s}) [GeV]", "y(#bar{s})", "#eta(#bar{s})", "#phi(#bar{s}) [rad]", "m^{reco}(#bar{s}) [GeV]", "p_{x}(u) [GeV]", "p_{y}(u) [GeV]", "p_{z}(u) [GeV]", "E(u) [GeV]", "m(u) [GeV]", "p_{T}(u) [GeV]", "y(u)", "#eta(u)", "#phi(u) [rad]", "m^{reco}(u) [GeV]", "p_{x}(#bar{u}) [GeV]", "p_{y}(#bar{u}) [GeV]", "p_{z}(#bar{u}) [GeV]", "E(#bar{u}) [GeV]", "m(#bar{u}) [GeV]", "p_{T}(#bar{u}) [GeV]", "y(#bar{u})", "#eta(#bar{u})", "#phi(#bar{u}) [rad]", "m^{reco}(#bar{u}) [GeV]", "p_{x}(c) [GeV]", "p_{y}(c) [GeV]", "p_{z}(c) [GeV]", "E(c) [GeV]", "m(c) [GeV]", "p_{T}(c) [GeV]", "y(c)", "#eta(c)", "#phi(c) [rad]", "m^{reco}(c) [GeV]", "p_{x}(#bar{c}) [GeV]", "p_{y}(#bar{c}) [GeV]", "p_{z}(#bar{c}) [GeV]", "E(#bar{c}) [GeV]", "m(#bar{c}) [GeV]", "p_{T}(#bar{c}) [GeV]", "y(#bar{c})", "#eta(#bar{c})", "#phi(#bar{c}) [rad]", "m^{reco}(#bar{c}) [GeV]",
        "p_{x}(l,#nu_{l}) [GeV]", "p_{y}(l,#nu_{l}) [GeV]", "p_{z}(l,#nu_{l}) [GeV]", "E(l,#nu_{l}) [GeV]", "m(l,#nu_{l}) [GeV]", "p_{T}(l,#nu_{l}) [GeV]", "y(l,#nu_{l})", "#eta(l,#nu_{l})", "#phi(l,#nu_{l}) [rad]", "p_{x}(l^{#pm}) [GeV]", "p_{y}(l^{#pm}) [GeV]", "p_{z}(l^{#pm}) [GeV]", "E(l^{#pm}) [GeV]", "m(l^{#pm}) [GeV]", "p_{T}(l^{#pm}) [GeV]", "y(l^{#pm})", "#eta(l^{#pm})", "#phi(l^{#pm}) [rad]", "p_{x}(#nu_{l}) [GeV]", "p_{y}(#nu_{l}) [GeV]", "p_{z}(#nu_{l}) [GeV]", "E(#nu_{l}) [GeV]", "m(#nu_{l}) [GeV]", "p_{T}(#nu_{l}) [GeV]", "y(#nu_{l})", "#eta(#nu_{l})", "#phi(#nu_{l}) [rad]", "p_{x}(e^{#pm}) [GeV]", "p_{y}(e^{#pm}) [GeV]", "p_{z}(e^{#pm}) [GeV]", "E(e^{#pm}) [GeV]", "m(e^{#pm}) [GeV]", "p_{T}(e^{#pm}) [GeV]", "y(e^{#pm})", "#eta(e^{#pm})", "#phi(e^{#pm}) [rad]", "p_{x}(e^{-}) [GeV]", "p_{y}(e^{-}) [GeV]", "p_{z}(e^{-}) [GeV]", "E(e^{-}) [GeV]", "m(e^{-}) [GeV]", "p_{T}(e^{-}) [GeV]", "y(e^{-})", "#eta(e^{-})", "#phi(e^{-}) [rad]", "p_{x}(e^{+}) [GeV]", "p_{y}(e^{+}) [GeV]", "p_{z}(e^{+}) [GeV]", "E(e^{+}) [GeV]", "m(e^{+}) [GeV]", "p_{T}(e^{+}) [GeV]", "y(e^{+})", "#eta(e^{+})", "#phi(e^{+}) [rad]", "p_{x}(#mu^{#pm}) [GeV]", "p_{y}(#mu^{#pm}) [GeV]", "p_{z}(#mu^{#pm}) [GeV]", "E(#mu^{#pm}) [GeV]", "m(#mu^{#pm}) [GeV]", "p_{T}(#mu^{#pm}) [GeV]", "y(#mu^{#pm})", "#eta(#mu^{#pm})", "#phi(#mu^{#pm}) [rad]", "p_{x}(#mu^{-}) [GeV]", "p_{y}(#mu^{-}) [GeV]", "p_{z}(#mu^{-}) [GeV]", "E(#mu^{-}) [GeV]", "m(#mu^{-}) [GeV]", "p_{T}(#mu^{-}) [GeV]", "y(#mu^{-})", "#eta(#mu^{-})", "#phi(#mu^{-}) [rad]", "p_{x}(#mu^{+}) [GeV]", "p_{y}(#mu^{+}) [GeV]", "p_{z}(#mu^{+}) [GeV]", "E(#mu^{+}) [GeV]", "m(#mu^{+}) [GeV]", "p_{T}(#mu^{+}) [GeV]", "y(#mu^{+})", "#eta(#mu^{+})", "#phi(#mu^{+}) [rad]", "p_{x}(#tau^{#pm}) [GeV]", "p_{y}(#tau^{#pm}) [GeV]", "p_{z}(#tau^{#pm}) [GeV]", "E(#tau^{#pm}) [GeV]", "m(#tau^{#pm}) [GeV]", "p_{T}(#tau^{#pm}) [GeV]", "y(#tau^{#pm})", "#eta(#tau^{#pm})", "#phi(#tau^{#pm}) [rad]", "p_{x}(#tau^{-}) [GeV]", "p_{y}(#tau^{-}) [GeV]", "p_{z}(#tau^{-}) [GeV]", "E(#tau^{-}) [GeV]", "m(#tau^{-}) [GeV]", "p_{T}(#tau^{-}) [GeV]", "y(#tau^{-})", "#eta(#tau^{-})", "#phi(#tau^{-}) [rad]", "p_{x}(#tau^{+}) [GeV]", "p_{y}(#tau^{+}) [GeV]", "p_{z}(#tau^{+}) [GeV]", "E(#tau^{+}) [GeV]", "m(#tau^{+}) [GeV]", "p_{T}(#tau^{+}) [GeV]", "y(#tau^{+})", "#eta(#tau^{+})", "#phi(#tau^{+}) [rad]", "p_{x}(#nu_{e}) [GeV]", "p_{y}(#nu_{e}) [GeV]", "p_{z}(#nu_{e}) [GeV]", "E(#nu_{e}) [GeV]", "m(#nu_{e}) [GeV]", "p_{T}(#nu_{e}) [GeV]", "y(#nu_{e})", "#eta(#nu_{e})", "#phi(#nu_{e}) [rad]", "p_{x}(#nu_{#mu}) [GeV]", "p_{y}(#nu_{#mu}) [GeV]", "p_{z}(#nu_{#mu}) [GeV]", "E(#nu_{#mu}) [GeV]", "m(#nu_{#mu}) [GeV]", "p_{T}(#nu_{#mu}) [GeV]", "y(#nu_{#mu})", "#eta(#nu_{#mu})", "#phi(#nu_{#mu}) [rad]", "p_{x}(#nu_{#tau}) [GeV]", "p_{y}(#nu_{#tau}) [GeV]", "p_{z}(#nu_{#tau}) [GeV]", "E(#nu_{#tau}) [GeV]", "m(#nu_{#tau}) [GeV]", "p_{T}(#nu_{#tau}) [GeV]", "y(#nu_{#tau})", "#eta(#nu_{#tau})", "#phi(#nu_{#tau}) [rad]",
    ]
    Nquan=len(quan_name)

    Nfiles = len(root_files)

    cms_color_0 = ROOT.TColor.GetColor(87,144,252) #blue
    cms_color_1 = ROOT.TColor.GetColor(248,156,32) #orange
    cms_color_2 = ROOT.TColor.GetColor(228,37,54) #red
    cms_color_3 = ROOT.TColor.GetColor(150,74,139) #purple
    cms_color_4 = ROOT.TColor.GetColor(156,156,161) #gray
    cms_color_5 = ROOT.TColor.GetColor(122,33,221) #purple


    # Colours
    color_0 = ROOT.TColor.GetColor(165,0,38)
    color_1 = ROOT.TColor.GetColor(215,48,39)
    color_2 = ROOT.TColor.GetColor(244,109,67)
    color_3 = ROOT.TColor.GetColor(253,174,97)
    color_4 = ROOT.TColor.GetColor(254,224,144)
    color_5 = ROOT.TColor.GetColor(255,255,191)
    color_6 = ROOT.TColor.GetColor(224,243,248)
    color_7 = ROOT.TColor.GetColor(171,217,233)
    color_8 = ROOT.TColor.GetColor(116,173,209)
    color_9 = ROOT.TColor.GetColor(69,117,180)
    color_10 = ROOT.TColor.GetColor(49,54,149)

    color_11 = ROOT.TColor.GetColor(55,195,82)
    color_12 = ROOT.TColor.GetColor(128,64,142) 

    color_13 = ROOT.TColor.GetColor(255,128,0)
    color_14 = ROOT.TColor.GetColor(100,149,237)

    file_color = [ROOT.kBlack, cms_color_0, cms_color_1, cms_color_2, cms_color_3, cms_color_4]

    # For legend
    plot_lines = ["#mu_{F}#scale[0.8]{#times}1.0", "#mu_{F}#scale[0.8]{#times}2.0", "#mu_{F}#scale[0.8]{#times}1.0 h#scale[0.65]{dampUP}", "#mu_{F}#scale[0.8]{#times}2.0 h#scale[0.65]{dampUP}", "#mu_{F}#scale[0.8]{#times}1.0 h#scale[0.65]{dampDOWN}",  "#mu_{F}#scale[0.8]{#times}2.0 h#scale[0.65]{dampDOWN}"]
    Nlines=len(plot_lines)

    # Open ROOT files
    file = [ROOT.TFile.Open(root_file, "READ") for root_file in root_files]

    # Load all the histograms from files
    hist = [[ROOT.TH1F() for _ in range(Nquan)] for _ in range(Nfiles)]

    for ifile in range(Nfiles):
        for quan in range(Nquan):
            hist[ifile][quan] = file[ifile].Get(quan_name[quan])
            hist[ifile][quan].SetDirectory(0)

    # Make plots and save
    output_file = ROOT.TFile.Open(output_directory + "/output.root", "RECREATE")

    # Define canvas
    canvas_ratio = [ROOT.TCanvas() for _ in range(Nquan)]


    for quan in range(Nquan):
        #Triple ratio
        canvas_ratio[quan] = ROOT.TCanvas("canvas_ratio_" + quan_name[quan], quan_name[quan] + " Ratio", 800, 600)
        
        #xlow, ylow, xup, yup
        pad1 = ROOT.TPad("pad1_" + quan_name[quan], quan_name[quan] + " pad1", 0, 0.29, 1, 1)
        pad1.SetBottomMargin(0.04)  # Set bottom margin for pad1
        # pad1.SetLogy()  # Set log scale for pad1
        pad1.SetFrameLineWidth(2)
        # pad1.SetGridx()
        # pad1.SetGridy()
        pad1.Draw()
        pad1.cd()

        # Make histograms
        for ifile in range(Nfiles):
            hist[ifile][quan].Scale(1.0/hist[ifile][quan].Integral())
            hist[ifile][quan].SetLineColor(file_color[ifile])
            #hist[ifile][quan].SetLineColor(ROOT.kBlack)
            hist[ifile][quan].SetLineWidth(2)
            hist[ifile][quan].SetMarkerColor(file_color[ifile])
            hist[ifile][quan].SetMarkerStyle(7)
            hist[ifile][quan].SetStats(0)
        
        hist[1][quan].SetLineStyle(1)
        hist[2][quan].SetLineStyle(2)
        hist[3][quan].SetLineStyle(2)
        hist[4][quan].SetLineStyle(7)
        hist[5][quan].SetLineStyle(7)

            #hist[ifile][quan].SetAxisRange(axes_limits[quan][0], axes_limits[quan][1], "X")

        # hist[0][quan].GetXaxis().SetTitle(quan_title[quan])
        hist[0][quan].GetYaxis().SetTitle("Entries (Normalized to 1)")
        hist[0][quan].GetYaxis().SetTitleSize(0.045)
        hist[0][quan].GetYaxis().SetLabelSize(0.045)
        hist[0][quan].GetXaxis().SetTitle("")
        hist[0][quan].GetXaxis().SetLabelSize(0)
        
        hist[0][quan].Draw("hist")
        for ifile in range(1, Nfiles):
            hist[ifile][quan].Draw("SAME hist")

        # Create a legend 
        legend = ROOT.TLegend(0.66, 0.62, 0.88, 0.88)
        legend.SetBorderSize(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.045)
        # legend.SetHeader(quan_title[quan], "C")
        for pline in range(Nlines):
            entry = legend.AddEntry(hist[pline][quan], plot_lines[pline], "l")
            entry.SetTextColor(file_color[pline])
        legend.Draw()

        latex = ROOT.TLatex()
        latex.SetTextSize(0.0585)
        latex.DrawLatexNDC(0.1, 0.91, "#font[61]{CMS}")
        latex.SetTextSize(0.045)
        latex.DrawLatexNDC(0.167, 0.91, "#font[52]{Simulation Work-in-Progress}")
        latex.SetTextSize(0.045)
        latex.DrawLatexNDC(0.79, 0.91, "#font[42]{SL channel}")


        # Update pad1
        pad1.Update()

        # Go back to the main canvas
        canvas_ratio[quan].cd()

        # Create a new pad for the histogram above the ratio plot
        pad2 = ROOT.TPad("pad2_" + quan_name[quan], quan_name[quan] + " pad2", 0, 0, 1, 0.29)
        pad2.SetTopMargin(0.03)  # Set top margin for pad2
        pad2.SetBottomMargin(0.3)  # Set bottom margin for pad2
        pad2.SetGridy()  # Add horizontal grid lines to pad2
        pad2.SetFrameLineWidth(2)
        pad2.Draw()
        pad2.cd()

        ratio1 = hist[0][quan].Clone()
        ratio1.Divide(hist[1][quan])

        ratio2 = hist[0][quan].Clone()
        ratio2.Divide(hist[2][quan])

        ratio3 = hist[0][quan].Clone()
        ratio3.Divide(hist[3][quan])

        ratio4 = hist[0][quan].Clone()
        ratio4.Divide(hist[4][quan])

        ratio5 = hist[0][quan].Clone()
        ratio5.Divide(hist[5][quan])

        ratio1Error = ratio1.Clone()
        ratio2Error = ratio2.Clone()
        ratio3Error = ratio3.Clone()
        ratio4Error = ratio4.Clone()
        ratio5Error = ratio5.Clone()


        ratio1.SetLineColor(cms_color_0)
        ratio1.SetMarkerStyle(2)
        ratio1.SetMarkerColor(cms_color_0)
        ratio1.SetMarkerSize(0)
        ratio1.SetLineWidth(2)
        ratio1.SetLineStyle(1)
        ratio1Error.SetFillColorAlpha(cms_color_0, 0.2)
        ratio1Error.SetMarkerSize(0)
        #ratio1Error.SetFillStyle(3003)

        ratio2.SetLineColor(cms_color_1)
        ratio2.SetMarkerStyle(2)
        ratio2.SetMarkerColor(cms_color_1)
        ratio2.SetMarkerSize(0)
        ratio2.SetLineWidth(2)
        ratio2.SetLineStyle(2)
        ratio2Error.SetFillColorAlpha(cms_color_1, 0.2)
        ratio2Error.SetMarkerSize(0)
        #ratio2Error.SetFillStyle(3007)

        ratio3.SetLineColor(cms_color_2)
        ratio3.SetMarkerStyle(2)
        ratio3.SetMarkerColor(cms_color_2)
        ratio3.SetMarkerSize(0)
        ratio3.SetLineWidth(2)
        ratio3.SetLineStyle(2)
        ratio3Error.SetFillColorAlpha(cms_color_2, 0.2)
        ratio3Error.SetMarkerSize(0)

        ratio4.SetLineColor(cms_color_3)
        ratio4.SetMarkerStyle(2)
        ratio4.SetMarkerColor(cms_color_3)
        ratio4.SetMarkerSize(0)
        ratio4.SetLineWidth(2)
        ratio4.SetLineStyle(7)
        ratio4Error.SetFillColorAlpha(cms_color_3, 0.2)
        ratio4Error.SetMarkerSize(0)
        #ratio4Error.SetFillStyle(3003)

        ratio5.SetLineColor(cms_color_4)
        ratio5.SetMarkerStyle(2)
        ratio5.SetMarkerColor(cms_color_4)
        ratio5.SetMarkerSize(0)
        ratio5.SetLineWidth(2)
        ratio5.SetLineStyle(7)
        ratio5Error.SetFillColorAlpha(cms_color_4, 0.2)
        ratio5Error.SetMarkerSize(0)
        #ratio5Error.SetFillStyle(3007)

        ratio1.GetYaxis().SetLabelSize(0.09)
        ratio1.GetYaxis().SetTitle("")
        ratio1.GetYaxis().SetTitleSize(0.06)
        ratio1.GetXaxis().SetLabelSize(0.09)
        ratio1.GetXaxis().SetTitle(quan_title[quan])
        ratio1.GetXaxis().SetTitleSize(0.12)
        
        ratio1.SetAxisRange(0.8, 1.2, "Y")
        ratio1.GetYaxis().SetNdivisions(505)

        ratio1.Draw("hist 9")
        ratio2.Draw("hist same 9")
        ratio3.Draw("hist same 9")
        ratio4.Draw("hist same 9")
        ratio5.Draw("hist same 9")
        ratio1Error.Draw("E2 same 9")
        ratio2Error.Draw("E2 same 9")
        ratio3Error.Draw("E2 same 9")
        ratio4Error.Draw("E2 same 9")
        ratio5Error.Draw("E2 same 9")


        # Draw a horizontal line at y=1 for reference
        line = ROOT.TLine(hist[0][quan].GetXaxis().GetXmin(), 1, hist[0][quan].GetXaxis().GetXmax(), 1)
        line.SetLineStyle(2)
        line.Draw()

        latex2 = ROOT.TLatex()
        latex2.SetTextSize(0.095)
        latex2.SetTextAngle(90)
        latex2.DrawLatexNDC(0.04, 0.24, "#font[42]{Nominal/Variations}")

        pad2.Update()

        # Update the canvas
        canvas_ratio[quan].Update()

        canvas_ratio[quan].SaveAs(output_directory + "/ratio_" + quan_name[quan] + ".png")
        canvas_ratio[quan].SaveAs(output_directory + "/ratio_" + quan_name[quan] + ".pdf")
        canvas_ratio[quan].Write()

    # Close all files
    output_file.Close()

    for ifile in range(Nfiles):
        file[ifile].Close()

if __name__ == "__main__":
    root_files = ["lhe_production_SL_r1.0_f1.0_m172.5_p320900_plots_linear.root", "lhe_production_SL_r1.0_f2.0_m172.5_p320900_plots_linear.root", "lhe_hdampUP_SL_r1.0_f1.0_m172.5_p320900_plots_linear.root", "lhe_hdampUP_SL_r1.0_f2.0_m172.5_p320900_plots_linear.root", "lhe_hdampDOWN_SL_r1.0_f1.0_m172.5_p320900_plots_linear.root",  "lhe_hdampDOWN_SL_r1.0_f2.0_m172.5_p320900_plots_linear.root"]
    output_overlay_directory = "ratio_all_linear_SL"

    if not os.path.exists(output_overlay_directory):
        os.makedirs(output_overlay_directory)

    overlay_and_ratio_plot(root_files, output_overlay_directory)
