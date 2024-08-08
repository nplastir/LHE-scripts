import ROOT
import os

def overlay_and_ratio_plot(root_files, output_directory):

    ROOT.gStyle.SetOptTitle(0) #removes the hist title
    # Include axes on all sides 
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    # Number of quantities
    quan_name=[
        "all_px", "all_py", "all_pz", "all_energy", "all_mass", "all_pt", "all_rapidity", "all_pseudorapidity", "all_phi", "all_pdgid", "initial_particles_px", "initial_particles_py", "initial_particles_pz", "initial_particles_energy", "initial_particles_mass", "initial_particles_pt", "initial_particles_rapidity", "initial_particles_pseudorapidity", "initial_particles_phi", "initial_particles_pdgid", "radiation_px", "radiation_py", "radiation_pz", "radiation_energy", "radiation_mass", "radiation_pt", "radiation_rapidity", "radiation_pseudorapidity", "radiation_phi", "radiation_pdgid", "top_px", "top_py", "top_pz", "top_energy", "top_mass", "top_pt", "top_rapidity", "top_phi", "antitop_px", "antitop_py", "antitop_pz", "antitop_energy", "antitop_mass", "antitop_pt", "antitop_rapidity", "antitop_phi", "ttbar_px", "ttbar_py", "ttbar_pz", "ttbar_energy", "ttbar_mass", "ttbar_pt", "ttbar_rapidity", "ttbar_phi", "Wp_px", "Wp_py", "Wp_pz", "Wp_energy", "Wp_mass", "Wp_pt", "Wp_rapidity", "Wp_phi", "Wm_px", "Wm_py", "Wm_pz", "Wm_energy", "Wm_mass", "Wm_pt", "Wm_rapidity", "Wm_phi", "W_px", "W_py", "W_pz", "W_energy", "W_mass", "W_pt", "W_rapidity", "W_phi", "b_all_px", "b_all_py", "b_all_pz", "b_all_energy", "b_all_mass", "b_all_pt", "b_all_rapidity", "b_all_pseudorapidity", "b_all_phi", "bottom_px", "bottom_py", "bottom_pz", "bottom_energy", "bottom_mass", "bottom_pt", "bottom_rapidity", "bottom_pseudorapidity", "bottom_phi", "antibottom_px", "antibottom_py", "antibottom_pz", "antibottom_energy", "antibottom_mass", "antibottom_pt", "antibottom_rapidity", "antibottom_pseudorapidity", "antibottom_phi", "b_from_top_px", "b_from_top_py", "b_from_top_pz", "b_from_top_energy", "b_from_top_mass", "b_from_top_pt", "b_from_top_rapidity", "b_from_top_pseudorapidity", "b_from_top_phi", "b_from_initial_px", "b_from_initial_py", "b_from_initial_pz", "b_from_initial_energy", "b_from_initial_mass", "b_from_initial_pt", "b_from_initial_rapidity", "b_from_initial_pseudorapidity", "b_from_initial_phi", "tt_system_pt", "tt_system_rapidity", "tt_system_HT", "tt_system_invariant_mass", "tt_system_dR", "tt_system_dR_2", "bb_system_pt", "bb_system_rapidity", "bb_system_pseudorapidity", "bb_system_HT", "bb_system_invariant_mass", "bb_system_dR", "prompt_bs_pt", "prompt_bs_rapidity","prompt_bs_pseudorapidity", "prompt_bs_HT", "prompt_bs_invariant_mass", "prompt_bs_dR", "bs_from_top_pt", "bs_from_top_rapidity", "bs_from_top_pseudorapidity", "bs_from_top_HT", "bs_from_top_invariant_mass", "bs_from_top_dR", "top_reco_mass", "antitop_reco_mass", "ttbar_reco_mass", "Wp_reco_mass", "Wm_reco_mass", "W_reco_mass", "quarks_system_pt", "quarks_system_rapidity", "quarks_system_pseudorapidity", "quarks_system_HT", "quarks_system_invariant_mass", "quarks_system_dR","dd_system_pt", "dd_system_rapidity", "dd_system_pseudorapidity", "dd_system_HT", "dd_system_invariant_mass", "dd_system_dR", "uu_system_pt", "uu_system_rapidity", "uu_system_pseudorapidity", "uu_system_HT", "uu_system_invariant_mass", "uu_system_dR", "ss_system_pt", "ss_system_rapidity", "ss_system_pseudorapidity", "ss_system_HT", "ss_system_invariant_mass", "ss_system_dR", "cc_system_pt", "cc_system_rapidity", "cc_system_pseudorapidity", "cc_system_HT", "cc_system_invariant_mass", "cc_system_dR", "ud_pairs_system_pt", "ud_pairs_system_rapidity", "ud_pairs_system_pseudorapidity", "ud_pairs_system_HT", "ud_pairs_system_invariant_mass", "ud_pairs_system_dR", "sc_pairs_system_pt", "sc_pairs_system_rapidity", "sc_pairs_system_pseudorapidity", "sc_pairs_system_HT", "sc_pairs_system_invariant_mass", "sc_pairs_system_dR", "cd_pairs_system_pt", "cd_pairs_system_rapidity", "cd_pairs_system_pseudorapidity", "cd_pairs_system_HT", "cd_pairs_system_invariant_mass", "cd_pairs_system_dR", "su_pairs_system_pt", "su_pairs_system_rapidity", "su_pairs_system_pseudorapidity", "su_pairs_system_HT", "su_pairs_system_invariant_mass", "su_pairs_system_dR", "dubar_pairs_system_pt", "dubar_pairs_system_rapidity", "dubar_pairs_system_pseudorapidity", "dubar_pairs_system_HT", "dubar_pairs_system_invariant_mass", "dubar_pairs_system_dR", "dbaru_pairs_system_pt", "dbaru_pairs_system_rapidity", "dbaru_pairs_system_pseudorapidity", "dbaru_pairs_system_HT", "dbaru_pairs_system_invariant_mass", "dbaru_pairs_system_dR", "scbar_pairs_system_pt", "scbar_pairs_system_rapidity", "scbar_pairs_system_pseudorapidity", "scbar_pairs_system_HT", "scbar_pairs_system_invariant_mass", "scbar_pairs_system_dR", "sbarc_pairs_system_pt", "sbarc_pairs_system_rapidity", "sbarc_pairs_system_pseudorapidity", "sbarc_pairs_system_HT", "sbarc_pairs_system_invariant_mass", "sbarc_pairs_system_dR", "subar_pairs_system_pt", "subar_pairs_system_rapidity", "subar_pairs_system_pseudorapidity", "subar_pairs_system_HT", "subar_pairs_system_invariant_mass", "subar_pairs_system_dR", "sbaru_pairs_system_pt", "sbaru_pairs_system_rapidity", "sbaru_pairs_system_pseudorapidity", "sbaru_pairs_system_HT", "sbaru_pairs_system_invariant_mass", "dcbar_pairs_system_pt", "dcbar_pairs_system_rapidity", "dcbar_pairs_system_pseudorapidity", "dcbar_pairs_system_HT", "dcbar_pairs_system_invariant_mass", "dcbar_pairs_system_dR", "sbaru_pairs_system_dR", "dbarc_pairs_system_pt", "dbarc_pairs_system_rapidity", "dbarc_pairs_system_pseudorapidity", "dbarc_pairs_system_HT", "dbarc_pairs_system_invariant_mass", "dbarc_pairs_system_dR", "pairs_from_Wm_pt", "pairs_from_Wm_rapidity", "pairs_from_Wm_pseudorapidity", "pairs_from_Wm_HT", "pairs_from_Wm_invariant_mass", "pairs_from_Wm_dR", "pairs_from_Wp_pt", "pairs_from_Wp_rapidity", "pairs_from_Wp_pseudorapidity", "pairs_from_Wp_HT", "pairs_from_Wp_invariant_mass", "pairs_from_Wp_dR", "quarks_px", "quarks_py", "quarks_pz", "quarks_energy", "quarks_mass", "quarks_pt", "quarks_rapidity", "quarks_pseudorapidity", "quarks_phi", "quarks_reco_mass", "d_quarks_px", "d_quarks_py", "d_quarks_pz", "d_quarks_energy", "d_quarks_mass", "d_quarks_pt", "d_quarks_rapidity", "d_quarks_pseudorapidity", "d_quarks_phi", "d_quarks_reco_mass", "u_quarks_px", "u_quarks_py", "u_quarks_pz", "u_quarks_energy", "u_quarks_mass", "u_quarks_pt", "u_quarks_rapidity", "u_quarks_pseudorapidity", "u_quarks_phi", "u_quarks_reco_mass", "c_quarks_px", "c_quarks_py", "c_quarks_pz", "c_quarks_energy", "c_quarks_mass", "c_quarks_pt", "c_quarks_rapidity", "c_quarks_pseudorapidity", "c_quarks_phi", "c_quarks_reco_mass", "s_quarks_px", "s_quarks_py", "s_quarks_pz", "s_quarks_energy", "s_quarks_mass", "s_quarks_pt", "s_quarks_rapidity", "s_quarks_pseudorapidity", "s_quarks_phi", "s_quarks_reco_mass", "d_px", "d_py", "d_pz", "d_energy", "d_mass", "d_pt", "d_rapidity", "d_pseudorapidity", "d_phi", "d_reco_mass", "dbar_px", "dbar_py", "dbar_pz", "dbar_energy", "dbar_mass", "dbar_pt", "dbar_rapidity", "dbar_pseudorapidity", "dbar_phi", "dbar_reco_mass", "s_px", "s_py", "s_pz", "s_energy", "s_mass", "s_pt", "s_rapidity", "s_pseudorapidity", "s_phi", "s_reco_mass", "sbar_px", "sbar_py", "sbar_pz", "sbar_energy", "sbar_mass", "sbar_pt", "sbar_rapidity", "sbar_pseudorapidity", "sbar_phi", "sbar_reco_mass", "u_px", "u_py", "u_pz", "u_energy", "u_mass", "u_pt", "u_rapidity", "u_pseudorapidity", "u_phi", "u_reco_mass", "ubar_px", "ubar_py", "ubar_pz", "ubar_energy", "ubar_mass", "ubar_pt", "ubar_rapidity", "ubar_pseudorapidity", "ubar_phi", "ubar_reco_mass", "c_px", "c_py", "c_pz", "c_energy", "c_mass", "c_pt", "c_rapidity", "c_pseudorapidity", "c_phi", "c_reco_mass", "cbar_px", "cbar_py", "cbar_pz", "cbar_energy", "cbar_mass", "cbar_pt", "cbar_rapidity", "cbar_pseudorapidity", "cbar_phi", "cbar_reco_mass",
    ]

    quan_title=[
        "p_{x} [GeV]", "p_{y} [GeV]", "p_{z} [GeV]", "E [GeV]", "m [GeV]", "p_{T} [GeV]", "y", "#eta", "#phi [rad]", "PDG ID", "p_{x} [GeV]", "p_{y} [GeV]", "p_{z} [GeV]", "E [GeV]", "m [GeV]", "p_{T} [GeV]", "y", "#eta", "#phi [rad]", "PDG ID", "p_{x}(l_{j}^{ extra}) [GeV]", "p_{y}(l_{j}^{ extra}) [GeV]", "p_{z}(l_{j}^{ extra}) [GeV]", "E(l_{j}^{ extra}) [GeV]", "m(l_{j}^{ extra}) [GeV]", "p_{T}(l_{j}^{ extra}) [GeV]", "y(l_{j}^{ extra})", "#eta(l_{j}^{ extra})", "#phi(l_{j}^{ extra}) [rad]", "PDG ID", "p_{x}(t) [GeV]", "p_{y}(t) [GeV]", "p_{z}(t) [GeV]", "E(t) [GeV]", "m(t) [GeV]", "p_{T}(t) [GeV]", "y(t)", "#phi(t) [rad]", "p_{x}(#bar{t}) [GeV]", "p_{y}(#bar{t}) [GeV]", "p_{z}(#bar{t}) [GeV]", "E(#bar{t}) [GeV]", "m(#bar{t}) [GeV]", "p_{T}(#bar{t}) [GeV]", "y(#bar{t})", "#phi(#bar{t}) [rad]", "p_{x}(t,#bar{t}) [GeV]", "p_{y}(t,#bar{t}) [GeV]", " p_{z}(t,#bar{t}) [GeV]", "E(t,#bar{t}) [GeV]", "m(t,#bar{t}) [GeV]", "p_{T}(t,#bar{t}) [GeV]", "y(t,#bar{t})", "#phi(t,#bar{t}) [rad]", "p_{x}(W^{+}) [GeV]", "p_{y}(W^{+}) [GeV]", "p_{z}(W^{+}) [GeV]", "E(W^{+}) [GeV]", "m(W^{+}) [GeV]", "p_{T}(W^{+}) [GeV]", "y(W^{+})", "#phi(W^{+}) [rad]", "p_{x}(W^{-}) [GeV]", "p_{y}(W^{-}) [GeV]", "p_{z}(W^{-}) [GeV]", "E(W^{-}) [GeV]", "m(W^{-}) [GeV]", "p_{T}(W^{-}) [GeV]", "y(W^{-})", "#phi(W^{-}) [rad]", "p_{x}(W^{#pm}) [GeV]", "p_{y}(W^{#pm}) [GeV]", "p_{z}(W^{#pm}) [GeV]", "E(W^{#pm}) [GeV]", "m(W^{#pm}) [GeV]", "p_{T}(W^{#pm}) [GeV]", "y(W^{#pm})", "#phi(W^{#pm}) [rad]", "p_{x}(b,#bar{b}^{ all}) [GeV]", "p_{y}(b,#bar{b}^{ all}) [GeV]", "p_{z}(b,#bar{b}^{ all}) [GeV]", "E(b,#bar{b}^{ all}) [GeV]", "m(b,#bar{b}^{ all}) [GeV]", "p_{T}(b,#bar{b}^{ all}) [GeV]", "y(b,#bar{b}^{ all})", "#eta(b,#bar{b}^{ all})", "#phi(b,#bar{b}^{ all}) [rad]", "p_{x}(b^{ all}) [GeV]", "p_{y}(b^{ all}) [GeV]", "p_{z}(b^{ all}) [GeV]", "E(b^{ all}) [GeV]", "m(b^{ all}) [GeV]", "p_{T}(b^{ all}) [GeV]", "y(b^{ all})", "#eta(b^{ all})", "#phi(b^{ all}) [rad]", "p_{x}(#bar{b}^{ all}) [GeV]", "p_{y}(#bar{b}^{ all}) [GeV]", "p_{z}(#bar{b}^{ all}) [GeV]", "E(#bar{b}^{ all}) [GeV]", "m(#bar{b}^{ all}) [GeV]", "p_{T}(#bar{b}^{ all}) [GeV]", "y(#bar{b}^{ all})", "#eta(#bar{b}^{ all})", "#phi(#bar{b}^{ all}) [rad]",  "p_{x}(b,#bar{b}) [GeV]", "p_{y}(b,#bar{b}) [GeV]", "p_{z}(b,#bar{b}) [GeV]", "E(b,#bar{b}) [GeV]", "m(b,#bar{b}) [GeV]", "p_{T}(b,#bar{b}) [GeV]", "y(b,#bar{b})", "#eta(b,#bar{b})", "#phi(b,#bar{b}) [rad]", "p_{x}(b,#bar{b}^{ prompt}) [GeV]", "p_{y}(b,#bar{b}^{ prompt}) [GeV]", "p_{z}(b,#bar{b}^{ prompt}) [GeV]", "E(b,#bar{b}^{ prompt}) [GeV]", "m(b,#bar{b}^{ prompt}) [GeV]", "p_{T}(b,#bar{b}^{ prompt}) [GeV]", "y(b,#bar{b}^{ prompt})", "#eta(b,#bar{b}^{ prompt})", "#phi(b,#bar{b}^{ prompt}) [rad]", "p_{T}(t#bar{t}) [GeV]", "y(t#bar{t})", "H_{T}(t#bar{t}) [GeV]", "m(t#bar{t})", "#DeltaR(t#bar{t})", "#DeltaR2(t#bar{t})", "p_{T}(b#bar{b}^{ all}) [GeV]", "y(b#bar{b}^{ all})","#eta(b#bar{b}^{ all})", "H_{T}(b#bar{b}^{ all}) [GeV]", "m(b#bar{b}^{ all})", "#DeltaR(b#bar{b}^{ all})", "p_{T}(b#bar{b}^{ prompt}) [GeV]", "y(b#bar{b}^{ prompt})", "#eta(b#bar{b}^{ prompt})", "H_{T}(b#bar{b}^{ prompt}) [GeV]", "m(b#bar{b}^{ prompt})", "#DeltaR(b#bar{b}^{ prompt})", "p_{T}(b#bar{b}) [GeV]","y(b#bar{b})", "#eta(b#bar{b})", "H_{T}(b#bar{b}) [GeV]", "m(b#bar{b})", "#DeltaR(b#bar{b})", "m^{reco}(t) [GeV]", "m^{reco}(#bar{t}) [GeV]", "m^{reco}(t,#bar{t}) [GeV]", "m^{reco}(W^{+}) [GeV]", "m^{reco}(W^{-}) [GeV]", "m^{reco}(W^{#pm}) [GeV]", "p_{T}(q'#bar{q}) [GeV]", "y(q'#bar{q})", "#eta(q'#bar{q})", "H_{T}(q'#bar{q}) [GeV]", "m^{reco}(q'#bar{q}) [GeV]", "#DeltaR(q'#bar{q})", "p_{T}(d#bar{d}) [GeV]", "y(d#bar{d})", "#eta(d#bar{d})", "H_{T}(d#bar{d}) [GeV]", "m^{reco}(d#bar{d}) [GeV]", "#DeltaR(d#bar{d})", "p_{T}(u#bar{u}) [GeV]", "y(u#bar{u})", "#eta(u#bar{u})", "H_{T}(u#bar{u}) [GeV]", "m^{reco}(u#bar{u}) [GeV]", "#DeltaR(u#bar{u})", "p_{T}(s#bar{s}) [GeV]", "y(s#bar{s})", "#eta(s#bar{s})", "H_{T}(s#bar{s}) [GeV]", "m^{reco}(s#bar{s}) [GeV]", "#DeltaR(s#bar{s})", "p_{T}(c#bar{c}) [GeV]", "y(c#bar{c})", "#eta(c#bar{c})", "H_{T}(c#bar{c}) [GeV]", "m^{reco}(c#bar{c}) [GeV]", "#DeltaR(c#bar{c})", "p_{T}(ud) [GeV]", "y(ud)", "#eta(ud)", "H_{T}(ud) [GeV]", "m^{reco}(ud) [GeV]", "#DeltaR(ud)", "p_{T}(cs) [GeV]", "y(cs)", "#eta(cs)", "H_{T}(cs) [GeV]", "m^{reco}(cs) [GeV]", "#DeltaR(cs)", "p_{T}(cd) [GeV]", "y(cd)", "#eta(cd)", "H_{T}(cd) [GeV]", "m^{reco}(cd) [GeV]", "#DeltaR(cd)", "p_{T}(us) [GeV]", "y(us)", "#eta(us)", "H_{T}(us) [GeV]", "m^{reco}(us)", "#DeltaR(us)", "p_{T}(#bar{u}d) [GeV]", "y(#bar{u}d)", "#eta(#bar{u}d)", "H_{T}(#bar{u}d) [GeV]", "m^{reco}(#bar{u}d) [GeV]", "#DeltaR(#bar{u}d)", "p_{T}(u#bar{d}) [GeV]", "y(u#bar{d})", "#eta(u#bar{d})", "H_{T}(u#bar{d}) [GeV]", "m^{reco}(u#bar{d}) [GeV]", "#DeltaR(u#bar{d})", "p_{T}(#bar{c}s) [GeV]", "y(#bar{c}s)", "#eta(#bar{c}s)", "H_{T}(#bar{c}s) [GeV]", "m^{reco}(#bar{c}s) [GeV]", "#DeltaR(#bar{c}s)", "p_{T}(c#bar{s}) [GeV]", "y(c#bar{s})", "#eta(c#bar{s})", "H_{T}(c#bar{s}) [GeV]", "m^{reco}(c#bar{s}) [GeV]", "#DeltaR(c#bar{s})", "p_{T}(#bar{u}s) [GeV]", "y(#bar{u}s)", "#eta(#bar{u}s)", "H_{T}(#bar{u}s) [GeV]", "m^{reco}(#bar{u}s) [GeV]", "#DeltaR(#bar{u}s)", "p_{T}(u#bar{s}) [GeV]", "y(u#bar{s})", "#eta(u#bar{s})", "H_{T}(u#bar{s}) [GeV]", "m^{reco}(u#bar{s}) [GeV]", "#DeltaR(u#bar{s})", "p_{T}(#bar{c}d) [GeV]", "y(#bar{c}d)", "#eta(#bar{c}d)", "H_{T}(#bar{c}d) [GeV]", "m^{reco}(#bar{c}d) [GeV]", "#DeltaR(#bar{c}d)", "p_{T}(c#bar{d}) [GeV]", "y(c#bar{d})", "#eta(c#bar{d})", "H_{T}(c#bar{d}) [GeV]", "m^{reco}(c#bar{d}) [GeV]", "#DeltaR(c#bar{d})", "p_{T}(W^{-}_{reco}) [GeV]", "y(W^{-}_{reco})", "#eta(W^{-}_{reco})", "H_{T}(W^{-}_{reco}) [GeV]", "m^{reco}(W^{-}_{reco}) [GeV]", "#DeltaR(W^{-}_{reco})", "p_{T}(W^{+}_{reco}) [GeV]", "y(W^{+}_{reco})", "#eta(W^{+}_{reco})", "H_{T}(W^{+}_{reco}) [GeV]", "m^{reco}(W^{+}_{reco}) [GeV]", "#DeltaR(W^{+}_{reco})", "p_{x}(q_{lf}) [GeV]", "p_{y}(q_{lf}) [GeV]", "p_{z}(q_{lf}) [GeV]", "E(q_{lf}) [GeV]", "m(q_{lf}) [GeV]", "p_{T}(q_{lf}) [GeV]", "y(q_{lf})", "#eta(q_{lf})", "#phi(q_{lf}) [rad]", "m^{reco}(q_{lf}) [GeV]", "p_{x}(d, #bar{d}) [GeV]", "p_{y}(d, #bar{d}) [GeV]", "p_{z}(d, #bar{d}) [GeV]", "E(d, #bar{d}) [GeV]", "m(d, #bar{d}) [GeV]", "p_{T}(d, #bar{d}) [GeV]", "y(d, #bar{d})", "#eta(d, #bar{d})", "#phi(d, #bar{d}) [rad]", "m^{reco}(d, #bar{d}) [GeV]", "p_{x}(u, #bar{u}) [GeV]", "p_{y}(u, #bar{u}) [GeV]", "p_{z}(u, #bar{u}) [GeV]", "E(u, #bar{u}) [GeV]", "m(u, #bar{u}) [GeV]", "p_{T}(u, #bar{u}) [GeV]", "y(u, #bar{u})", "#eta(u, #bar{u})", "#phi(u, #bar{u}) [rad]", "m^{reco}(u, #bar{u}) [GeV]", "p_{x}(c, #bar{c}) [GeV]", "p_{y}(c, #bar{c}) [GeV]", "p_{z}(c, #bar{c}) [GeV]", "E(c, #bar{c}) [GeV]", "m(c, #bar{c}) [GeV]", "p_{T}(c, #bar{c}) [GeV]", "y(c, #bar{c})", "#eta(c, #bar{c})", "#phi(c, #bar{c}) [rad]", "m^{reco}(c, #bar{c}) [GeV]", "p_{x}(s, #bar{s}) [GeV]", "p_{y}(s, #bar{s}) [GeV]", "p_{z}(s, #bar{s}) [GeV]", "E(s, #bar{s}) [GeV]", "m(s, #bar{s}) [GeV]", "p_{T}(s, #bar{s}) [GeV]", "y(s, #bar{s})", "#eta(s, #bar{s})", "#phi(s, #bar{s}) [rad]", "m^{reco}(s, #bar{s}) [GeV]", "p_{x}(d) [GeV]", "p_{y}(d) [GeV]", "p_{z}(d) [GeV]", "E(d) [GeV]", "m(d) [GeV]", "p_{T}(d) [GeV]", "y(d)", "#eta(d)", "#phi(d) [rad]", "m^{reco}(d) [GeV]", "p_{x}(#bar{d}) [GeV]", "p_{y}(#bar{d}) [GeV]", "p_{z}(#bar{d}) [GeV]", "E(#bar{d}) [GeV]", "m(#bar{d}) [GeV]", "p_{T}(#bar{d}) [GeV]", "y(#bar{d})", "#eta(#bar{d})", "#phi(#bar{d}) [rad]", "m^{reco}(#bar{d}) [GeV]", "p_{x}(s) [GeV]", "p_{y}(s) [GeV]", "p_{z}(s) [GeV]", "E(s) [GeV]", "m(s) [GeV]", "p_{T}(s) [GeV]", "y(s)", "#eta(s)", "#phi(s) [rad]", "m^{reco}(s) [GeV]", "p_{x}(#bar{s}) [GeV]", "p_{y}(#bar{s}) [GeV]", "p_{z}(#bar{s}) [GeV]", "E(#bar{s}) [GeV]", "m(#bar{s}) [GeV]", "p_{T}(#bar{s}) [GeV]", "y(#bar{s})", "#eta(#bar{s})", "#phi(#bar{s}) [rad]", "m^{reco}(#bar{s}) [GeV]", "p_{x}(u) [GeV]", "p_{y}(u) [GeV]", "p_{z}(u) [GeV]", "E(u) [GeV]", "m(u) [GeV]", "p_{T}(u) [GeV]", "y(u)", "#eta(u)", "#phi(u) [rad]", "m^{reco}(u) [GeV]", "p_{x}(#bar{u}) [GeV]", "p_{y}(#bar{u}) [GeV]", "p_{z}(#bar{u}) [GeV]", "E(#bar{u}) [GeV]", "m(#bar{u}) [GeV]", "p_{T}(#bar{u}) [GeV]", "y(#bar{u})", "#eta(#bar{u})", "#phi(#bar{u}) [rad]", "m^{reco}(#bar{u}) [GeV]", "p_{x}(c) [GeV]", "p_{y}(c) [GeV]", "p_{z}(c) [GeV]", "E(c) [GeV]", "m(c) [GeV]", "p_{T}(c) [GeV]", "y(c)", "#eta(c)", "#phi(c) [rad]", "m^{reco}(c) [GeV]", "p_{x}(#bar{c}) [GeV]", "p_{y}(#bar{c}) [GeV]", "p_{z}(#bar{c}) [GeV]", "E(#bar{c}) [GeV]", "m(#bar{c}) [GeV]", "p_{T}(#bar{c}) [GeV]", "y(#bar{c})", "#eta(#bar{c})", "#phi(#bar{c}) [rad]", "m^{reco}(#bar{c}) [GeV]",
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
        latex.DrawLatexNDC(0.79, 0.91, "#font[42]{FH channel}")


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
    root_files = ["lhe_production_FH_r1.0_f1.0_m172.5_p320900_plots_linear.root", "lhe_production_FH_r1.0_f2.0_m172.5_p320900_plots_linear.root", "lhe_hdampUP_FH_r1.0_f1.0_m172.5_p320900_plots_linear.root", "lhe_hdampUP_FH_r1.0_f2.0_m172.5_p320900_plots_linear.root", "lhe_hdampDOWN_FH_r1.0_f1.0_m172.5_p320900_plots_linear.root",  "lhe_hdampDOWN_FH_r1.0_f2.0_m172.5_p320900_plots_linear.root"]
    output_overlay_directory = "ratio_all_linear_FH"

    if not os.path.exists(output_overlay_directory):
        os.makedirs(output_overlay_directory)

    overlay_and_ratio_plot(root_files, output_overlay_directory)
