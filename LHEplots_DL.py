import ROOT
import argparse
import pathlib
import os
import math

def calculate_pt(px, py):
    return math.sqrt(px**2 + py**2)

def calculate_rapidity(energy, pz):
    if energy - pz <= 0 or (energy + pz) <= 0:
        return float('nan')
    
    return 0.5 * math.log((energy + pz) / (energy - pz))

def calculate_pseudorapidity(px,py,pz):
    p_total = math.sqrt(px**2 + py**2 + pz**2)
    if p_total - pz<=0 or (p_total + pz) <=0:
        return float('nan')

    return 0.5 * math.log((p_total + pz)/(p_total - pz))

def calculate_phi(px, py):
    return math.atan2(py, px)

def create_histograms(particle_name, particles_info):
    histograms = {} 
    for var, info in particles_info.items():
        hist = ROOT.TH1F(f"{particle_name}_{var}", info["title"], info["bins"], info["xmin"], info["xmax"])
        hist.GetXaxis().SetTitle(info.get("xlabel", ""))
        histograms[var] = hist
    return histograms

def add_overflow(hist):
    nbins = hist.GetNbinsX()+1
    e1 = hist.GetBinError(nbins-1)
    e2 = hist.GetBinError(nbins)
    hist.AddBinContent(nbins-1, hist.GetBinContent(nbins))
    hist.SetBinError(nbins-1, math.sqrt(e1*e1 + e2*e2))
    hist.SetBinContent(nbins, 0)
    hist.SetBinError(nbins, 0)
    return hist

def fill_histograms(histograms, values):
    for hist, value in zip(histograms, values):
        hist.Fill(value)

def fill_particle_histograms(histograms, event, j):
    for var, hist in histograms.items():
        if var == "px":
            fill_histograms([hist], [event.px[j]])
        elif var == "py":
            fill_histograms([hist], [event.py[j]])
        elif var == "pz":
            fill_histograms([hist], [event.pz[j]])
        elif var == "energy":
            fill_histograms([hist], [event.energy[j]])
        elif var == "mass":
            fill_histograms([hist], [event.mass[j]])
        elif var == "pt":
            pt = calculate_pt(event.px[j], event.py[j])
            fill_histograms([hist], [pt])
        elif var == "rapidity":
            rapidity = calculate_rapidity(event.energy[j],event.pz[j])
            fill_histograms([hist], [rapidity])
        elif var =="pseudorapidity":
            pseudorapidity = calculate_pseudorapidity(event.px[j], event.py[j], event.pz[j])
            fill_histograms([hist], [pseudorapidity])
        elif var == "phi":
            phi = calculate_phi(event.px[j], event.py[j])
            fill_histograms([hist], [phi])
        elif var =="pdgid":
            fill_histograms([hist], [event.pid[j]])

def overflow_Write(histograms):
    for hist in histograms.values():
        add_overflow(hist)
        hist.Write()

def create_canvas(hist, canvas_title, output_file,log_scale=False):
    canvas = ROOT.TCanvas(canvas_title, canvas_title, 800, 600)

    hist.SetStats(0)
    stats = ROOT.TPaveStats(0.6, 0.75, 0.9, 0.9, "NDC")
    stats.SetFillColor(0)
    stats.SetTextColor(ROOT.kBlack)  # Set text color
    stats.SetTextFont(42)            # Set font (42 is the default font in ROOT)
    stats.SetTextSize(0.033)
    stats.SetBorderSize(1)
    stats.SetLineWidth(1)
    hist.GetListOfFunctions().Add(stats)
    ROOT.gStyle.SetOptStat("RMe")
    hist.SetStats(1)

    # For 2D histograms, use "colz" option
    if "TH2" in hist.ClassName():
        hist.Draw("colz")
        hist.GetXaxis().SetTitle("[GeV]")
        hist.GetYaxis().SetTitle("[GeV]")
    else:
        hist.Draw("HIST")
        hist.GetYaxis().SetTitle("Entries")

    if log_scale:
        canvas.SetLogy()

    latex = ROOT.TLatex()
    latex.SetTextSize(0.04)
    latex.DrawLatexNDC(0.6, 0.91, "#bf{ #font[22]{CMS} #font[72]{Work In Progress} }")

    canvas.SaveAs(output_file)

def draw(histograms, name, output_directory):
    for var, hist in histograms.items():
        create_canvas(hist, f"canvas_{name}_{var}", f"{output_directory}/{name}_{var}.png", False)

def draw_log(histograms, name, output_directory):
    for var, hist in histograms.items():
        create_canvas(hist, f"canvas_{name}_{var}", f"{output_directory}/{name}_{var}.png", True)

def create_plots(input_root_file, output_directory):
    ROOT.gROOT.SetBatch(True) #not show canvas when draw 

    root_file = ROOT.TFile.Open(input_root_file)

    if not root_file or root_file.IsZombie():
        print(f"Error: Unable to open or read file: {input_root_file}")
        return

    tree = root_file.Get("events")

    if not tree:
        print("Error: TTree 'events' not found in the ROOT file.")
        return

    #====================================================================================================================
    particles_info = {
        "px": {"title": "Particle Momentum (px)", "bins": 100, "xmin": -500, "xmax": 500, "xlabel": "px [GeV]"},
        "py": {"title": "Particle Momentum (py)", "bins": 100, "xmin": -500, "xmax": 500, "xlabel": "py [GeV]"},
        "pz": {"title": "Particle Momentum (pz)", "bins": 100, "xmin": -1500, "xmax": 1500, "xlabel": "pz [GeV]"},
        "energy": {"title": "Particle Energy", "bins": 100, "xmin": 0, "xmax": 2500, "xlabel": "Energy [GeV]"},
        "mass": {"title": "Particle Mass", "bins": 100, "xmin": 0, "xmax": 200, "xlabel": "Mass [GeV]"},
        "pt": {"title": "Particle Transverse Momentum (pT)", "bins": 100, "xmin": 0, "xmax": 800, "xlabel": "pt [GeV]"},
        "rapidity": {"title": "Particle Rapidity", "bins": 100, "xmin": -6, "xmax": 6, "xlabel": "Rapidity"},
        "pseudorapidity": {"title": "Particle Pseudorapidity", "bins": 100, "xmin": -6, "xmax": 6, "xlabel": "Pseudorapidity"},
        "phi": {"title": "Particle Phi Angle", "bins": 100, "xmin": -math.pi, "xmax": math.pi, "xlabel": "Phi [rad]"},
        "pdgid": {"title": "Particle ID", "bins": 60, "xmin": -30, "xmax": 30, "xlabel": "PDG ID"},
    }

    top_quark_info = {
        "px": {"title": "Top Quark Momentum (px)", "bins": 100, "xmin": -500, "xmax": 500, "xlabel": "px [GeV]"},
        "py": {"title": "Top Quark Momentum (py)", "bins": 100, "xmin": -500, "xmax": 500, "xlabel": "py [GeV]"},
        "pz": {"title": "Top Quark Momentum (pz)", "bins": 100, "xmin": -2500, "xmax": 2500, "xlabel": "pz [GeV]"},
        "energy": {"title": "Top Quark Energy", "bins": 100, "xmin": 0, "xmax": 2500, "xlabel": "Energy [GeV]"},
        "mass": {"title": "Particle Mass", "bins": 100, "xmin": 150, "xmax": 200, "xlabel": "Mass [GeV]"},
        "pt": {"title": "Top Quark Transverse Momentum (pT)", "bins": 100, "xmin": 0, "xmax": 800, "xlabel": "pt [GeV]"},
        "rapidity": {"title": "Top Quark Rapidity", "bins": 100, "xmin": -5, "xmax": 5, "xlabel": "Rapidity"},
        "phi": {"title": "Top Quark Phi Angle", "bins": 100, "xmin": -math.pi, "xmax": math.pi, "xlabel": "Phi [rad]"},
    }

    W_boson_info={
        "px": {"title": "W boson Momentum (px)", "bins": 100, "xmin": -500, "xmax": 500, "xlabel": "px [GeV]"},
        "py": {"title": "W boson Momentum (py)", "bins": 100, "xmin": -500, "xmax": 500, "xlabel": "py [GeV]"},
        "pz": {"title": "W boson Momentum (pz)", "bins": 100, "xmin": -2000, "xmax": 2000, "xlabel": "pz [GeV]"},
        "energy": {"title": "W boson Energy", "bins": 100, "xmin": 0, "xmax": 2500, "xlabel": "Energy [GeV]"},
        "mass": {"title": "Particle Mass", "bins": 100, "xmin": 60, "xmax": 100, "xlabel": "Mass [GeV]"},
        "pt": {"title": "W boson Transverse Momentum (pT)", "bins": 100, "xmin": 0, "xmax": 800, "xlabel": "pt [GeV]"},
        "rapidity": {"title": "W boson Rapidity", "bins": 100, "xmin": -5, "xmax": 5, "xlabel": "Rapidity"},
        "phi": {"title": "W boson Phi Angle", "bins": 100, "xmin": -math.pi, "xmax": math.pi, "xlabel": "Phi [rad]"},
    }

    bottom_quark_info = {
        "px": {"title": "Bottom Quark Momentum (px)", "bins": 100, "xmin": -600, "xmax": 600, "xlabel": "px [GeV]"},
        "py": {"title": "Bottom Quark Momentum (py)", "bins": 100, "xmin": -600, "xmax": 600, "xlabel": "py [GeV]"},
        "pz": {"title": "Bottom Quark Momentum (pz)", "bins": 100, "xmin": -2500, "xmax": 2500, "xlabel": "pz [GeV]"},
        "energy": {"title": "Bottom Quark Energy", "bins": 100, "xmin": 0, "xmax": 2500, "xlabel": "Energy [GeV]"},
        "mass": {"title": "Bottom Quark Mass", "bins": 100, "xmin": 0, "xmax": 8, "xlabel": "Mass [GeV]"},
        "pt": {"title": "Bottom Quark Transverse Momentum (pT)", "bins": 100, "xmin": 0, "xmax": 800, "xlabel": "pt [GeV]"},
        "rapidity": {"title": "Bottom Quark Rapidity", "bins": 100, "xmin": -6, "xmax": 6, "xlabel": "Rapidity"},
        "pseudorapidity": {"title": "Bottom Quark Pseudorapidity", "bins": 100, "xmin": -6, "xmax": 6, "xlabel": "Pseudorapidity"},
        "phi": {"title": "Bottom Quark Phi Angle", "bins": 100, "xmin": -math.pi, "xmax": math.pi, "xlabel": "Phi [rad]"},
    }

    leptons_info = {
        "px": {"title": "Momentum (px)", "bins": 100, "xmin": -500, "xmax": 500, "xlabel": "px [GeV]"},
        "py": {"title": "Momentum (py)", "bins": 100, "xmin": -500, "xmax": 500, "xlabel": "py [GeV]"},
        "pz": {"title": "Momentum (pz)", "bins": 100, "xmin": -2000, "xmax": 2000, "xlabel": "pz [GeV]"},
        "energy": {"title": "Energy", "bins": 100, "xmin": 0, "xmax": 1500, "xlabel": "Energy [GeV]"},
        "mass": {"title": "Mass", "bins": 100, "xmin": 0, "xmax": 2, "xlabel": "Mass [GeV]"},
        "pt": {"title": "Transverse Momentum (pT)", "bins": 100, "xmin": 0, "xmax": 800, "xlabel": "pt [GeV]"},
        "rapidity": {"title": "Rapidity", "bins": 100, "xmin": -5, "xmax": 5, "xlabel": "Rapidity"},
        "pseudorapidity": {"title": "Pseudorapidity", "bins": 100, "xmin": -5, "xmax": 5, "xlabel": "Pseudorapidity"},
        "phi": {"title": "Phi Angle", "bins": 100, "xmin": -math.pi, "xmax": math.pi, "xlabel": "Phi [rad]"},
    }

    #====================================================================================================================

    all_particles_histograms = create_histograms("all", particles_info) 

    initial_particles_histograms = create_histograms("initial_particles", particles_info)
    radiation_histograms = create_histograms("radiation", particles_info)

    top_quark_histograms = create_histograms("top", top_quark_info)
    antitop_quark_histograms = create_histograms("antitop", top_quark_info)
    ttbar_histograms = create_histograms("ttbar", top_quark_info)

    W_boson_histograms = create_histograms("Wp", W_boson_info)
    Wm_boson_histograms = create_histograms("Wm", W_boson_info)
    Wpm_histograms = create_histograms("W", W_boson_info)

    b_all_histograms = create_histograms("b_all", bottom_quark_info)
    bottom_quark_histograms = create_histograms("bottom", bottom_quark_info)
    antibottom_quark_histograms = create_histograms("antibottom", bottom_quark_info)
    b_from_top_histograms = create_histograms("b_from_top", bottom_quark_info)
    b_from_initial_histograms = create_histograms("b_from_initial", bottom_quark_info)
    # maybe add bottom antibottom separately for the last two

    leptons_histograms = create_histograms("leptons", leptons_info)
    charged_leptons_histograms = create_histograms("charged_leptons", leptons_info)
    neutrinos_histograms = create_histograms("neutrinos", leptons_info)
    electron_pairs_histograms = create_histograms("electron_pairs", leptons_info)
    electrons_histograms = create_histograms("electrons", leptons_info)
    positrons_histograms = create_histograms("positrons", leptons_info)
    muon_pairs_histograms = create_histograms("muon_pairs", leptons_info)
    muons_histograms = create_histograms("muons", leptons_info)
    antimuons_histograms = create_histograms("antimuons", leptons_info)
    tau_pairs_histograms = create_histograms("tau_pairs", leptons_info)
    taus_histograms = create_histograms("taus", leptons_info)
    antitaus_histograms = create_histograms("antitaus", leptons_info)
    electron_neutrinos_histograms = create_histograms("electron_neutrinos", leptons_info)
    muon_neutrinos_histograms = create_histograms("muon_neutrinos", leptons_info)
    tau_neutrinos_histograms = create_histograms("tau_neutrinos", leptons_info)
    #maybe add neutrino and antineutrino saparately

    #============ Include all the above histograms in the list 
    histogram_dict = {
        "all_particles": all_particles_histograms,
        "initial_particles": initial_particles_histograms,
        "radiation": radiation_histograms,
        "top_quark": top_quark_histograms,
        "antitop_quark": antitop_quark_histograms,
        "ttbar": ttbar_histograms,
        "W_boson": W_boson_histograms,
        "Wm_boson": Wm_boson_histograms,
        "Wpm": Wpm_histograms,
        "b_all": b_all_histograms,
        "bottom_quark": bottom_quark_histograms,
        "antibottom_quark": antibottom_quark_histograms,
        "b_from_top": b_from_top_histograms,
        "b_from_initial": b_from_initial_histograms,
        "leptons": leptons_histograms,
        "charged_leptons": charged_leptons_histograms,
        "neutrinos": neutrinos_histograms,
        "electron_pairs": electron_pairs_histograms,
        "electrons": electrons_histograms,
        "positrons": positrons_histograms,
        "muon_pairs": muon_pairs_histograms,
        "muons": muons_histograms,
        "antimuons": antimuons_histograms,
        "tau_pairs": tau_pairs_histograms,
        "taus": taus_histograms,
        "antitaus": antitaus_histograms,
        "electron_neutrinos": electron_neutrinos_histograms,
        "muon_neutrinos": muon_neutrinos_histograms,
        "tau_neutrinos": tau_neutrinos_histograms
    }

    for event in tree:
        for i in range(len(event.px)):
            fill_particle_histograms(all_particles_histograms, event, i)

            if event.mother1[i] == 0:
                fill_particle_histograms(initial_particles_histograms, event, i)
            if event.pid[i] == 6 or event.pid[i] == -6:
                fill_particle_histograms(ttbar_histograms, event, i)
            if event.pid[i] == 6:
                fill_particle_histograms(top_quark_histograms, event, i)
            if event.pid[i] == -6:
                fill_particle_histograms(antitop_quark_histograms, event, i)
            
            if event.pid[i] == 24 or event.pid[i] == -24:
                fill_particle_histograms(Wpm_histograms, event, i)
            if event.pid[i] == 24:
                fill_particle_histograms(W_boson_histograms, event, i)
            if event.pid[i] == -24:
                fill_particle_histograms(Wm_boson_histograms, event, i)

            if event.pid[i] == 5 or event.pid[i] == -5:
                fill_particle_histograms(b_all_histograms, event, i)
                if event.mother1[i] > 2: 
                    fill_particle_histograms(b_from_top_histograms, event, i)
                else:
                    fill_particle_histograms(b_from_initial_histograms, event, i)
            if event.pid[i] == 5:
                fill_particle_histograms(bottom_quark_histograms, event, i)
            if event.pid[i] == -5:
                fill_particle_histograms(antibottom_quark_histograms, event, i)
            
            if event.pid[i] == 11 or event.pid[i] == -11 or event.pid[i] == 13 or event.pid[i] == -13 or event.pid[i] == 12 or event.pid[i] == -12 or event.pid[i] == 14 or event.pid[i] == -14 or event.pid[i] == 15 or event.pid[i] == -15 or event.pid[i] == 16 or event.pid[i] == -16:
                fill_particle_histograms(leptons_histograms, event, i)
            if event.pid[i] == 11 or event.pid[i] == -11 or event.pid[i] == 13 or event.pid[i] == -13 or event.pid[i] == 15 or event.pid[i] == -15:
                fill_particle_histograms(charged_leptons_histograms, event, i)
            if event.pid[i] == 11 or event.pid[i] == -11:
                fill_particle_histograms(electron_pairs_histograms, event, i)
            if event.pid[i] == 11:
                fill_particle_histograms(electrons_histograms, event, i)
            if event.pid[i] == -11:
                fill_particle_histograms(positrons_histograms,event, i)
            if event.pid[i] == 13 or event.pid[i] == -13:
                fill_particle_histograms(muon_pairs_histograms, event, i)
            if event.pid[i] == 13:
                fill_particle_histograms(muons_histograms, event, i)
            if event.pid[i] == -13:
                fill_particle_histograms(antimuons_histograms,event, i)
            if event.pid[i] == 15 or event.pid[i] == -15:
                fill_particle_histograms(tau_pairs_histograms, event, i)
            if event.pid[i] == 15:
                fill_particle_histograms(taus_histograms, event, i)
            if event.pid[i] == -15:
                fill_particle_histograms(antitaus_histograms,event, i)

            if event.pid[i] == 12 or event.pid[i] == -12 or event.pid[i] == 14 or event.pid[i] == -14 or event.pid[i] == 16 or event.pid[i] == -16:
                fill_particle_histograms(neutrinos_histograms, event, i)
            if event.pid[i] == 12 or event.pid[i] == -12:
                fill_particle_histograms(electron_neutrinos_histograms, event, i)
            if event.pid[i] == 14 or event.pid[i] == -14:
                fill_particle_histograms(muon_neutrinos_histograms, event, i)
            if event.pid[i] == 16 or event.pid[i] == -16:
                fill_particle_histograms(tau_neutrinos_histograms, event, i)

            if event.status[i] == 1:
                if (event.pid[i] <5 and event.pid[i] >-5) or event.pid[i] == 21:
                    fill_particle_histograms(radiation_histograms, event, i)

    output_root_file = f"{output_directory}.root"
    output_file = ROOT.TFile.Open(output_root_file, "RECREATE")
    output_file.cd()

    for name, histograms in histogram_dict.items():
        overflow_Write(histograms)

    output_file.Close()

    for name, histograms in histogram_dict.items():
        draw(histograms, name, output_directory)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, help='Input file path and Name')
    parser.add_argument("--output", type=str, help='Output directory path')
    args = parser.parse_args()

    if pathlib.PurePosixPath(args.input).suffix.lower() != ".root":
        print('Error: Only .root files are accepted as input. Exiting ...')
    else:
        input_root_file = args.input
        output_directory = args.output or "output"

        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        create_plots(input_root_file, output_directory)
