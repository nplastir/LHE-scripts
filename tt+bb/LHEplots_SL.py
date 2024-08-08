import ROOT
import argparse
import pathlib
import os
import math

#---------------------------------------------------------------------------------------------------------------------
## Functions necessary for the calculation of physical quantities


# Function to calculate the HT of a particle system
def calculate_HT(pTs):
    return sum(pTs)

# Function to calculate the invariant mass of a particle
def calculate_mass(E, px, py, pz):
    p_total = math.sqrt(px**2 + py**2 + pz**2)
    mass_sq =abs(E**2 - p_total**2)
    #mass = math.sqrt(max(0, mass_sq))  # Ensure non-negative square root
    mass = math.sqrt(mass_sq)
    return mass

# Function to calculate the invariant mass of a two particle system
def calculate_invariant_mass(px1, py1, pz1, energy1, px2, py2, pz2, energy2):
    # m^2 = (E1 + E2)^2 - (px1 + px2)^2 - (py1 + py2)^2 - (pz1 + pz2)^2
    invariant_mass_sq = (energy1 + energy2)**2 - (px1 + px2)**2 - (py1 + py2)**2 - (pz1 + pz2)**2
    invariant_mass = math.sqrt(max(0, invariant_mass_sq))  # Ensure non-negative square root
    return invariant_mass

# Function to calculate the dphi between two particles
def calculate_dphi(px1,py1,px2,py2):
    dphi = abs(calculate_phi(px1, py1) - calculate_phi(px2, py2))
    if dphi > math.pi: dphi = 2*math.pi-dphi
    ##Alternative way of calculating dphi
    # if abs((px1*px2 +py1*py2)/(calculate_pt(px1, py1)*calculate_pt(px2,py2))) >1:
    #     return float('nan')
    # dphi = math.acos((px1*px2 +py1*py2)/(calculate_pt(px1, py1)*calculate_pt(px2,py2)))
    return dphi

# Function to calculate the deta between two particles
def calculate_deta(px1,py1,pz1,px2,py2,pz2):
    deta = abs(calculate_pseudorapidity(px1, py1, pz1) - calculate_pseudorapidity(px2, py2, pz2))
    return deta

# Function to calculate the dR between two particles
def calculate_dR(px1, py1, pz1, energy1, px2, py2, pz2, energy2):
    # Definition: dR = sqrt(deta^2 + dphi^2)
    deta = calculate_deta(px1,py1,pz1,px2,py2,pz2)
    dphi = calculate_dphi(px1,py1,px2,py2)
    dR = math.sqrt(deta**2 + dphi**2)
    return dR

# Function to calculate the dR between two particles (Here we use dy instead of the deta. This is only for tests. Do not use this for the dR calculation)
def calculate_dR_2(px1, py1, pz1, energy1, px2, py2, pz2, energy2):
    dy = abs(calculate_rapidity(energy1, pz1) - calculate_rapidity(energy2, pz2))
    dphi = calculate_dphi(px1,py1,px2,py2)
    dR = math.sqrt(dy**2 + dphi**2)
    return dR

# Function to calculate the transverse momentum of a particle
def calculate_pt(px, py):
    return math.sqrt(px**2 + py**2)

# Function to calculate the rapidity of a particle
def calculate_rapidity(energy, pz):
    if energy - pz <= 0 or (energy + pz) <= 0:
        return float('nan')
    
    return 0.5 * math.log((energy + pz) / (energy - pz))

# Function to calculate the pseudorapidity of a particle (Use this only for massless particles, where Ε~p )
def calculate_pseudorapidity(px,py,pz):
    p_total = math.sqrt(px**2 + py**2 + pz**2)
    if p_total - pz<=0 or (p_total + pz) <=0:
        return float('nan')

    return 0.5 * math.log((p_total + pz)/(p_total - pz))

# Function to calculate the azimuthal angle of a particle
def calculate_phi(px, py):
    return math.atan2(py, px)

#---------------------------------------------------------------------------------------------------------------------
## Functions that automate the process of the creation of canvases and histograms

# Function that automatically defines all the 1D-histograms of a dictionary
def create_histograms(particle_name, particles_info):
    histograms = {} 
    for var, info in particles_info.items():
        hist = ROOT.TH1F(f"{particle_name}_{var}", info["title"], info["bins"], info["xmin"], info["xmax"])
        hist.GetXaxis().SetTitle(info.get("xlabel", ""))
        histograms[var] = hist
    return histograms

# Function that automatically defines all the 2D-histograms of a dictionary
def create_histograms2(particle_name, particles_info):
    histograms = {} 
    for var, info in particles_info.items():
        hist = ROOT.TH2F(f"{particle_name}_{var}", info["title"], info["bins1"], info["xmin1"], info["xmax1"], info["bins2"], info["xmin2"], info["xmax2"])
        hist.GetXaxis().SetTitle(info.get("xlabel", ""))
        hist.GetYaxis().SetTitle(info.get("ylabel", ""))
        histograms[var] = hist
    return histograms

# Function that adds the overflow in the last bin
def add_overflow(hist):
    nbins = hist.GetNbinsX()+1
    e1 = hist.GetBinError(nbins-1)
    e2 = hist.GetBinError(nbins)
    hist.AddBinContent(nbins-1, hist.GetBinContent(nbins))
    hist.SetBinError(nbins-1, math.sqrt(e1*e1 + e2*e2))
    hist.SetBinContent(nbins, 0)
    hist.SetBinError(nbins, 0)
    return hist

# Function that adds the underflow in the first bin
def add_underflow(hist):
    e1 = hist.GetBinError(1)
    e0 = hist.GetBinError(0)
    hist.AddBinContent(1, hist.GetBinContent(0))
    hist.SetBinError(1, math.sqrt(e1 * e1 + e0 * e0))
    hist.SetBinContent(0, 0)
    hist.SetBinError(0, 0)
    return hist

# Function that fills the histograms 
def fill_histograms(histograms, values):
    for hist, value in zip(histograms, values):
        hist.Fill(value)

# Function that fills all the histograms of particles
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
        elif var =="reco_mass":
            m = calculate_mass(event.energy[j],event.px[j], event.py[j], event.pz[j])
            fill_histograms([hist], [m])

# Function that fills all the histograms of two particle systems
def fill_system_histograms(histograms, event, k, l): # k: particle index, l: antiparticle index 
    for var, hist in histograms.items():
        if var == "pt":
            pt = calculate_pt(event.px[k] + event.px[l], event.py[k] + event.py[l])
            fill_histograms([hist], [pt])
        elif var == "rapidity":
            rapidity = calculate_rapidity(event.energy[k] + event.energy[l], event.pz[k] + event.pz[l])
            fill_histograms([hist], [rapidity])
        elif var == "pseudorapidity":
            pseudorapidity = calculate_pseudorapidity(event.px[k] + event.px[l], event.py[k] + event.py[l], event.pz[k] + event.pz[l])
            fill_histograms([hist], [pseudorapidity])
        elif var == "phi":
            phi = calculate_phi(event.px[k] + event.px[l], event.py[k] + event.py[l])
            fill_histograms([hist], [phi])
        elif var == "dR":
            dR = calculate_dR(event.px[k], event.py[k], event.pz[k], event.energy[k],event.px[l], event.py[l], event.pz[l], event.energy[l])
            fill_histograms([hist], [dR])
        elif var == "dR_2":
            dR = calculate_dR_2(event.px[k], event.py[k], event.pz[k], event.energy[k],event.px[l], event.py[l], event.pz[l], event.energy[l])
            fill_histograms([hist], [dR])
        elif var == "invariant_mass":
            invariant_mass = calculate_invariant_mass(event.px[k], event.py[k], event.pz[k], event.energy[k], event.px[l], event.py[l], event.pz[l], event.energy[l])
            fill_histograms([hist], [invariant_mass])
        elif var == "HT":
            HT = calculate_HT([calculate_pt(event.px[k], event.py[k]),calculate_pt(event.px[l], event.py[l])])
            fill_histograms([hist], [HT])

# Function that adds the overflow bin and writes the histogram
def overflow_Write(histograms):
    for hist in histograms.values():
        add_overflow(hist)
        hist.Write()

# Function that automatically creates the canvas needed 
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

    ROOT.gStyle.SetOptTitle(0)
    ROOT.TGaxis.SetExponentOffset(-0.05, 0.01, "Y")
    latex = ROOT.TLatex()
    latex.SetTextSize(0.043)
    latex.DrawLatexNDC(0.1, 0.91, "#font[61]{CMS}")
    latex.SetTextSize(0.033)
    latex.DrawLatexNDC(0.17, 0.91, "#font[52]{Simulation Work-in-Progress}")

    canvas.SaveAs(f"{output_file}.png")
    canvas.SaveAs(f"{output_file}.pdf")

# Function that automatically draws and saves all the histograms
def draw(histograms, name, output_directory):
    for var, hist in histograms.items():
        create_canvas(hist, f"canvas_{name}_{var}", f"{output_directory}/{name}_{var}", False)

# Function that automatically draws and saves all the histograms in log scale 
def draw_log(histograms, name, output_directory):
    for var, hist in histograms.items():
        create_canvas(hist, f"canvas_{name}_{var}", f"{output_directory}/{name}_{var}", True)

#---------------------------------------------------------------------------------------------------------------------
## Main function 

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
        "px": {"title": "Particle Momentum (px)", "bins": 50, "xmin": -200, "xmax": 200, "xlabel": "p_{x} [GeV]"},
        "py": {"title": "Particle Momentum (py)", "bins": 50, "xmin": -200, "xmax": 200, "xlabel": "p_{y} [GeV]"},
        "pz": {"title": "Particle Momentum (pz)", "bins": 50, "xmin": -1100, "xmax": 1100, "xlabel": "p_{z} [GeV]"},
        "energy": {"title": "Particle Energy", "bins": 50, "xmin": 0, "xmax": 1500, "xlabel": "E [GeV]"},
        "mass": {"title": "Particle Mass", "bins": 50, "xmin": 0, "xmax": 200, "xlabel": "m [GeV]"},
        "pt": {"title": "Particle Transverse Momentum (pT)", "bins": 50, "xmin": 0, "xmax": 400, "xlabel": "p_{T} [GeV]"},
        "rapidity": {"title": "Particle Rapidity", "bins": 50, "xmin": -7, "xmax": 7, "xlabel": "y"},
        "pseudorapidity": {"title": "Particle Pseudorapidity", "bins": 50, "xmin": -7, "xmax": 7, "xlabel": "#eta"},
        "phi": {"title": "Particle Phi Angle", "bins": 50, "xmin": -math.pi, "xmax": math.pi, "xlabel": "#phi [rad]"},
        "pdgid": {"title": "Particle ID", "bins": 30, "xmin": -30, "xmax": 30, "xlabel": "PDG ID"},
    }

    top_quark_info = {
        "px": {"title": "Top Quark Momentum (px)", "bins": 50, "xmin": -500, "xmax": 500, "xlabel": "p_{x}(t) [GeV]"},
        "py": {"title": "Top Quark Momentum (py)", "bins": 50, "xmin": -500, "xmax": 500, "xlabel": "p_{y}(t) [GeV]"},
        "pz": {"title": "Top Quark Momentum (pz)", "bins": 50, "xmin": -1500, "xmax": 1500, "xlabel": "p_{z}(t) [GeV]"},
        "energy": {"title": "Top Quark Energy", "bins": 50, "xmin": 0, "xmax": 1500, "xlabel": "E(t) [GeV]"},
        "mass": {"title": "Particle Mass", "bins": 50, "xmin": 150, "xmax": 200, "xlabel": "m(t) [GeV]"},
        "pt": {"title": "Top Quark Transverse Momentum (pT)", "bins": 50, "xmin": 0, "xmax": 600, "xlabel": "p_{T}(t) [GeV]"},
        "rapidity": {"title": "Top Quark Rapidity", "bins": 50, "xmin": -5, "xmax": 5, "xlabel": "y(t)"},
        "phi": {"title": "Top Quark Phi Angle", "bins": 50, "xmin": -math.pi, "xmax": math.pi, "xlabel": "#phi(t) [rad]"},
        "reco_mass": {"title": "Mass", "bins": 50, "xmin": 150, "xmax": 200, "xlabel": "m^{reco}(t) [GeV]"},
    }

    W_boson_info={
        "px": {"title": "W boson Momentum (px)", "bins": 50, "xmin": -300, "xmax": 300, "xlabel": "p_{x}(W) [GeV]"},
        "py": {"title": "W boson Momentum (py)", "bins": 50, "xmin": -300, "xmax": 300, "xlabel": "p_{y}(W) [GeV]"},
        "pz": {"title": "W boson Momentum (pz)", "bins": 50, "xmin": -1000, "xmax": 1000, "xlabel": "p_{z}(W) [GeV]"},
        "energy": {"title": "W boson Energy", "bins": 50, "xmin": 0, "xmax": 1000, "xlabel": "E(W) [GeV]"},
        "mass": {"title": "Particle Mass", "bins": 50, "xmin": 60, "xmax": 100, "xlabel": "m(W) [GeV]"},
        "pt": {"title": "W boson Transverse Momentum (pT)", "bins": 50, "xmin": 0, "xmax": 400, "xlabel": "p_{T}(W) [GeV]"},
        "rapidity": {"title": "W boson Rapidity", "bins": 50, "xmin": -5, "xmax": 5, "xlabel": "y(W)"},
        "phi": {"title": "W boson Phi Angle", "bins": 50, "xmin": -math.pi, "xmax": math.pi, "xlabel": "#phi(W) [rad]"},
        "reco_mass": {"title": "Mass", "bins": 50, "xmin": 60, "xmax": 100, "xlabel": "m^{reco}(W) [GeV]"},
    }

    bottom_quark_info = {
        "px": {"title": "Bottom Quark Momentum (px)", "bins": 50, "xmin": -200, "xmax": 200, "xlabel": "p_{x}(b) [GeV]"},
        "py": {"title": "Bottom Quark Momentum (py)", "bins": 50, "xmin": -200, "xmax": 200, "xlabel": "p_{y}(b) [GeV]"},
        "pz": {"title": "Bottom Quark Momentum (pz)", "bins": 50, "xmin": -800, "xmax": 800, "xlabel": "p_{z}(b) [GeV]"},
        "energy": {"title": "Bottom Quark Energy", "bins": 50, "xmin": 0, "xmax": 1000, "xlabel": "E(b) [GeV]"},
        "mass": {"title": "Bottom Quark Mass", "bins": 50, "xmin": 0, "xmax": 8, "xlabel": "m(b) [GeV]"},
        "pt": {"title": "Bottom Quark Transverse Momentum (pT)", "bins": 50, "xmin": 0, "xmax": 300, "xlabel": "p_{T}(b) [GeV]"},
        "rapidity": {"title": "Bottom Quark Rapidity", "bins": 50, "xmin": -7, "xmax": 7, "xlabel": "y(b)"},
        "pseudorapidity": {"title": "Bottom Quark Pseudorapidity", "bins": 50, "xmin": -7, "xmax": 7, "xlabel": "#eta(b)"},
        "phi": {"title": "Bottom Quark Phi Angle", "bins": 50, "xmin": -math.pi, "xmax": math.pi, "xlabel": "#phi(b) [rad]"},
        "reco_mass": {"title": "Mass", "bins": 50, "xmin": 0, "xmax": 5, "xlabel": "m^{reco}(b) [GeV]"},
    }

    leptons_info = {
        "px": {"title": "Momentum (px)", "bins": 50, "xmin": -200, "xmax": 200, "xlabel": "p_{x} [GeV]"},
        "py": {"title": "Momentum (py)", "bins": 50, "xmin": -200, "xmax": 200, "xlabel": "p_{y} [GeV]"},
        "pz": {"title": "Momentum (pz)", "bins": 50, "xmin": -600, "xmax": 600, "xlabel": "p_{z} [GeV]"},
        "energy": {"title": "Energy", "bins": 50, "xmin": 0, "xmax": 650, "xlabel": "E [GeV]"},
        "mass": {"title": "Mass", "bins": 50, "xmin": 0, "xmax": 0.2, "xlabel": "m [GeV]"},
        "pt": {"title": "Transverse Momentum (pT)", "bins": 50, "xmin": 0, "xmax": 300, "xlabel": "p_{T} [GeV]"},
        "rapidity": {"title": "Rapidity", "bins": 50, "xmin": -6, "xmax": 6, "xlabel": "y"},
        "pseudorapidity": {"title": "Pseudorapidity", "bins": 50, "xmin": -6, "xmax": 6, "xlabel": "#eta"},
        "phi": {"title": "Phi Angle", "bins": 50, "xmin": -math.pi, "xmax": math.pi, "xlabel": "#phi [rad]"},
        "reco_mass": {"title": "Mass", "bins": 50, "xmin": 0, "xmax": 0.2, "xlabel": "m^{reco} [GeV]"},
    }

    quarks_info = {
        "px": {"title": "Momentum (px)", "bins": 50, "xmin": -200, "xmax": 200, "xlabel": "p_{x} [GeV]"},
        "py": {"title": "Momentum (py)", "bins": 50, "xmin": -200, "xmax": 200, "xlabel": "p_{y} [GeV]"},
        "pz": {"title": "Momentum (pz)", "bins": 50, "xmin": -600, "xmax": 600, "xlabel": "p_{z} [GeV]"},
        "energy": {"title": "Energy", "bins": 50, "xmin": 0, "xmax": 700, "xlabel": "Energy [GeV]"},
        "mass": {"title": "Mass", "bins": 50, "xmin": 0, "xmax": 2, "xlabel": "Mass [GeV]"},
        "pt": {"title": "Transverse Momentum (pT)", "bins": 50, "xmin": 0, "xmax": 300, "xlabel": "p_{T} [GeV]"},
        "rapidity": {"title": "Rapidity", "bins": 50, "xmin": -5, "xmax": 5, "xlabel": "Rapidity"},
        "pseudorapidity": {"title": "Pseudorapidity", "bins": 50, "xmin": -5, "xmax": 5, "xlabel": "Pseudorapidity"},
        "phi": {"title": "Phi Angle", "bins": 50, "xmin": -math.pi, "xmax": math.pi, "xlabel": "#phi [rad]"},
        "reco_mass": {"title": "Mass", "bins": 100, "xmin": 0, "xmax": 0.2, "xlabel": "Mass [GeV]"},
    }

    tt_system_info ={
        "pt": {"title": "Transverse Momentum of tt System", "bins": 50, "xmin": 0, "xmax": 600, "xlabel": "p_{T} [GeV]"},
        "rapidity": {"title": "Rapidity of tt System", "bins": 50, "xmin": -5, "xmax": 5, "xlabel": "y"},
        "HT": {"title": "HT", "bins": 50, "xmin": 0, "xmax": 1000, "xlabel": "H_{T} [GeV]"},
        "invariant_mass": {"title": "Invariant Mass", "bins": 50, "xmin": 0, "xmax": 1500, "xlabel": "m^{reco} [GeV]"},
        "dR": {"title": "Angular Separation (dR)", "bins": 50, "xmin": 0, "xmax": 6, "xlabel": "#DeltaR"},
        "dR_2": {"title": "Angular Separation (dR)", "bins": 50, "xmin": 0, "xmax": 6, "xlabel": "#DeltaR"},
    }
    
    bb_system_info={
        "pt": {"title": "Transverse Momentum of bb System", "bins": 50, "xmin": 0, "xmax": 1000, "xlabel": "p_{T} [GeV]"},
        "rapidity": {"title": "Rapidity of bb System", "bins": 50, "xmin": -5, "xmax": 5, "xlabel": "y"},
        "pseudorapidity": {"title": "Pseudorapidity of bb System", "bins": 50, "xmin": -6, "xmax": 6, "xlabel": "#eta"},
        "HT": {"title": "HT", "bins": 50, "xmin": 0, "xmax": 400, "xlabel": "H_{T} [GeV]"},
        "invariant_mass": {"title": "Invariant Mass", "bins": 50, "xmin": 0, "xmax": 600, "xlabel": "m^{reco} [GeV]"},
        "dR": {"title": "Angular Separation (dR)", "bins": 50, "xmin": 0, "xmax": 6, "xlabel": "#DeltaR"},
    }

    quarks_system_info={
        "pt": {"title": "Transverse Momentum ", "bins": 50, "xmin": 0, "xmax": 500, "xlabel": "p_{T} [GeV]"},
        "rapidity": {"title": "Rapidity ", "bins": 50, "xmin": -5, "xmax": 5, "xlabel": "Rapidity"},
        "pseudorapidity": {"title": "Pseudorapidity of ", "bins": 50, "xmin": -6, "xmax": 6, "xlabel": "Pseudorapidity"},
        "HT": {"title": "HT", "bins": 50, "xmin": 0, "xmax": 600, "xlabel": "H_{T} [GeV]"},
        "invariant_mass": {"title": "Invariant Mass", "bins": 50, "xmin": 0, "xmax": 600, "xlabel": "Mass [GeV]"},
        "dR": {"title": "Angular Separation (dR)", "bins": 50, "xmin": 0, "xmax": 6, "xlabel": "#DeltaR"},  
    }

    quark_pairs_system_info={
        "pt": {"title": "Transverse Momentum ", "bins": 50, "xmin": 0, "xmax": 500, "xlabel": "p_{T} [GeV]"},
        "rapidity": {"title": "Rapidity ", "bins": 50, "xmin": -5, "xmax": 5, "xlabel": "Rapidity"},
        "pseudorapidity": {"title": "Pseudorapidity of ", "bins": 50, "xmin": -6, "xmax": 6, "xlabel": "Pseudorapidity"},
        "HT": {"title": "HT", "bins": 50, "xmin": 0, "xmax": 600, "xlabel": "H_{T} [GeV]"},
        "invariant_mass": {"title": "Invariant Mass", "bins": 50, "xmin": 40, "xmax": 120, "xlabel": "Mass [GeV]"},
        "dR": {"title": "Angular Separation (dR)", "bins": 50, "xmin": 0, "xmax": 6, "xlabel": "#DeltaR"},  
    }
    ll_system_info={
        "pt": {"title": "Transverse Momentum of ll System", "bins": 50, "xmin": 0, "xmax": 300, "xlabel": "p_{T} [GeV]"},
        "rapidity": {"title": "Rapidity of ll System", "bins": 50, "xmin": -5, "xmax": 5, "xlabel": "y"},
        "pseudorapidity": {"title": "Pseudorapidity of ll System", "bins": 50, "xmin": -6, "xmax": 6, "xlabel": "#eta"},
        "HT": {"title": "HT", "bins": 50, "xmin": 0, "xmax": 400, "xlabel": "H_{T} [GeV]"},
        "invariant_mass": {"title": "Invariant Mass", "bins": 50, "xmin": 0, "xmax": 600, "xlabel": "m^{reco} [GeV]"},
        "dR": {"title": "Angular Separation (dR)", "bins": 50, "xmin": 0, "xmax": 6, "xlabel": "#DeltaR"},       
    }

    tt_system_2d_info={
        "phivseta": {"title": "DphivsDeta", "bins1": 100, "xmin1": 0, "xmax1": 6, "bins2": 30, "xmin2": 0, "xmax2": 3.14,  "xlabel": "|#Delta#eta|", "ylabel": "|#Delta#phi| [rad]" }
    }
    #====================================================================================================================

    # Particle systems 
    tt_system_histograms = create_histograms("tt_system", tt_system_info)
    tt_system_2d_histograms = create_histograms2("tt system 2d", tt_system_2d_info)
    bb_system_histograms = create_histograms("bb_system", bb_system_info)
    b_all_system_histograms = create_histograms("b_all_system", bb_system_info)
    prompt_bs_histograms = create_histograms("prompt_bs", bb_system_info)
    bs_from_top_histograms = create_histograms("bs_from_top", bb_system_info)

    quark_system_histograms = create_histograms("quarks_system", quark_pairs_system_info)
    ud_pairs_system_histograms = create_histograms("ud_pairs_system", quark_pairs_system_info)
    sc_pairs_system_histograms = create_histograms("sc_pairs_system", quark_pairs_system_info)
    cd_pairs_system_histograms = create_histograms("cd_pairs_system", quark_pairs_system_info)
    su_pairs_system_histograms = create_histograms("su_pairs_system", quark_pairs_system_info)

    # lnu_system_histograms = create_histograms("lnu_system", ll_system_info)
    e_system_histograms = create_histograms("e_system", ll_system_info)
    m_system_histograms = create_histograms("m_system", ll_system_info)
    t_system_histograms = create_histograms("t_system", ll_system_info)
    em_system_histograms = create_histograms("em_system", ll_system_info)
    mm_system_histograms = create_histograms("mm_system", ll_system_info)
    tm_system_histograms = create_histograms("tm_system", ll_system_info)
    ep_system_histograms = create_histograms("ep_system", ll_system_info)
    mp_system_histograms = create_histograms("mp_system", ll_system_info)
    tp_system_histograms = create_histograms("tp_system", ll_system_info)


    # Particles 
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

    quarks_histograms = create_histograms("quarks", quarks_info)
    d_quarks_histograms = create_histograms("d_quarks", quarks_info)
    u_quarks_histograms = create_histograms("u_quarks", quarks_info)
    s_quarks_histograms = create_histograms("s_quarks", quarks_info)
    c_quarks_histograms = create_histograms("c_quarks", quarks_info)

    d_histograms = create_histograms("d", quarks_info)
    dbar_histograms = create_histograms("dbar", quarks_info)
    u_histograms = create_histograms("u", quarks_info)
    ubar_histograms = create_histograms("ubar", quarks_info)
    s_histograms = create_histograms("s", quarks_info)
    sbar_histograms = create_histograms("sbar", quarks_info)
    c_histograms = create_histograms("c", quarks_info)
    cbar_histograms = create_histograms("cbar", quarks_info)


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
    electron_neutrino_histograms = create_histograms("electron_neutrino", leptons_info)
    muon_neutrino_histograms = create_histograms("muon_neutrino", leptons_info)
    tau_neutrino_histograms = create_histograms("tau_neutrino", leptons_info)
    electron_antineutrino_histograms = create_histograms("electron_antineutrino", leptons_info)
    muon_antineutrino_histograms = create_histograms("muon_antineutrino", leptons_info)
    tau_antineutrino_histograms = create_histograms("tau_antineutrino", leptons_info)

    #====================================================================================================================

    # Include all the above histograms in the list 
    histogram_dict = {
        "tt_system": tt_system_histograms,
        "tt_system_2d": tt_system_2d_histograms,
        "bb_system": bb_system_histograms,
        "b_all_system": b_all_system_histograms,
        "prompt_bs": prompt_bs_histograms,
        "bs_from_top": bs_from_top_histograms,
        "quark_system": quark_system_histograms,
        "ud_pairs_system": ud_pairs_system_histograms,
        "sc_pairs_system": sc_pairs_system_histograms,
        "cd_pairs_system": cd_pairs_system_histograms,
        "su_pairs_system": su_pairs_system_histograms,        
        "e_system": e_system_histograms,
        "m_system": m_system_histograms,
        "t_system": t_system_histograms,
        "em_system": em_system_histograms,
        "mm_system": mm_system_histograms,
        "tm_system": tm_system_histograms,
        "ep_system": ep_system_histograms,
        "mp_system": mp_system_histograms,
        "tp_system": tp_system_histograms,

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
        "quarks": quarks_histograms,
        "d_quarks": d_quarks_histograms,
        "u_quarks": u_quarks_histograms,
        "s_quarks": s_quarks_histograms,
        "c_quarks": c_quarks_histograms,
        "d": d_histograms,
        "s": s_histograms,
        "dbar": dbar_histograms,
        "sbar": sbar_histograms,
        "u": u_histograms,
        "c": c_histograms,
        "ubar": ubar_histograms,
        "cbar": cbar_histograms,

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
        "tau_neutrinos": tau_neutrinos_histograms,
        "electron_neutrino": electron_neutrino_histograms,
        "muon_neutrino": muon_neutrino_histograms,
        "tau_neutrino": tau_neutrino_histograms,
        "electron_antineutrino": electron_antineutrino_histograms,
        "muon_antineutrino": muon_antineutrino_histograms,
        "tau_antineutrino": tau_antineutrino_histograms
    }

    #====================================================================================================================

    # Loop through events and fill histograms
    for event in tree:
        temp = 0 # avoid double counting for lepton pairs 
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

            if event.mother1[i] > 2:
                if event.pid[i] == 1 or event.pid[i] == -1 or event.pid[i] == 2 or event.pid[i] == -2 or event.pid[i] == 3 or event.pid[i] == -3 or event.pid[i] == 4 or event.pid[i] == -4:
                    fill_particle_histograms(quarks_histograms, event, i)
                if event.pid[i] == 1 or event.pid[i] == -1: 
                    fill_particle_histograms(d_quarks_histograms, event, i)
                if event.pid[i] == 1:
                    fill_particle_histograms(d_histograms, event, i)
                if event.pid[i] == -1:
                    fill_particle_histograms(dbar_histograms, event, i)
                if event.pid[i] == 2 or event.pid[i] == -2: 
                    fill_particle_histograms(u_quarks_histograms, event, i)
                if event.pid[i] == 2:
                    fill_particle_histograms(u_histograms, event, i)
                if event.pid[i] == -2:
                    fill_particle_histograms(ubar_histograms, event, i)
                if event.pid[i] == 3 or event.pid[i] == -3:
                    fill_particle_histograms(s_quarks_histograms, event, i)
                if event.pid[i] == 3:
                    fill_particle_histograms(s_histograms, event, i)
                if event.pid[i] == -3:
                    fill_particle_histograms(sbar_histograms, event, i)
                if event.pid[i] == 4 or event.pid[i] == -4: 
                    fill_particle_histograms(c_quarks_histograms, event, i)
                if event.pid[i] == 4:
                    fill_particle_histograms(c_histograms, event, i)
                if event.pid[i] == -4:
                    fill_particle_histograms(cbar_histograms, event, i)

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
            if event.pid[i] == 12:
                fill_particle_histograms(electron_neutrino_histograms, event, i)
            if event.pid[i] == -12:
                fill_particle_histograms(electron_antineutrino_histograms, event, i)
            if event.pid[i] == 14:
                fill_particle_histograms(muon_neutrino_histograms, event, i)
            if event.pid[i] == -14:
                fill_particle_histograms(muon_antineutrino_histograms, event, i)
            if event.pid[i] == 16:
                fill_particle_histograms(tau_neutrino_histograms, event, i)
            if event.pid[i] == -16:
                fill_particle_histograms(tau_antineutrino_histograms, event, i)

            if event.status[i] == 1 and event.mother1[i] < 2 and event.pid[i] != 5 and event.pid[i] != -5:
                #if (event.pid[i] <5 and event.pid[i] >-5) or event.pid[i] == 21:
                fill_particle_histograms(radiation_histograms, event, i)

            # Fill histograms for particle systems

            for j in range(len(event.px)):
                if event.pid[i] == 6 and event.pid[j] == -6:
                    fill_system_histograms(tt_system_histograms, event, i, j )
                    fill_system_histograms(tt_system_2d_histograms, event, i, j )
                if event.pid[i] == 5 and event.pid[j] == -5:
                    if i < j:
                        fill_system_histograms(b_all_system_histograms, event, i, j)
                    
                    if event.mother1[i] > 2 and event.mother1[j] > 2:
                        fill_system_histograms(bs_from_top_histograms, event, i, j)
                    elif event.mother1[i] < 3 and event.mother1[j] < 3:
                        fill_system_histograms(prompt_bs_histograms, event, i, j)
                    
                    if event.mother1[i] == event.mother1[j] or event.mother1[i] + event.mother1[j] > 6:
                        fill_system_histograms(bb_system_histograms, event, i, j)
                

                if event.mother1[i] > 2 and event.mother1[j] > 2:              
                    if (event.pid[i] == 1 and event.pid[j] == -2) or (event.pid[i] == -1 and event.pid[j] == 2):
                        fill_system_histograms(ud_pairs_system_histograms, event, i,j)
                    
                    if (event.pid[i] == 3 and event.pid[j] == -4) or (event.pid[i] == -3 and event.pid[j] == 4):
                        fill_system_histograms(sc_pairs_system_histograms, event, i,j)

                    if (event.pid[i] == 3 and event.pid[j] == -2) or (event.pid[i] == -3 and event.pid[j] == 2):
                        fill_system_histograms(su_pairs_system_histograms, event, i,j)
                    
                    if (event.pid[i] == 1 and event.pid[j] == -4) or (event.pid[i] == -1 and event.pid[j] == 4):
                        fill_system_histograms(cd_pairs_system_histograms, event, i,j)

                    if (event.pid[i] == 1 and event.pid[j] == -2) or (event.pid[i] == -1 and event.pid[j] == 2) or (event.pid[i] == 3 and event.pid[j] == -4) or (event.pid[i] == -3 and event.pid[j] == 4) or (event.pid[i] == 3 and event.pid[j] == -2) or (event.pid[i] == -3 and event.pid[j] == 2) or (event.pid[i] == 1 and event.pid[j] == -4) or (event.pid[i] == -1 and event.pid[j] == 4):
                        fill_system_histograms(quark_system_histograms, event, i, j)

                
                if event.pid[i] == 11 or event.pid[i] == -11:
                    if event.pid[j] == 12 or event.pid[j] == -12:
                        fill_system_histograms(e_system_histograms, event, i, j)

                if event.pid[i] == 13 or event.pid[i] == -13:
                    if event.pid[j] == 14 or event.pid[j] == -14:
                        fill_system_histograms(m_system_histograms, event, i, j)

                if event.pid[i] == 15 or event.pid[i] == -15:
                    if event.pid[j] == 16 or event.pid[j] == -16:
                        fill_system_histograms(t_system_histograms, event, i, j)

                if event.pid[i] == 11 and event.pid[j] == -12:
                    fill_system_histograms(em_system_histograms, event, i, j)

                if event.pid[i] == -11 and event.pid[j] == 12:
                    fill_system_histograms(ep_system_histograms, event, i, j)

                if event.pid[i] == 13 and event.pid[j] == -14:
                    fill_system_histograms(mm_system_histograms, event, i, j)

                if event.pid[i] == -13 and event.pid[j] == 14:
                    fill_system_histograms(mp_system_histograms, event, i, j)

                if event.pid[i] == 15 and event.pid[j] == -16:
                    fill_system_histograms(tm_system_histograms, event, i, j)

                if event.pid[i] == -15 and event.pid[j] == 16:
                    fill_system_histograms(tp_system_histograms, event, i, j)
    
    #====================================================================================================================

    # Save all the histograms in a root file and draw them 
    output_root_file = f"{output_directory}.root"
    output_file = ROOT.TFile.Open(output_root_file, "RECREATE")
    output_file.cd()

    for name, histograms in histogram_dict.items():
        overflow_Write(histograms)

    output_file.Close()

    for name, histograms in histogram_dict.items():
        draw(histograms, name, output_directory)  # use the draw_log function if you want log scale in the y axis
    
    #====================================================================================================================

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