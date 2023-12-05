import ROOT
import argparse
import pathlib
import sys
import os

def create_canvas(hist, canvas_title, output_file):
    canvas = ROOT.TCanvas(canvas_title, canvas_title, 800, 600)

    # For 2D histograms, use "colz" option
    if "TH2" in hist.ClassName():
        hist.Draw("colz")
    else:
        hist.Draw()

    canvas.SaveAs(output_file)

def create_plots(input_root_file, output_directory):
    
    ROOT.gROOT.SetBatch(True) # not show canvas when draw

    # Open the ROOT file and check if it is valid
    root_file = ROOT.TFile.Open(input_root_file)

    if not root_file or root_file.IsZombie():
        print(f"Error: Unable to open or read file: {input_root_file}")
        sys.exit(1)

    # Access the TTree and check if it is valid 
    tree = root_file.Get("events")

    if not tree:
        print("Error: TTree 'events' not found in the ROOT file.")
        sys.exit(1)

    # Create histograms for selected branches
    hist_numParticles = ROOT.TH1F("hist_numParticles", "Number of Particles", 100, 0, 100)
    hist_eventweight = ROOT.TH1F("hist_eventweight", "Event Weight", 100, -1, 1)
    hist_px = ROOT.TH1F("hist_px", "Particle Momentum (px)", 100, -500, 500)
    hist_py = ROOT.TH1F("hist_py", "Particle Momentum (py)", 100, -500, 500)
    hist_pz = ROOT.TH1F("hist_pz", "Particle Momentum (pz)", 100, -500, 500)
    hist_pid = ROOT.TH1F("hist_pid", "Particle ID (pid)", 100, -50, 50)
    hist_mother1 = ROOT.TH1F("hist_mother1", "Mother ID (mother1)", 100, -50, 50)
    hist_mother2 = ROOT.TH1F("hist_mother2", "Mother ID (mother2)", 100, -50, 50)
    hist_energy = ROOT.TH1F("hist_energy", "Particle Energy", 100, 0, 1000)
    hist_mass = ROOT.TH1F("hist_mass", "Particle Mass", 100, 0, 200)
    hist_pt = ROOT.TH1F("hist_pt", "Particle Transverse Momentum (pT)", 100, 0, 500)

    # Create 2D histograms
    hist_px_py = ROOT.TH2F("hist_px_py", "Particle Momentum (px) vs (py)", 100, -500, 500, 100, -500, 500)
    hist_px_pz = ROOT.TH2F("hist_px_pz", "Particle Momentum (px) vs (pz)", 100, -500, 500, 100, -500, 500)

    # Loop through the TTree and fill histograms
    for event in tree:
        hist_numParticles.Fill(event.numParticles)
        hist_eventweight.Fill(event.eventweight)
        for i in range(len(event.px)):
            hist_px.Fill(event.px[i])
            hist_py.Fill(event.py[i])
            hist_pz.Fill(event.pz[i])
            hist_pid.Fill(event.pid[i])
            hist_mother1.Fill(event.mother1[i])
            hist_mother2.Fill(event.mother2[i])
            hist_energy.Fill(event.energy[i])
            hist_mass.Fill(event.mass[i])

            # Calculate pT
            pt = (event.px[i]**2 + event.py[i]**2)**0.5
            hist_pt.Fill(pt)

            # Fill 2D histograms
            hist_px_py.Fill(event.px[i], event.py[i])
            hist_px_pz.Fill(event.px[i], event.pz[i])

    # Create canvas and draw the histograms
    canvas = ROOT.TCanvas("canvas", "Plots", 800, 600)
    canvas.Divide(3, 3)

    canvas.cd(1)
    hist_px.Draw()
    canvas.cd(2)
    hist_py.Draw()
    canvas.cd(3)
    hist_pz.Draw()
    canvas.cd(4)
    hist_pid.Draw()
    canvas.cd(5)
    hist_mother1.Draw()
    canvas.cd(6)
    hist_mother2.Draw()
    canvas.cd(7)
    hist_energy.Draw()
    canvas.cd(8)
    hist_mass.Draw()
    canvas.cd(9)
    hist_pt.Draw()

    # Save the canvas as PNG and as PDF
    output_png = f"{output_directory}/all.png"
    canvas.SaveAs(output_png)
    output_pdf = f"{output_directory}/all.pdf"
    canvas.SaveAs(output_pdf)

    # Create and save individual histograms
    create_canvas(hist_px, "canvas_px", f"{output_directory}/px.png")
    create_canvas(hist_py, "canvas_py", f"{output_directory}/py.png")
    create_canvas(hist_pz, "canvas_pz", f"{output_directory}/pz.png")
    create_canvas(hist_pid, "canvas_pid", f"{output_directory}/pid.png")
    create_canvas(hist_mother1, "canvas_mother1", f"{output_directory}/mother1.png")
    create_canvas(hist_mother2, "canvas_mother2", f"{output_directory}/mother2.png")
    create_canvas(hist_energy, "canvas_energy", f"{output_directory}/energy.png")
    create_canvas(hist_mass, "canvas_mass", f"{output_directory}/mass.png")
    create_canvas(hist_pt, "canvas_pt", f"{output_directory}/pt.png")
    create_canvas(hist_px_py, "canvas_px_py", f"{output_directory}/px_py.png")
    create_canvas(hist_px_pz, "canvas_px_pz", f"{output_directory}/px_pz.png")

    # Close the ROOT file
    root_file.Close()

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, help='Input file path and Name')
    parser.add_argument("--output", type=str, help='Output directory path')
    args = parser.parse_args()

    # Check if the input file has a .root extension
    if pathlib.PurePosixPath(args.input).suffix.lower() != ".root":
        print('Error: Only .root files are accepted as input. Exiting ...')
        sys.exit(1)

    input_root_file = args.input
    output_directory = args.output

    # Check if the output directory exists, create it if not
    if output_directory and not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    create_plots(input_root_file, output_directory)
