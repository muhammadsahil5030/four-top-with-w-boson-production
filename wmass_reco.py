import ROOT

#---------------- *Open the ROOT files and get the histogram from each file* ----------

fIn1 = ROOT.TFile.Open("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_tttxtx_wminus_fullyHadronic_signal/Events/run_01_decayed_1/unweighted_events.root")

h1 = fIn1.Get("wboson0_mass")
h2 = fIn1.Get("wboson1_mass")
h3 = fIn1.Get("wboson2_mass")
h4 = fIn1.Get("wboson3_mass")
h5 = fIn1.Get("wboson4_mass")

h1.SetLineWidth(2)
h2.SetLineWidth(2)
h3.SetLineWidth(2)
h4.SetLineWidth(2)
h5.SetLineWidth(2)

#--------------------------- *Scale the histograms by the cross sections* ------------
factor = 1.0
h1.Scale(factor/h1.Integral())
h2.Scale(factor/h2.Integral())
h3.Scale(factor/h3.Integral())
h4.Scale(factor/h4.Integral())
h5.Scale(factor/h5.Integral())

#--------------------------- *Create a new canvas and draw the histograms* ------------

canvas = ROOT.TCanvas("canvas", "W boson Mass", 800, 600)
h1.SetLineColor(ROOT.kRed)
h2.SetLineColor(ROOT.kBlue)
h3.SetLineColor(ROOT.kGreen)
h4.SetLineColor(ROOT.kBlack)
h5.SetLineColor(ROOT.kOrange)
h1.SetStats(0)

#------------------------------*Normalizing Histogram*-----------------------------------------
factor = 1.0
h1.Scale(factor/h1.Integral())
h2.Scale(factor/h2.Integral())
h3.Scale(factor/h3.Integral())
h4.Scale(factor/h4.Integral())
h5.Scale(factor/h5.Integral())

#h1.GetXaxis().SetRangeUser(0, 100)  # Set x-axis limits for hist1
#h1.GetYaxis().SetRangeUser(0,2600)  # Set y-axis limits for hist1

h1.Draw("hist")
h2.Draw("hist same")
h3.Draw("hist same")
h4.Draw("hist same")
h5.Draw("hist same")

#----------------------- *Add the title and axis labels* ---------------------

h1.SetTitle("W Boson mass")
h1.GetXaxis().SetTitle("m_{w} [GeV]")
h1.GetYaxis().SetTitle("Normalized")

tex1 = ROOT.TLatex(0.153,0.92,"t#bar{t}t#bar{t}W signal")
tex1.SetNDC()
tex1.SetTextAngle(0)
tex1.SetTextFont(42)
tex1.SetTextSize(0.04)
tex1.SetTextAlign(11)
tex1.Draw()

tex2 = ROOT.TLatex(0.730,0.92,"\sqrt{s} = {13TeV}")
tex2.SetNDC()
tex2.SetTextAngle(0)
tex2.SetTextFont(42)
tex2.SetTextAlign(11)
tex2.SetTextSize(0.04)
tex2.Draw()


#--------------------------------------- *Add a legend* -------------------------------

legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.85)
legend.AddEntry(h1, "w1 mass", "f")
legend.AddEntry(h2, "w2 mass", "f")
legend.AddEntry(h3, "w3 mass", "f")
legend.AddEntry(h4, "w4 mass", "f")
legend.AddEntry(h5, "w5 mass", "f")
legend.SetBorderSize(0)
legend.Draw()

#-------------------------------------- *Show the canvas* -----------------------------
#canvas.SetLogy()
canvas.Draw()
canvas.SaveAs("/home/msahil/work/analysis/fourtop_wminus/presentation_1/presentation_images/partons_plots/Wboson_Mass.png")
input("press enter to exit....")

#-------------------------------------------- *End* -----------------------------------
