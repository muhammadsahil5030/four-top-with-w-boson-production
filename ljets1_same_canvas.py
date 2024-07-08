import ROOT

#---------------- *Open the ROOT files and get the histogram from each file* ----------

fIn1 = ROOT.TFile.Open("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_tttxtx_wminus_fullyHadronic_signal/Events/run_01_decayed_1/signal_bdt_variables.root")
fIn2 = ROOT.TFile.Open("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttx_bg_tttxtxwm/Events/run_02_decayed_1/ttbar_bg_bdt_variables.root")
fIn3 = ROOT.TFile.Open("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxwm_bg_tttxtxwm/Events/run_01_decayed_1/ttbarwm_bg_bdt_variables.root")
fIn4 = ROOT.TFile.Open("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxz_bg_tttxtxwm/Events/run_01_decayed_1/ttbarz_bg_bdt_variables.root")
fIn5 = ROOT.TFile.Open("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxH_bg_tttxtxwm/Events/run_01_decayed_1/ttbarH_bg_bdt_variables.root")

h1 = fIn1.Get("jet_pt_1")
h2 = fIn2.Get("jet_pt_1")
h3 = fIn3.Get("jet_pt_1")
h4 = fIn4.Get("jet_pt_1")
h5 = fIn5.Get("jet_pt_1")

#--------------------------- *Scale the histograms by the cross sections* ------------
factor = 1.0
h1.Scale(factor/h1.Integral())
h2.Scale(factor/h2.Integral())
h3.Scale(factor/h3.Integral())
h4.Scale(factor/h4.Integral())
h5.Scale(factor/h5.Integral())

#--------------------------- *Create a new canvas and draw the histograms* ------------

canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)


h1.SetLineColor(ROOT.kRed)
h2.SetFillColor(ROOT.kCyan-10)
h3.SetFillColor(ROOT.kCyan-8)
h4.SetFillColor(ROOT.kCyan-5)
h5.SetFillColor(ROOT.kYellow-8)
h1.SetLineWidth(2)
#h2.SetLineWidth(2)
#h3.SetLineWidth(2)
#h4.SetLineWidth(2)
#h5.SetLineWidth(2)
h2.SetStats(0)
h2.SetTitle("2nd Leading Jet P_{T}")
h2.GetXaxis().SetTitle("Pt_{j} [GeV]")
h2.GetYaxis().SetTitle("Normalized")
h2.GetYaxis().SetRangeUser(0,0.30)

h2.Draw("hist")
h3.Draw("hist same")
h4.Draw("hist same")
h5.Draw("hist same")
h1.Draw("hist same")

tex3 = ROOT.TLatex(0.153,0.92,"t#bar{t}t#bar{t}w- Signal")
tex3.SetNDC()
tex3.SetTextAngle(0)
tex3.SetTextFont(42)
tex3.SetTextSize(0.04)
tex3.SetTextAlign(11)
tex3.Draw()

tex4 = ROOT.TLatex(0.730,0.92,"\sqrt{S} = {13TeV}")
tex4.SetNDC()
tex4.SetTextAngle(0)
tex4.SetTextFont(42)
tex4.SetTextAlign(11)
tex4.SetTextSize(0.04)
tex4.Draw()

legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.85)
legend.AddEntry(h1, "t#bar{t}t#bar{t}w", "l")
legend.AddEntry(h2, "t#bar{t}", "f")
legend.AddEntry(h3, "t#bar{t}W", "f")
legend.AddEntry(h4, "t#bar{t}Z", "f")
legend.AddEntry(h5, "t#bar{t}H", "f")
legend.SetBorderSize(0)
legend.Draw()

#-------------------------------------- *Show the canvas* -----------------------------
#canvas.SetLogy()
canvas.Draw()
#canvas.SaveAs("plots/uquarks.png")
input("press enter to exit....")

#-------------------------------------------- *End* -----------------------------------
