import ROOT

#---------------- *Open the ROOT files and get the histogram from each file* ----------

file = ROOT.TFile("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_tttxtx_wminus_fullyHadronic_signal/Events/run_01_decayed_1/TMVAMulticlassification.root", "READ")
directory = file.Get("dataset/Method_BDT/BDT")


#------------------------------------leading Jets-------------------------------------------
signal_h1 = directory.Get("MVA_BDT_Test_Signal_prob_for_Signal")
ttbar_h2 = directory.Get("MVA_BDT_Test_tt~_prob_for_Signal")
ttbarH_h3 = directory.Get("MVA_BDT_Test_tt~H_prob_for_Signal")
ttbarW_h4 = directory.Get("MVA_BDT_Test_tt~W_prob_for_Signal")
ttbarZ_h5 = directory.Get("MVA_BDT_Test_tt~Z_prob_for_Signal")
tttbbar_h6 = directory.Get("MVA_BDT_Test_tt~tb~_prob_for_Signal")
tttWm_h7 = directory.Get("MVA_BDT_Test_tt~tW-_prob_for_Signal")
tttWp_h8 = directory.Get("MVA_BDT_Test_tt~t~w+_prob_for_Signal")

#Scale the histograms by the cross sections

signal_h1.Scale(0.00019884/signal_h1.Integral())
ttbar_h2.Scale(6060.0/ttbar_h2.Integral())
ttbarH_h3.Scale(4.68/ttbarH_h3.Integral())
ttbarW_h4.Scale(1.0122/ttbarW_h4.Integral())
ttbarZ_h5.Scale(3.8640/ttbarZ_h5.Integral())
tttbbar_h6.Scale(0.0006162/tttbbar_h6.Integral())
tttWm_h7.Scale(0.005777999/tttWm_h7.Integral())
tttWp_h8.Scale(0.005751/tttWp_h8.Integral())

#Create a new canvas and draw the histograms
canvas1 = ROOT.TCanvas("canvas1", "Leading Jet Pt", 800, 600)
signal_h1.SetLineColor(ROOT.kBlack)
ttbar_h2.SetFillColor(ROOT.kCyan-10)
ttbarH_h3.SetFillColor(ROOT.kCyan-8)
ttbarW_h4.SetFillColor(ROOT.kCyan-5)
ttbarZ_h5.SetFillColor(ROOT.kYellow-8)
tttbbar_h6.SetFillColor(ROOT.kRed-4)
tttWm_h7.SetFillColor(ROOT.kRed-5)
tttWp_h8.SetFillColor(ROOT.kYellow+1)

signal_h1.SetLineWidth(2)
ttbar_h2.SetStats(0)
ttbar_h2.SetTitle("BDT response")
ttbar_h2.GetXaxis().SetTitle("BDT response for Signal")
ttbar_h2.GetYaxis().SetTitle("(1/N)dN/dx")
#jet0_h2.GetYaxis().SetRangeUser(0,0.25)
#jet0_h2.GetXaxis().SetRangeUser(0, 700)

ttbar_h2.Draw("hist")
ttbarH_h3.Draw("hist same")
ttbarW_h4.Draw("hist same")
ttbarZ_h5.Draw("hist same")
tttbbar_h6.Draw("hist same")
tttWm_h7.Draw("hist same")
tttWp_h8.Draw("hist same")
signal_h1.Draw("hist same")

tex1 = ROOT.TLatex(0.153,0.92,"t#bar{t}t#bar{t}W Signal")
tex1.SetNDC()
tex1.SetTextAngle(0)
tex1.SetTextFont(42)
tex1.SetTextSize(0.04)
tex1.SetTextAlign(11)
tex1.Draw()

tex2 = ROOT.TLatex(0.730,0.92,"\sqrt{S} = {13TeV}")
tex2.SetNDC()
tex2.SetTextAngle(0)
tex2.SetTextFont(42)
tex2.SetTextAlign(11)
tex2.SetTextSize(0.04)
tex2.Draw()

legend1 = ROOT.TLegend(0.7, 0.7, 0.9, 0.85)
legend1.AddEntry(signal_h1, "t#bar{t}t#bar{t}W", "l")
legend1.AddEntry(ttbar_h2, "t#bar{t}", "f")
legend1.AddEntry(ttbarH_h3, "t#bar{t}H", "f")
legend1.AddEntry(ttbarW_h4, "t#bar{t}W", "f")
legend1.AddEntry(ttbarZ_h5, "t#bar{t}Z", "f")
legend1.AddEntry(tttbbar_h6, "t#bar{t}t#bar{b}", "f")
legend1.AddEntry(tttWm_h7, "t#bar{t}tW^{-}", "f")
legend1.AddEntry(tttWp_h8, "t#bar{t}#bar{t}W^{+}", "f")
legend1.SetBorderSize(0)
legend1.Draw()

#canvas1.SetLogy()
canvas1.Draw()
#canvas1.SaveAs("tmva_bdt_variables_plots/bdt_response_signal.png")
input("press enter to exit....")
#-------------------------------------------- *End* -----------------------------------

