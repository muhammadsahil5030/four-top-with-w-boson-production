import ROOT

#---------------- *Open the ROOT files and get the histogram from each file* ----------

fIn1 = ROOT.TFile.Open("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_tttxtx_wminus_fullyHadronic_signal/Events/run_01_decayed_1/unweighted_events.root")

h4 = fIn1.Get("bquark0_pt")
h3 = fIn1.Get("bquark1_pt")
h2 = fIn1.Get("bquark2_pt")
h1 = fIn1.Get("bquark3_pt")
hPt = fIn1.Get("hPt_b")
hEta = fIn1.Get("hEta_b")
hPhi = fIn1.Get("hPhi_b")

#--------------------------- *Scale the histograms by the cross sections* ------------
factor = 1.0
h1.Scale(factor/h1.Integral())
h2.Scale(factor/h2.Integral())
h3.Scale(factor/h3.Integral())
h4.Scale(factor/h4.Integral())
hPt.Scale(factor/hPt.Integral())
hEta.Scale(factor/hEta.Integral())
hPhi.Scale(factor/hPhi.Integral())

#--------------------------- *Create a new canvas and draw the histograms* ------------

#canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
#canvas.Divide(2,2)

#----------------------------------*plot Pt*--------------------------------------------
#canvas.cd(1)
canvas1 = ROOT.TCanvas("canvas1", "Pt b-quarks", 800, 600)
hPt.SetLineWidth(2)
hPt.SetLineColor(ROOT.kBlack)
hPt.SetStats(0)
hPt.SetTitle("B Quark")
hPt.GetXaxis().SetTitle("P_{T} [GeV]")
hPt.GetYaxis().SetTitle("Normalized")
hPt.Draw("hist")

tex1 = ROOT.TLatex(0.153,0.92,"t#bar{t}t#bar{t}W Signal")
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
canvas1.SaveAs("/home/msahil/work/analysis/fourtop_wminus/presentation_1/presentation_images/partons_plots/pt_bquarks.png")

#----------------------------*plot leading subleading Pt*-------------------------
#canvas.cd(2)
canvas2 = ROOT.TCanvas("canvas2", "leading_subleading b-jets", 800, 600)
h1.SetLineColor(ROOT.kRed)
h2.SetLineColor(ROOT.kBlue)
h3.SetLineColor(ROOT.kGreen)
h4.SetLineColor(ROOT.kOrange)

h1.SetLineWidth(2)
h2.SetLineWidth(2)
h3.SetLineWidth(2)
h4.SetLineWidth(2)

h4.SetStats(0)
h4.SetTitle("B Quarks")
h4.GetXaxis().SetTitle("P_{T} [GeV]")
h4.GetYaxis().SetTitle("Normalized")
h4.Draw("hist")
h2.Draw("hist same")
h3.Draw("hist same")
h1.Draw("hist same")

tex3 = ROOT.TLatex(0.153,0.92,"t#bar{t}t#bar{t}W Signal")
tex3.SetNDC()
tex3.SetTextAngle(0)
tex3.SetTextFont(42)
tex3.SetTextSize(0.04)
tex3.SetTextAlign(11)
tex3.Draw()

tex4 = ROOT.TLatex(0.730,0.92,"\sqrt{s} = {13TeV}")
tex4.SetNDC()
tex4.SetTextAngle(0)
tex4.SetTextFont(42)
tex4.SetTextAlign(11)
tex4.SetTextSize(0.04)
tex4.Draw()

legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.85)
legend.AddEntry(h1, "leading", "f")
legend.AddEntry(h2, "2nd leading", "f")
legend.AddEntry(h3, "3rd leading", "f")
legend.AddEntry(h4, "4th leading", "f")
legend.SetBorderSize(0)
legend.Draw()
canvas2.SaveAs("/home/msahil/work/analysis/fourtop_wminus/presentation_1/presentation_images/partons_plots/leading_bquarks.png")

#---------------------------------*plot Eta*-------------------------------------
#canvas.cd(3)
canvas3 = ROOT.TCanvas("canvas3", "b-jets Eta", 800, 600)
hEta.SetLineWidth(2)
hEta.SetLineColor(ROOT.kRed)
hEta.SetStats(0)
hEta.SetTitle("B Quark")
hEta.GetXaxis().SetTitle("#eta")
hEta.GetYaxis().SetTitle("Normalized")
hEta.Draw("E")

tex5 = ROOT.TLatex(0.153,0.92,"t#bar{t}t#bar{t}W Signal")
tex5.SetNDC()
tex5.SetTextAngle(0)
tex5.SetTextFont(42)
tex5.SetTextSize(0.04)
tex5.SetTextAlign(11)
tex5.Draw()

tex6 = ROOT.TLatex(0.730,0.92,"\sqrt{s} = {13TeV}")
tex6.SetNDC()
tex6.SetTextAngle(0)
tex6.SetTextFont(42)
tex6.SetTextAlign(11)
tex6.SetTextSize(0.04)
tex6.Draw()
canvas3.SaveAs("/home/msahil/work/analysis/fourtop_wminus/presentation_1/presentation_images/partons_plots/eta_bquarks.png")

#------------------------------*plot Phi*----------------------------------
#canvas.cd(4)
canvas4 = ROOT.TCanvas("canvas4", "b-jets Phi", 800, 600)
hPhi.SetLineWidth(2)
hPhi.SetLineColor(ROOT.kRed)
hPhi.SetStats(0)
hPhi.SetTitle("B Quark")
hPhi.GetXaxis().SetTitle("#phi")
hPhi.GetYaxis().SetTitle("Normalized")
#h1.GetXaxis().SetRangeUser(0, 100)  # Set x-axis limits for hist1
hPhi.GetYaxis().SetRangeUser(0.015,0.025)  # Set y-axis limits for hist1
hPhi.Draw("E")

tex7 = ROOT.TLatex(0.153,0.92,"t#bar{t}t#bar{t}W Signal")
tex7.SetNDC()
tex7.SetTextAngle(0)
tex7.SetTextFont(42)
tex7.SetTextSize(0.04)
tex7.SetTextAlign(11)
tex7.Draw()

tex8 = ROOT.TLatex(0.730,0.92,"\sqrt{s} = {13TeV}")
tex8.SetNDC()
tex8.SetTextAngle(0)
tex8.SetTextFont(42)
tex8.SetTextAlign(11)
tex8.SetTextSize(0.04)
tex8.Draw()
canvas4.SaveAs("/home/msahil/work/analysis/fourtop_wminus/presentation_1/presentation_images/partons_plots/phi_bquarks.png")

#-------------------------------------- *Show the canvas* -----------------------------

#canvas.Draw()
#canvas.SaveAs("plots/bquarks.png")
input("press enter to exit....")

#-------------------------------------------- *End* -----------------------------------
