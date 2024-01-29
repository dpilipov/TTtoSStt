import ROOT
f = ROOT.TFile.Open('rootfiles/THselection_HT0_StoAA_16.root','read')
# get the histos
h1 = f.Get("MtpvMs_TvsQCD_mvaID_SR_pass__nominal")
h2 = f.Get("MtpvMs_TvsQCD_mvaID_SR_fail__nominal")
h3 = f.Get("MtpvMs_TvsQCD_mvaID_CR_pass__nominal")
h4 = f.Get("MtpvMs_TvsQCD_mvaID_CR_fail__nominal")

# make a canvas
c = ROOT.TCanvas('c','c')
c.cd()
c.Print('outputSelection.pdf[')  # this tells the canvas to open file for multi-page printing
c.Clear()

# loop over the histos u want to draw
for h in [h1, h2,h3,h4]:
    c.Clear()   # clear any previous draws from the canvas
    h.SetStats(0)
    h.SetMaximum(3.00)
    h.Draw("colz")  # draw the histo with colz
    c.Print("outputSelection.pdf")  # save it to the pdf

c.Print("outputSelection.pdf]")  # finished, close the multipage file
