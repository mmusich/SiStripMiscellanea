import json, os, sys
from pprint import pprint
from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TH1F, TLegend, TString
from ROOT import gROOT
from ROOT import gStyle

################################################################
def makeNicePlot(hist):
################################################################
    hist.SetLineWidth(2);
    hist.GetXaxis().SetNdivisions(505);
    hist.GetXaxis().CenterTitle(True);
    hist.GetYaxis().CenterTitle(True);
    hist.GetXaxis().SetTitleFont(42); 
    hist.GetYaxis().SetTitleFont(42);  
    hist.GetXaxis().SetTitleSize(0.05);
    hist.GetYaxis().SetTitleSize(0.05);
    hist.GetXaxis().SetTitleOffset(0.9);
    hist.GetYaxis().SetTitleOffset(1.2);
    hist.GetXaxis().SetLabelFont(42);
    hist.GetXaxis().SetLabelSize(0.03);
    hist.GetYaxis().SetLabelFont(42);
 
################################################################
if len(sys.argv) < 3:
  print "Error, expecting 2 arguments (two input json files)!"
  sys.exit()

oldInputFile     = sys.argv[2]
newInputFile     = sys.argv[1]

json_data=open(newInputFile).read()
data = json.loads(json_data)

json_data2=open(oldInputFile).read()
data2 = json.loads(json_data2)

treeIndex=0
for tree in (data[0]["Trees"]):
    print "======>",tree["Name"],tree["Size"]

    for tt in (data2[0]["Trees"]):
        if(tt["Name"]==tree["Name"]):
            tree2 = tt

    c1 = TCanvas('c1'+tree["Name"], 'Histogram Drawing Options',1000,800)
    c1.SetBottomMargin(0.32)
    c1.SetTopMargin(0.08)
    c1.SetRightMargin(0.05)
    c1.SetLeftMargin(0.12)

    c1.cd()
    h1 = TH1F("h1"+tree["Name"],tree["Name"].replace("/"," ")+";;#LT bytes/entry #GT",len(tree["Branches"]),0.,len(tree["Branches"]))
    h2 = TH1F("h2"+tree["Name"],tree["Name"].replace("/"," ")+";;#LT bytes/entry #GT",len(tree["Branches"]),0.,len(tree["Branches"]))
    h1.SetStats(False)
    h2.SetStats(False)

    bincounter=0;

    print "first tree",len(tree["Branches"])," second tree",len(tree2["Branches"])

    if(tree["Entries"]>0.):
        for entry in tree["Branches"]:

            the2treeBin=-1
            for bin in range(0,len(tree2["Branches"])):
                if (tree2["Branches"][bin]["Name"]==entry["Name"]):
                    the2treeBin=bin
                    break

            if(the2treeBin==-1):
                print "=======>   ",entry["Name"]," not found in the second tree" 

            entry2 = tree2["Branches"][the2treeBin]
            
            bincounter+=1

            if(entry["Entries"]>0.): 
                print "tree1:",entry["Name"],entry["Size"],entry["Entries"],(float(entry["Size"])/float(entry["Entries"]))
                if(the2treeBin>=0):
                    print "tree2:",entry2["Name"],entry2["Size"],entry2["Entries"],(entry2["Size"]/entry2["Entries"])

                h1.SetBinContent(bincounter,float(entry["Size"])/float(entry["Entries"]))
                h1.GetXaxis().SetBinLabel(bincounter,entry["Name"])
                h1.GetXaxis().LabelsOption("v");

                if(the2treeBin>=0):
                    h2.SetBinContent(bincounter,float(entry2["Size"])/float(entry2["Entries"]))
                h2.GetXaxis().SetBinLabel(bincounter,entry["Name"])
                h2.GetXaxis().LabelsOption("v");

            else:
                print entry["Name"],entry["Size"],entry["Entries"]

    treeIndex+=1

    gStyle.SetPaintTextFormat("4.1f B");

    makeNicePlot(h1)
    h1.SetLineColor(4);
    h1.SetFillColor(4);
    h1.SetLineStyle(9);
    h1.SetMarkerColor(4);
    
    h1.SetBarWidth(0.45);
    h1.SetBarOffset(0.1);
    h1.GetYaxis().SetRangeUser(0.,h1.GetMaximum()*1.30);

    h1.Draw("bar2");
    h1.Draw("TEXTsame");

    makeNicePlot(h2)
    h2.SetLineColor(2);
    h2.SetFillColor(2);
    h2.SetLineStyle(9);
    h2.SetMarkerColor(2);
    
    h2.SetBarWidth(0.4);
    h2.SetBarOffset(0.55);
    h2.Draw("bar2same");
    #h2.Draw("TEXT45same");

    # Draw the legend
    infoBox = TLegend(0.65, 0.82, 0.95, 0.92, "");
    if(tree["Entries"]>0.): 
        infoBox.AddEntry(h1,"New LA: %.5s kB/entry" % str(float(tree["Size"])/float(1024*tree["Entries"])),"F");
        infoBox.AddEntry(h2,"original:  %.5s kB/entry" % str(float(tree2["Size"])/float(1024*tree2["Entries"])),"F");
    infoBox.SetShadowColor(0);  # 0 = transparent
    infoBox.SetFillColor(0); 
    infoBox.Draw("same");

    #h1.Draw("BAR")
    c1.Update()
    c1.SaveAs(tree["Name"].replace('/','_')+".png")
    c1.SaveAs(tree["Name"].replace('/','_')+".pdf")


#pprint(data)
