import ROOT as rt

f = rt.TFile("flatTuple.root","OPEN")
tree = f.Get("demo/tree")
N = tree.GetEntries()

tree.GetEntry(0)
