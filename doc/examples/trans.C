void trans()
{
TFile f("kaskada.root");
treeout->Draw("out[0].t","");
treeout->Draw("out[0].t","@all.size()==1","same");
c1->Print("trans.eps");
}
