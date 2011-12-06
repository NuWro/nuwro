void par()
{
TFile f("ev.root");
event *e;
treeout->SetBranchAddress("e",&e);
treeout->GetEntry(0);
e->par.list();
//c1->Print("trans.eps");
}
