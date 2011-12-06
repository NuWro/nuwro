//taken from Justyna, simple file creating simple geometry
// --> Header for unit conversion to CLHEP standards
#include "SystemOfUnits.h"

void make_block()
{
	gSystem->Load("libGeom");

	TGeoManager* gGeoManager = new TGeoManager("SingleBlock","Single Block");

	// new world volume which is vacuum, in which we embed everything else.
	double worldheight = 4000.0*CLHEP::cm; 
	double worldwidth  = 4000.0*CLHEP::cm;
	double worldlength = 3000.0*CLHEP::cm;

	double height = 400.*CLHEP::cm;	// height of the world, was 5000, but the detector was not in the middle
	double width  = 400.*CLHEP::cm;	// width of the world, was 7000, 
	double length = 200.*CLHEP::cm;	// length of the world, was 6000, but the detector was not in the middle

	/* materials */
	TGeoElementTable *ElTable = gGeoManager->GetElementTable();
	TGeoElement *elO  = ElTable->FindElement("Oxygen");
	TGeoElement *elSi = ElTable->FindElement("Silicon");
//	TGeoElement *elO  = new TGeoElement("Oxygen" , "O" , 8., 16.00);	// g/mole
//	TGeoElement *elSi = new TGeoElement("Silicon", "Si",14., 28.0855);	// g/mole

	TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum",0,0,0);

	TGeoMixture *matEarth = new TGeoMixture("Earth", 2, 2.15*CLHEP::gram/CLHEP::cm3);
	matEarth->SetState(TGeoMaterial::kMatStateSolid);
	matEarth->AddElement(elO, 2);
	matEarth->AddElement(elSi,1);

	TGeoMixture *matConcrete = new TGeoMixture("Concrete", 2, 2.3*CLHEP::gram/CLHEP::cm3); // g/cm3
	matConcrete->SetState(TGeoMaterial::kMatStateSolid);
	matConcrete->AddElement(elO, 2);
	matConcrete->AddElement(elSi,1);

	TGeoMedium *Vacuum   = new TGeoMedium("Vacuum",1,matVacuum);
	TGeoMedium *Earth    = new TGeoMedium("Earth",2,matEarth);
	TGeoMedium *Concrete = new TGeoMedium("Concrete",3,matConcrete);

	double halfheight, halflength, halfwidth;

	// World
	TGeoVolume *world = gGeoManager->MakeBox("World",Vacuum,worldwidth/2.,worldheight/2.,worldlength/2.);
	gGeoManager->SetTopVolume(world);
	
	TGeoBBox*    shapeBlock = new TGeoBBox(width/2., height/2., length/2.);
	TGeoVolume*    volBlock = new TGeoVolume("vBlock", shapeBlock, Concrete);
	TGeoTranslation *tBlock = new TGeoTranslation(0., 0., -2000.);
	tBlock->RegisterYourself();
	world->AddNode(volBlock,1,tBlock);

	TGeoBBox*    shapeBlock1 = new TGeoBBox(width/2., height/2., length/2.);
	TGeoVolume*    volBlock1 = new TGeoVolume("vBlock1", shapeBlock1, Concrete);
	TGeoTranslation *tBlock1 = new TGeoTranslation(0., 0., 2000.);
	tBlock1->RegisterYourself();
	world->AddNode(volBlock1,1,tBlock1);

 
	gGeoManager->CloseGeometry();
	gGeoManager->GetGeomPainter();
	gGeoManager->SetTopVisible();
	world->Draw();
//	world->Raytrace();
	gGeoManager->Export("block.root");
}
