

#include "beamHist.h"



BeamHist::BeamHist( int dim_res[EDims], string stats_fname, string histout_fname )
		:_hist( new NArray<int>(dim_res, EDims) ),
		_stats( stats_fname ),
		_inf( EDims ),
		_f( EDims )
		//,_rootout( new RootFSaver< Nd280Element >(string("mytreeout.root"), string("treeout")) )
{
	for( int i = 0; i < _hist->FlatCount(); ++i )
		_hist->Flat( i ) = 0;
	
	_hist->Load( histout_fname ); 
	_hist->DoSum();
	
	_inf[0] = (_stats.GetMax().xnu - _stats.GetMin().xnu)  /  _hist->Count(0);
	_inf[1] = (_stats.GetMax().ynu - _stats.GetMin().ynu)  /  _hist->Count(1);
	_inf[2] = (_stats.GetMax().nnu[0] - _stats.GetMin().nnu[0])  /  _hist->Count(2);
	_inf[3] = (_stats.GetMax().nnu[1] - _stats.GetMin().nnu[1])  /  _hist->Count(3);
	_inf[4] = (_stats.GetMax().nnu[2] - _stats.GetMin().nnu[2])  /  _hist->Count(4);

      _f[0] = _hist->Count(0); 
	for(int i=1; i<EDims; ++i)
	{
		_f[i] = _f[i-1] * _hist->Count(i);
	}

	//_hist = auto_hist->DoUniqueElemsArray();
}



/// create new particle for simulation.
/// each one is defined using the histogram
particle BeamHist::shoot()
{
	particle part( 14, 0.0 );

	int idx[7];
	_hist->RandomIdx(idx);
	
	/// position and time 
	part.r.t = 0;
	part.r.x = _stats.GetMin().xnu + (idx[0]+frandom()) * _inf[0] * 10;
	part.r.y = _stats.GetMin().ynu + (idx[1]+frandom()) * _inf[1] * 10;
	part.r.z = 0;

	/// particle energy and momentum
	part.x = (_stats.GetMin().nnu[0] + (idx[2]+frandom()) * _inf[2]);
	part.y = (_stats.GetMin().nnu[1] + (idx[3]+frandom()) * _inf[3]);
	part.z = (_stats.GetMin().nnu[2] + (idx[4]+frandom()) * _inf[4]);
	part.t = part.momentum();

	/// the young particle
	return part;
}



BeamHist::~BeamHist()
{
	delete _hist;
}
