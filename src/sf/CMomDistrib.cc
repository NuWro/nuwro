#include "CMomDistrib.h"

#include "CMomDistrib_C12_Ben.h"
#include "CMomDistrib_O16_Ben.h"
#include "CMomDistrib_O16_GCo.h"
#include "CMomDistrib_O16_CdA.h"
#include "CMomDistrib_Ca40_GCo.h"
#include "CMomDistrib_Ca40_CdA.h"
#include "CMomDistrib_Ca48_GCo.h"

#include "CMomDistrib_Fe56_Ben.h"

CMomDistrib* createDistrib(const MomDistribs i_md, const IsospinOfSF i_isospin)
{
	switch (i_md)
	{
		case md_C12_Ben:
			return new CMomDistrib_C12_Ben();

		case md_O16_Ben:
			return new CMomDistrib_O16_Ben();

		case md_O16_GCo:
			return new CMomDistrib_O16_GCo();
			
		case md_O16_CdA:
			return new CMomDistrib_O16_CdA();

		case md_Ca40_GCo:
			return new CMomDistrib_Ca40_GCo(i_isospin);

		case md_Ca40_CdA:
			return new CMomDistrib_Ca40_CdA();

		case md_Ca48_GCo:
			return new CMomDistrib_Ca48_GCo(i_isospin);
				
		case md_Fe56_Ben:
			return new CMomDistrib_Fe56_Ben();

		default:
			return 0;
	}
}

