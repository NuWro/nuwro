./xy-normalize <BNL_flux_precise.dat  >BNL_flux_precise.norm
./xy-normalize <BNL_flux.dat  >BNL_flux.norm
./xy2hist <BNL_flux_precise.norm  >BNL_flux_precise.norm.hist
./xy2hist <BNL_flux.norm  >BNL_flux.norm.hist
./hist2xy <BNL_flux_precise.norm.hist  >BNL_flux_precise.norm.xy
./hist2xy <BNL_flux.norm.hist  >BNL_flux.norm.xy


xmgrace    BNL_flux.norm  BNL_flux_precise.norm  BNL_flux.norm.xy  BNL_flux_precise.norm.xy
