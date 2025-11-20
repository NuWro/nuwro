#ifndef C_SPECTRAL_FUNCTIONS_H
#define C_SPECTRAL_FUNCTIONS_H

#include <fstream>
#include <iomanip>
#include <limits>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include "GConstants.h"
#include "pdg.h"
#include "CLorentzDistrib2.h"
#include "CInterpolatedData.h"
#include "CInterpolatedData2D.h"

/// Class: CSpectralFunctions
/// Encapsulates proton/neutron spectral functions and momentum distributions for different nuclei
class CSpectralFunctions
{
    public:
    /// Constructor
        /// i_target           = target nucleus enum
        /// i_switchSeparation = use separated (MF+corr) SF (only for Ar40)
        /// i_switchFSI        = load FSI-related inputs
        CSpectralFunctions(TargetNucleus i_target, const bool i_switchSeparation, const bool i_switchFSI)
        :m_target( i_target ),
        m_switchSeparation( m_target == Ar40_SF ? i_switchSeparation : false ),
        m_switchFSI( i_switchFSI )
        {
            std::ostringstream pHoleSFName;///proton hole spectral function
            std::ostringstream nHoleSFName;///neutron hole spectral function

            std::ostringstream pHoleMFSFName;///proton hole mean-field spectral function
            std::ostringstream nHoleMFSFName;///neutron hole mean-field spectral function

            std::ostringstream pHoleCorrSFName;///proton hole correlated spectral function
            std::ostringstream nHoleCorrSFName;///neutron hole correlated spectral function

            std::ostringstream transparencyNameA, transparencyNameB;
            std::ostringstream realOPName;
            std::ostringstream foldingFName;

            double foldingFNormalization(0.0);

            // thresholds for SRC additional nucleon generation
            m_threshold_E2_p = 0.0;
            m_threshold_E2_n = 0.0;

            m_switchIsospinAsymmetry = false;

            /// Target-specific initialization
            switch ( m_target )
            {
               case C12_newSF:
                    m_Z = carbonZ;
                    m_N = carbonN;
                    pHoleSFName<<"data/sf/pke_12C_new.dat";
                    // Lorentz distribution used for MC sampling of momenta
                    m_pLorentz = new CLorentzDistrib(1.23, 149.5*MeV, 82.0*MeV, 0.0, 800*MeV);
                    transparencyNameA<<"data/sf/transp_12C.dat";
                    transparencyNameB<<"data/sf/transp_12C_MC.dat";
                    realOPName<<"data/sf/realOP_12C_EDAI.dat";
                    foldingFName<<"data/sf/foldingF_12C.dat";
                    m_targetMass = carbon12Mass;
                    m_neutronEnergyShift = carbon11Mass - boron11Mass + nMass - pMass - 0.510998910*MeV;
                    m_CoulombAvEnergy  = 3.45273*MeV;
                    m_CoulombMaxEnergy = 4.58276*MeV;
                    m_pFermiMomentum = 211.0*MeV;
                    m_pNormalization = 5.99978;
                    foldingFNormalization = 0.120103;
                    break;

                case C12_SF:
                    m_Z = carbonZ;
                    m_N = carbonN;
                    pHoleSFName<<"data/sf/pke_12C.dat";
                    m_pLorentz = new CLorentzDistrib(1.24, 150.0*MeV, 82.0*MeV, 0.0, 800*MeV);
                    transparencyNameA<<"data/sf/transp_12C.dat";
                    transparencyNameB<<"data/sf/transp_12C_MC.dat";
                    realOPName<<"data/sf/realOP_12C_EDAI.dat";
                    foldingFName<<"data/sf/foldingF_12C.dat";
                    m_targetMass = carbon12Mass;
                    m_neutronEnergyShift = carbon11Mass - boron11Mass + nMass - pMass - 0.510998910*MeV;
                    m_CoulombAvEnergy  = 3.45273*MeV;
                    m_CoulombMaxEnergy = 4.58276*MeV;
                    m_pFermiMomentum = 211.0*MeV;
                    m_pNormalization = 6.02599;
                    foldingFNormalization = 0.120103;
                    break;

                case O16_SF:
                    m_Z = oxygenZ;
                    m_N = oxygenN;
                    pHoleSFName<<"data/sf/pke_16O.dat";
                    transparencyNameA<<"data/sf/transp_16O.dat";
                    transparencyNameB<<"data/sf/transp_16O_MC.dat";
                    realOPName<<"data/sf/realOP_16O_EDAI.dat";
                    foldingFName<<"data/sf/foldingF_12C.dat";
                    m_pLorentz = new CLorentzDistrib(1.28, 149.0*MeV, 87.0*MeV, 0.0, 800*MeV);
                    m_targetMass = oxygen16Mass;
                    m_neutronEnergyShift = oxygen15Mass - nitrogen15Mass + nMass - pMass - 0.510998910*MeV;
                    m_CoulombAvEnergy = 4.17276*MeV;
                    m_CoulombMaxEnergy = 5.48945*MeV;
                    m_pFermiMomentum = 212.0*MeV;
                    m_pNormalization = 16.0;
                    foldingFNormalization = 0.120103;
                    break;

                case Ar40_SF:
                    m_switchIsospinAsymmetry = true;
                    m_Z = argonZ;
                    m_N = argonN;
                    ///total
                    pHoleSFName<<"data/sf/pke_40Ar_exp.dat";
                    nHoleSFName<<"data/sf/pke_48Ti_exp.dat";
                    m_pLorentz = new CLorentzDistrib(1.21, 162.5*MeV, 85.0*MeV, 0.0, 800*MeV);
                    m_nLorentz = new CLorentzDistrib(1.21, 167.5*MeV, 83.5*MeV, 0.0, 800*MeV);
                    m_pNormalization = 17.0217;///17.0217 is the fraction of 18.0 that was actually measured (94.57%)
                    m_nNormalization = 20.3208;///20.3208 is the fraction of 22.0 that was actually measured (92.37%)
                    ///MF
                    pHoleMFSFName<<"data/sf/pke_40ArP_MF_exp.dat";
                    nHoleMFSFName<<"data/sf/pke_48TiP_MF_exp.dat";
                    m_pLorentz_MF = new CLorentzDistrib(1.226, 162.0*MeV, 80.0*MeV, 0.0, 400*MeV);
                    m_nLorentz_MF = new CLorentzDistrib(1.226, 167.5*MeV, 78.0*MeV, 0.0, 400*MeV);
                    m_pNormalization_MF = 17.0217;///17.0217 is the fraction of 18.0 that was actually measured (94.57%)
                    m_nNormalization_MF = 20.3208;///20.3208 is the fraction of 22.0 that was actually measured (92.37%)
                    ///corr
                    pHoleCorrSFName<<"data/sf/pke_40ArP_corr_exp.dat";
                    nHoleCorrSFName<<"data/sf/pke_48TiP_corr_exp.dat";
                    m_pLorentz_corr = new CLorentzDistrib(1.376, 160.0*MeV, 313.0*MeV, 0.0, 1000*MeV);
                    m_nLorentz_corr = new CLorentzDistrib(1.376, 160.0*MeV, 313.0*MeV, 0.0, 1000*MeV);
                    //pHoleCorrSFName<<"data/sf/pke_40ArP_corr_exp_test.dat";
                    //nHoleCorrSFName<<"data/sf/pke_48TiP_corr_exp_test.dat";
                    //m_pLorentz_corr = new CLorentzDistrib(1.347, 160.0*MeV, 310.0*MeV, 0.0, 800*MeV);
                    //m_nLorentz_corr = new CLorentzDistrib(1.347, 160.0*MeV, 310.0*MeV, 0.0, 800*MeV);
                    m_pNormalization_corr = 17.0217;///17.0217 is the fraction of 18.0 that was actually measured (94.57%)
                    m_nNormalization_corr = 20.3208;///20.3208 is the fraction of 22.0 that was actually measured (92.37%)
                    m_threshold_E2_p = 24.100*MeV;
                    m_threshold_E2_n = 25.199*MeV;
                    ///
                    transparencyNameA<<"data/sf/transp_40Ar.dat";
                    transparencyNameB<<"data/sf/transp_40Ar_MC.dat";
                    realOPName<<"data/sf/realOP_40Ar_EDAD_fit3.dat";
                    foldingFName<<"data/sf/foldingF_12C.dat";
                    m_targetMass = argon40Mass;
                    m_neutronEnergyShift = (argon39Mass + nMass - argon40Mass) - (scandium47Mass + pMass + 0.510998910*MeV - titanium48Mass);
                    m_CoulombAvEnergy = 7.26739*MeV;
                    m_CoulombMaxEnergy = 9.55497*MeV;
                    m_pFermiMomentum = 220.0*MeV;
                    m_nFermiMomentum = 225.0*MeV;
                    foldingFNormalization = 0.120103;
                    break;

                case Fe56_SF:
                    m_Z = ironZ;
                    m_N = ironN;
                    pHoleSFName<<"data/sf/pke_56Fe.dat";
                    transparencyNameA<<"data/sf/transp_56Fe.dat";
                    transparencyNameB<<"data/sf/transp_56Fe_MC.dat";
                    realOPName<<"data/sf/realOP_56Fe_EDAD_fit1.dat";
                    foldingFName<<"data/sf/foldingF_12C.dat";
                    m_pLorentz = new CLorentzDistrib(1.23, 176.5*MeV, 68.0*MeV, 0.0, 800*MeV);
                    m_targetMass = iron56Mass;
                    m_neutronEnergyShift = iron55Mass - manganese55Mass + nMass - pMass - 0.510998910*MeV;
                    m_CoulombAvEnergy = 9.58122*MeV;
                    m_CoulombMaxEnergy = 12.42035*MeV;
                    m_pFermiMomentum = 226.0*MeV;
                    m_pNormalization = 26.0964;
                    foldingFNormalization = 0.120103;
                    break;

                default:
                    std::cout<<"No such target implemented"<<std::endl;
                    return;
            };

            // If not separating MF and corr parts: load combined SFs
            if ( not m_switchSeparation )
            {
                m_pNormalization /= pow(MeV,4);
                m_nNormalization /= pow(MeV,4);

                // proton SF table (p,E) -> S(p,E)
                m_pHoleSF =  new CInterpolatedData2D( pHoleSFName.str(), MeV, MeV, 1.0/m_pNormalization, false );
                if (!m_pHoleSF->isValid())
                   {
                   throw std::runtime_error("Failed to initialize m_pHoleSF.");
                   }
                // momentum distribution (from integrating SF over E)
                m_pMomDistrib = new CInterpolatedData( m_pHoleSF->get_xDistribution(), m_pHoleSF->get_xRes(), MeV, 4*Pi, m_pHoleSF->get_xStart(), m_pHoleSF->get_xStep(), m_pHoleSF->get_xStop(), false);

                if ( m_switchIsospinAsymmetry )
                {
                    m_nHoleSF =  new CInterpolatedData2D( nHoleSFName.str(), MeV, MeV, 1.0/m_nNormalization, false );
                    if (!m_nHoleSF->isValid())
                        {
                        throw std::runtime_error("Failed to initialize m_nHoleSF.");
                        }
                    m_nMomDistrib = new CInterpolatedData( m_nHoleSF->get_xDistribution(), m_nHoleSF->get_xRes(), MeV, 4*Pi, m_nHoleSF->get_xStart(), m_nHoleSF->get_xStep(), m_nHoleSF->get_xStop(), false);
                }

            }

            else
            {
                // Separated MF and correlated parts
                m_pNormalization_MF /= pow(MeV,4);
                m_pHoleSF_MF =  new CInterpolatedData2D( pHoleMFSFName.str(), MeV, MeV, 1.0/m_pNormalization_MF, false );
                m_pMomDistrib_MF = new CInterpolatedData( m_pHoleSF_MF->get_xDistribution(), m_pHoleSF_MF->get_xRes(), MeV, 4*Pi, m_pHoleSF_MF->get_xStart(), m_pHoleSF_MF->get_xStep(), m_pHoleSF_MF->get_xStop(), false);

                m_pNormalization_corr /= pow(MeV,4);
                m_pHoleSF_corr =  new CInterpolatedData2D( pHoleCorrSFName.str(), MeV, MeV, 1.0/m_pNormalization_corr, false );
                m_pMomDistrib_corr = new CInterpolatedData( m_pHoleSF_corr->get_xDistribution(), m_pHoleSF_corr->get_xRes(), MeV, 4*Pi, m_pHoleSF_corr->get_xStart(), m_pHoleSF_corr->get_xStep(), m_pHoleSF_corr->get_xStop(), false);

                if ( m_switchIsospinAsymmetry )
                {
                    m_nNormalization_MF /= pow(MeV,4);
                    m_nHoleSF_MF =  new CInterpolatedData2D( nHoleMFSFName.str(), MeV, MeV, 1.0/m_nNormalization_MF, false );
                    m_nMomDistrib_MF = new CInterpolatedData( m_nHoleSF_MF->get_xDistribution(), m_nHoleSF_MF->get_xRes(), MeV, 4*Pi, m_nHoleSF_MF->get_xStart(), m_nHoleSF_MF->get_xStep(), m_nHoleSF_MF->get_xStop(), false);

                    m_nNormalization_corr /= pow(MeV,4);
                    m_nHoleSF_corr =  new CInterpolatedData2D( nHoleCorrSFName.str(), MeV, MeV, 1.0/m_nNormalization_corr, false );
                    m_nMomDistrib_corr = new CInterpolatedData( m_nHoleSF_corr->get_xDistribution(), m_nHoleSF_corr->get_xRes(), MeV, 4*Pi, m_nHoleSF_corr->get_xStart(), m_nHoleSF_corr->get_xStep(), m_nHoleSF_corr->get_xStop(), false);
                }

            }

            // If FSI enabled, load transparency, real OP, folding function
            if ( m_switchFSI )
            {
                m_transparency[0] = new CInterpolatedData(transparencyNameA.str(), GeV, 1.0, false);
                m_transparency[1] = new CInterpolatedData(transparencyNameB.str(), GeV, 1.0, false);
                m_realOP = new CInterpolatedData( realOPName.str(), MeV, MeV, false );
                m_foldingF = new CInterpolatedData( foldingFName.str(), MeV, 1.0/MeV/foldingFNormalization, false );
            }

        };


        /// Destructor
        ~CSpectralFunctions()
        {
            if ( not m_switchSeparation )
            {
                delete m_pLorentz;
                delete m_pHoleSF;
                delete m_pMomDistrib;

                if ( m_switchIsospinAsymmetry )
                {
                    delete m_nLorentz;
                    delete m_nHoleSF;
                    delete m_nMomDistrib;
                }
            }

            else
            {
                delete m_pLorentz_MF;
                delete m_pHoleSF_MF;
                delete m_pMomDistrib_MF;

                delete m_pLorentz_corr;
                delete m_pHoleSF_corr;
                delete m_pMomDistrib_corr;

                if ( m_switchIsospinAsymmetry )
                {
                    delete m_nLorentz_MF;
                    delete m_nHoleSF_MF;
                    delete m_nMomDistrib_MF;

                    delete m_nLorentz_corr;
                    delete m_nHoleSF_corr;
                    delete m_nMomDistrib_corr;
                }
            }

            if ( m_switchFSI )
            {
                for (auto ptr : m_transparency)
                delete ptr;
                delete m_realOP;
                delete m_foldingF;
            }
        };

        /// proton SF value at (p,E)
        inline double eval_protonSF(const double p, const double removE) const
        {
            return ( removE > m_pHoleSF->get_yStart() and removE < m_pHoleSF->get_yStop() ) ? m_Z*m_pHoleSF->linearI(p, removE) : 0.0;
        }
        /// neutron SF value at (p,E) with energy shift
        inline double eval_neutronSF(const double p, const double removE) const
        {
            const double shiftedRemovalEnergy( removE - m_neutronEnergyShift );

            if ( m_switchIsospinAsymmetry )
            {
                return ( shiftedRemovalEnergy > m_nHoleSF->get_yStart() and shiftedRemovalEnergy < m_nHoleSF->get_yStop() ) ? m_N*m_nHoleSF->linearI(p, shiftedRemovalEnergy) : 0.0;
            }

            return ( shiftedRemovalEnergy > m_pHoleSF->get_yStart() and shiftedRemovalEnergy < m_pHoleSF->get_yStop() ) ? m_N*m_pHoleSF->linearI(p, shiftedRemovalEnergy) : 0.0;
        }
        /// Proton total momentum distribution n_p(p)
        inline double eval_protonMomentumDistribution(const double p) const { return m_pMomDistrib->linearI(p); }
        /// Proton mean-field component of momentum distribution
        inline double eval_protonMomentumDistribution_MF(const double p) const { return m_pMomDistrib_MF->linearI(p);}
        /// Proton correlated component of momentum distribution
        inline double eval_protonMomentumDistribution_corr(const double p) const { return m_pMomDistrib_corr->linearI(p); }
        /// Neutron total momentum distribution n_n(p).
        inline double eval_neutronMomentumDistribution(const double p) const
        {
            if ( m_switchIsospinAsymmetry )
                return m_nMomDistrib->linearI(p);

            return m_pMomDistrib->linearI(p);
        }
        /// Neutron MF momentum distribution.
        inline double eval_neutronMomentumDistribution_MF(const double p) const
        {
            if ( m_switchIsospinAsymmetry )
                return m_nMomDistrib_MF->linearI(p);

            return m_pMomDistrib_MF->linearI(p);
        }
        /// Neutron correlated momentum distribution.
        inline double eval_neutronMomentumDistribution_corr(const double p) const
        {
            if ( m_switchIsospinAsymmetry )
                return m_nMomDistrib_corr->linearI(p);

            return m_pMomDistrib_corr->linearI(p);
        }
        ///generates initial proton momentum according to the total momentum distribution
        double generateProtonMomentum() const
        {
            while ( true )
            {
                const double p( m_pLorentz->generate() );

                const double lor ( m_pLorentz->eval(p) );

                if ( frandom11()*lor < p*p*eval_protonMomentumDistribution(p) )
                    return p;
            }

            return 0.0;
        }
        ///generates initial proton momentum according to the mean-field momentum distribution
        double generateProtonMomentum_MF() const
        {
            while ( true )
            {
                const double p( m_pLorentz_MF->generate() );

                const double lor ( m_pLorentz_MF->eval(p) );

                if ( frandom11()*lor < p*p*eval_protonMomentumDistribution_MF(p) )
                    return p;
            }

            return 0.0;
        }
        ///generates initial proton momentum according to the correlated momentum distribution
        double generateProtonMomentum_corr() const
        {
            while ( true )
            {
                const double p( m_pLorentz_corr->generate() );

                const double lor ( m_pLorentz_corr->eval(p) );

                if ( frandom11()*lor < p*p*eval_protonMomentumDistribution_corr(p) )
                    return p;
            }

            return 0.0;
        }
        ///generates initial neutron momentum according to the total momentum distribution
        double generateNeutronMomentum() const
        {
            while ( m_switchIsospinAsymmetry )
            {
                const double p( m_nLorentz->generate() );

                const double lor ( m_nLorentz->eval(p) );

                if ( frandom11()*lor < p*p*eval_neutronMomentumDistribution(p) )
                    return p;
            }

            return generateProtonMomentum();
        }
        ///generates initial neutron momentum according to the mean-field momentum distribution
        double generateNeutronMomentum_MF() const
        {
            while ( m_switchIsospinAsymmetry )
            {
                const double p( m_nLorentz_MF->generate() );

                const double lor ( m_nLorentz_MF->eval(p) );

                if ( frandom11()*lor < p*p*eval_neutronMomentumDistribution_MF(p) )
                    return p;
            }

            return generateProtonMomentum_MF();
        }
        ///generates initial proton momentum according to the correlated momentum distribution
        double generateNeutronMomentum_corr() const
        {
            while ( m_switchIsospinAsymmetry )
            {
                const double p( m_nLorentz_corr->generate() );

                const double lor ( m_nLorentz_corr->eval(p) );

                if ( frandom11()*lor < p*p*eval_neutronMomentumDistribution_corr(p) )
                    return p;
            }

            return generateProtonMomentum_corr();
        }
        ///generates proton's removal energy according to its distribution calculated from the spectral function
        inline double generateProtonRemovalEnergy(const double p) const { return m_pHoleSF->find_yFromCDF( p, frandom11() ); }
        ///generates proton's removal energy according to its distribution calculated from the mean-field spectral function
        inline double generateProtonRemovalEnergy_MF(const double p) const { return m_pHoleSF_MF->find_yFromCDF( p, frandom11() ); }
        ///generates proton's removal energy according to its distribution calculated from the correlated spectral function
        inline double generateProtonRemovalEnergy_corr(const double p) const { return m_pHoleSF_corr->find_yFromCDF( p, frandom11() ); }
        ///generates neutron's removal energy by shifting the proton's energy
        inline double generateNeutronRemovalEnergy(const double p) const
        {
            if ( m_switchIsospinAsymmetry )
                return m_neutronEnergyShift + m_nHoleSF->find_yFromCDF( p, frandom11() );
            return m_neutronEnergyShift + m_pHoleSF->find_yFromCDF( p, frandom11() );
        }
        ///generates neutron's removal energy according to its distribution calculated from the mean-field spectral function
        inline double generateNeutronRemovalEnergy_MF(const double p) const
        {
            if ( m_switchIsospinAsymmetry )
                return m_neutronEnergyShift + m_nHoleSF_MF->find_yFromCDF( p, frandom11() );
            return m_neutronEnergyShift + m_pHoleSF_MF->find_yFromCDF( p, frandom11() );
        }
        ///generates neutron's removal energy according to its distribution calculated from the correlated spectral function
        inline double generateNeutronRemovalEnergy_corr(const double p) const
        {
            if ( m_switchIsospinAsymmetry )
                return m_neutronEnergyShift + m_nHoleSF_corr->find_yFromCDF( p, frandom11() );
            return m_neutronEnergyShift + m_pHoleSF_corr->find_yFromCDF( p, frandom11() );
        }
        ///calculates the correlated spectral function's integrand
        double eval_correlatedSFIntegrand_Ca40(const double p, const double root, const double numA, const double kTwo, const double cosTheta) const
        {
            ///Ciofi degli Atti & Simula, PRC 53, 1689 (1996)
            const double m_alpha_CM( 0.98 * fm2 );
            const double coef_CM( std::pow(m_alpha_CM/Pi, 3.0/2.0) );

            const double kTwoSq( kTwo*kTwo );
            const double pSq( p*p );

            const double kCMSq( pSq + 2*p*kTwo*cosTheta + kTwoSq );
            const double momDistrib_CM( coef_CM*exp(-m_alpha_CM*kCMSq) );

            const double kRelativeSq( 0.25*(pSq - 2*p*kTwo*cosTheta + kTwoSq) );

            const double coef_CA( 4.4 );
            const double coef_rel( coef_CA*fm3/(4*Pi) );
            const double momDistrib_rel(   coef_rel*( 0.23444*exp(-3.2272*fm2*kRelativeSq) + 0.006989*exp(-0.23308*fm2*kRelativeSq) )   );

            const double coef( (numA - 2)/(numA - 1)*M );

            return coef * momDistrib_rel * momDistrib_CM * kTwoSq / root;
        }
        ///calculates possible values for the additional nucleon momentum. There are two solutions in general.
        double find_pSolutions(const double p, const double kappaSq, const double numA, double& weightRatio, double& kTwoPlus, double& kTwoMinus, const double cosTheta) const
        {
            const double pDiv(   p/(numA - 1.0)   );

            double radicand(   pDiv*pDiv*(cosTheta*cosTheta - 1.0) + kappaSq   );
            // std::cout << "radicand = " << radicand << std::endl;

            ///if the kinematics cannot be solved, return 0, i.e., the value for which a cosine cannot be determined
            if (radicand < 0.0) return 0.0;

            const double root( sqrt(radicand) );

            kTwoPlus =  -pDiv*cosTheta + root;
            kTwoMinus = -pDiv*cosTheta - root;

            const double weightPlus(    eval_correlatedSFIntegrand_Ca40(p, root, numA, kTwoPlus, cosTheta)   );
            const double weightMinus(   (kTwoMinus > 0.0) ? eval_correlatedSFIntegrand_Ca40(p, root, numA, kTwoMinus, cosTheta) : 0.0   );

            const double result( weightPlus + weightMinus );

            weightRatio = weightMinus/result;
            // std::cout << "weightPlus = " << weightPlus << ", weightMinus = " << weightMinus << ", result = " << result << std::endl;

            if (result == 0.0) {
            weightRatio = 0.0;
            return 0.0;
            }

            return result;
        }
        ///generate the norm of the additional nucleon's momentum and the cosine of its angle with respect to the struck nucleon's momentum
        void generateNormAndCosine(const double p, const double kappaSq, const double numA, double& kTwo, double& cosTheta) const
        {
            double kTwoPlus(0.0);
            double kTwoMinus(0.0);
            double weightRatio(0.0);

            const double yMax( find_pSolutions(p, kappaSq, numA, weightRatio, kTwoPlus, kTwoMinus, -1.0) );
            const double yMid( find_pSolutions(p, kappaSq, numA, weightRatio, kTwoPlus, kTwoMinus,  0.0) );
            const double yMin( find_pSolutions(p, kappaSq, numA, weightRatio, kTwoPlus, kTwoMinus,  1.0) );

            const double ySumN( yMax + yMid );
            const double ySumP( yMid + yMin );
            const double yDifferenceN( yMax - yMid );
            const double yDifferenceP( yMid - yMin );

            const double pRatio( ySumP/(ySumP + ySumN) );

            double cosine(0.0), r(0.0);

            while (true)
            {
                r = frandom();

                if( r < pRatio )
                {
                    r /= pRatio;///recycling of the random variable

                    cosine =  ( yMid - sqrt(yMid*yMid - r*ySumP*yDifferenceP) )/yDifferenceP;

                    if( frandom()*(yMid - yDifferenceP*cosine) < find_pSolutions(p, kappaSq, numA, weightRatio, kTwoPlus, kTwoMinus, cosine) )
                    {
                        cosTheta = cosine;
                        kTwo = ( frandom() < weightRatio ) ? kTwoMinus : kTwoPlus;
                        break;
                    }
                }

                else
                {
                    r =  (r - pRatio)/(1.0 - pRatio);///recycling of the random variable

                    cosine = ( yMid - sqrt(yMax*yMax - r*ySumN*yDifferenceN) )/yDifferenceN;

                    if( frandom()*(yMid - yDifferenceN*cosine) < find_pSolutions(p, kappaSq, numA, weightRatio, kTwoPlus, kTwoMinus, cosine) )
                    {
                        cosTheta = cosine;
                        kTwo = ( frandom() < weightRatio ) ? kTwoMinus : kTwoPlus;
                        break;
                    }
                }
            }
        }
        ///generates additional protons's 3-momentum for a struck neutron
        vec generateAdditionalProtonsMomentum(const double p, const double removalE, const vec& p_direction) const
        {
            ///event impossible, return zero momentum
            if ( removalE < m_threshold_E2_n)
            {
             //std::cout << "removalE (" << removalE << ") < m_threshold_E2_n (" << m_threshold_E2_n << ")" << std::endl;
             return vec(0.0, 0.0, 0.0);
            }

            double kNorm(0.0);   ///momentum magnitude
            double cosTheta(0.0);///cosine of the angle between p and k

            const double numA( titaniumA );
            const double kappaSq(   2*M*(numA - 2.0)/(numA - 1.0)*(removalE - m_threshold_E2_n)   );

            //std::cout << "kappaSq = " << kappaSq << std::endl;
            if (kappaSq <= 0.0) {
            std::cout << "kappaSq is non-positive." << std::endl;
            return vec(0.0, 0.0, 0.0);
            }

            generateNormAndCosine(p, kappaSq, numA, kNorm, cosTheta);

            if (kNorm == 0.0) {
            std::cout << "kNorm is zero after generateNormAndCosine." << std::endl;
            return vec(0.0, 0.0, 0.0);
            }

            // Sample azimuthal angle phi uniformly between 0 and 2*Pi
            double phi = 2.0 * Pi * frandom11();

            // Compute sinTheta, ensuring numerical stability
            double sinThetaSq = 1.0 - cosTheta * cosTheta;
            double sinTheta = (sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0;

            // Components in the coordinate system where p_vector is along the z-axis
            double kx_local = kNorm * sinTheta * cos(phi);
            double ky_local = kNorm * sinTheta * sin(phi);
            double kz_local = kNorm * cosTheta;

            double k2 = kx_local*kx_local + ky_local*ky_local + kz_local*kz_local;

            if (k2 == 0.0) {
            std::cout << "k2 is zero after generateNormAndCosine." << std::endl;
            return vec(0.0, 0.0, 0.0);
            }

            // Local direction vector of k in p_vector's coordinate system
            vec k_local(kx_local, ky_local, kz_local);
            //  std::cout << k_local.length() << std::endl;
            vec k_rotated = k_local.fromZto(p_direction);
            if (k_rotated.length() < 1e-6 || !std::isfinite(k_rotated.length())) {
                std::cerr << "[DBG][SF] Zero-mom (after rotation): "
                          << "removalE=" << removalE
                          << " kNorm=" << kNorm
                          << " p_direction=(" << p_direction.x << "," << p_direction.y << "," << p_direction.z << ")\n";
            }

            return k_rotated;
        }
        ///generates additional neutron's 3-momentum for a struck proton
        vec generateAdditionalNeutronsMomentum(const double p, const double removalE, const vec& p_direction) const
        {
            ///event impossible, return zero momentum
            if ( removalE < m_threshold_E2_p)
            {
             //std::cout << "removalE (" << removalE << ") < m_threshold_E2_p (" << m_threshold_E2_p << ")" << std::endl;
             return vec(0.0, 0.0, 0.0);
            }

            double kNorm(0.0);   ///momentum magnitude
            double cosTheta(0.0);///cosine of the angle between p and k

            const double numA( argonA );
            const double kappaSq(   2*M*(numA - 2.0)/(numA - 1.0)*(removalE - m_threshold_E2_p)   );

            //std::cout << "kappaSq = " << kappaSq << std::endl;
            if (kappaSq <= 0.0) {
            std::cout << "kappaSq is non-positive." << std::endl;
            return vec(0.0, 0.0, 0.0);
            }

            generateNormAndCosine(p, kappaSq, numA, kNorm, cosTheta);
            if (kNorm == 0.0) {
            std::cout << "kNorm is zero after generateNormAndCosine." << std::endl;
            return vec(0.0, 0.0, 0.0);
            }

            // Sample azimuthal angle phi uniformly between 0 and 2*Pi
            double phi = 2.0 * Pi * frandom11();

            // Compute sinTheta, ensuring numerical stability
            double sinThetaSq = 1.0 - cosTheta * cosTheta;
            double sinTheta = (sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0;

            // Components in the coordinate system where p_vector is along the z-axis
            double kx_local = kNorm * sinTheta * cos(phi);
            double ky_local = kNorm * sinTheta * sin(phi);
            double kz_local = kNorm * cosTheta;

            double k2 = kx_local*kx_local + ky_local*ky_local + kz_local*kz_local;

            if (k2 == 0.0) {
            std::cout << "k2 is zero after generateNormAndCosine." << std::endl;
            return vec(0.0, 0.0, 0.0);
            }

            // Local direction vector of k in p_vector's coordinate system
            vec k_local(kx_local, ky_local, kz_local);

            vec k_rotated = k_local.fromZto(p_direction);

            return k_rotated;
        }
        /// Evaluate particle SF
        inline double eval_particleSF(const int nucleon_pdg, const double pPrime) const
        {
            const double pF (    ( (not m_switchIsospinAsymmetry) or nucleon_pdg == pdg_proton ) ? m_pFermiMomentum : m_nFermiMomentum );
            return (pPrime < pF) ? 0.0 : 1.0;
        }
        /// Real part of the ntial U as a function of outgoing kinetic energy)
        inline double eval_realOP(const double tPPrime) const { return m_realOP->linearI(tPPrime); }
        /// scaled transparency table from selected transparency table (0=theory, 1=MC)
        inline double eval_sqrtOfTransparency(const double tPPrime, const double scale, int table_idx) const { return sqrt(scale * m_transparency[table_idx]->linearI(tPPrime)); }
        /// transparency table from selected transparency table (0=theory, 1=MC)
        inline double eval_sqrtOfTransparency(const double tPPrime, int table_idx) const { return sqrt(m_transparency[table_idx]->linearI(tPPrime)); }
        /// Folding function
        inline double eval_foldingF(const double differ) const { return m_foldingF->linearI(differ); }

        // Basic getters
        inline TargetNucleus get_target() const { return m_target; }
        inline int get_N() const { return m_N; }
        inline int get_Z() const { return m_Z; }
        inline double get_targetMass() const { return m_targetMass; }
        inline double get_CoulombAvEnergy() const { return m_CoulombAvEnergy; }
        inline double get_CoulombMaxEnergy() const { return m_CoulombMaxEnergy; }

        /// Momentum table coverage [p_min, p_max] for hole SF
        inline double get_pMin() const { return m_pHoleSF->get_xStart(); }
        inline double get_pMax() const { return m_pHoleSF->get_xStop(); }

        /// Removal-energy coverage for protons/neutrons (taking energy shift into account)
        inline double get_eMin() const
        {
            return m_switchIsospinAsymmetry
            ? std::min(m_pHoleSF->get_yStart(), m_nHoleSF->get_yStart() + m_neutronEnergyShift)
            : std::min(m_pHoleSF->get_yStart(), m_pHoleSF->get_yStart() + m_neutronEnergyShift);
        }

        inline double get_protonEMax() const { return m_pHoleSF->get_yStop(); }

        inline double get_eMax() const
        {
            return m_switchIsospinAsymmetry
            ? std::max( m_nHoleSF->get_yStop() + m_neutronEnergyShift, m_pHoleSF->get_yStop() )
            : std::max(m_pHoleSF->get_yStop(), m_pHoleSF->get_yStop() + m_neutronEnergyShift);
        }
        /// Global correlated fraction for protons: ∫corr / (∫mf + ∫corr)
        double get_corrProtonFraction() const
        {
            if ( not m_switchSeparation ) return 0.0;

            const double mf( m_pHoleSF_MF->get_normalization() );
            const double corr( m_pHoleSF_corr->get_normalization() );

            return corr/(mf + corr);
        }
        /// correlated fraction for protons at momentum p
        double eval_corrProtonFraction(const double p) const
        {
            if ( not m_switchSeparation )
                return 0.0;

            if ( p > m_pMomDistrib_MF->get_xStop() )
                return 1.0;

            const double mf( m_pMomDistrib_MF->eval(p) );
            const double corr( m_pMomDistrib_corr->eval(p) );
            //std::cout<<p<<" "<<mf<<" "<<corr<<std::endl;
            return corr/(mf + corr);
        }
        /// Global correlated fraction for neutrons
        inline double get_corrNeutronFraction() const
        {
            if ( not m_switchIsospinAsymmetry ) return get_corrProtonFraction();

            if ( not m_switchSeparation ) return 0.0;

            const double mf( m_nHoleSF_MF->get_normalization() );
            const double corr( m_nHoleSF_corr->get_normalization() );

            return corr/(mf + corr);
        }
        /// correlated fraction for neutrons at momentum p
        double eval_corrNeutronFraction(const double p) const
        {
            if ( not m_switchIsospinAsymmetry )
                return get_corrProtonFraction();

            if ( not m_switchSeparation )
                return 0.0;

            if ( p > m_nMomDistrib_MF->get_xStop() )
                return 1.0;

            const double mf( m_nMomDistrib_MF->eval(p) );
            const double corr( m_nMomDistrib_corr->eval(p) );

            return corr/(mf + corr);
        }

        inline double get_protonFermiMomentum() const { return m_pFermiMomentum; }
        inline double get_neutronFermiMomentum() const
        {
            if ( m_switchIsospinAsymmetry ) return m_nFermiMomentum;
                return m_pFermiMomentum;
        }

        const CInterpolatedData* m_foldingF;

    private:
        const TargetNucleus m_target;     // which nucleus this instance represents
        const bool m_switchSeparation;    // MF+corr split (Ar40 only)
        const bool m_switchFSI;           // whether to load FSI-related tables
        bool m_switchIsospinAsymmetry;    // Ar: treat p and n with different inputs

        int m_N = 0;                      // neutron number
        int m_Z = 0;                      // proton number
        double m_targetMass = 0.0;        // nuclear mass
        double m_neutronEnergyShift = 0.0;// E shift applied to neutron SF
        double m_CoulombAvEnergy = 0.0;   // average Coulomb energy
        double m_CoulombMaxEnergy = 0.0;  // max Coulomb energy
        double m_pFermiMomentum = 0.0;    // proton Fermi momentum
        double m_nFermiMomentum = 0.0;    // neutron Fermi momentum (if asymmetry)

        double m_pNormalization = 0.0;    // converts table units to MeV^-4 (proton)
        double m_nNormalization = 0.0;    // same for neutron
        const CLorentzDistrib* m_pLorentz = nullptr; // envelope for p
        const CLorentzDistrib* m_nLorentz = nullptr; // envelope for n
        const CInterpolatedData2D* m_pHoleSF = nullptr; // S_p(p,E)
        const CInterpolatedData2D* m_nHoleSF = nullptr; // S_n(p,E)
        const CInterpolatedData* m_pMomDistrib = nullptr; // n_p(p)
        const CInterpolatedData* m_nMomDistrib = nullptr; // n_n(p)

        double m_pNormalization_MF;
        double m_nNormalization_MF;
        const CLorentzDistrib* m_pLorentz_MF;
        const CLorentzDistrib* m_nLorentz_MF;
        const CInterpolatedData2D* m_pHoleSF_MF;
        const CInterpolatedData2D* m_nHoleSF_MF;
        const CInterpolatedData* m_pMomDistrib_MF;
        const CInterpolatedData* m_nMomDistrib_MF;

        double m_pNormalization_corr;
        double m_nNormalization_corr;
        const CLorentzDistrib* m_pLorentz_corr;
        const CLorentzDistrib* m_nLorentz_corr;
        const CInterpolatedData2D* m_pHoleSF_corr;
        const CInterpolatedData2D* m_nHoleSF_corr;
        const CInterpolatedData* m_pMomDistrib_corr;
        const CInterpolatedData* m_nMomDistrib_corr;
        double m_threshold_E2_p;
        double m_threshold_E2_n;

        std::array<const CInterpolatedData*, 2> m_transparency;
        const CInterpolatedData* m_realOP;
};

#endif // C_SPECTRAL_FUNCTION_H
