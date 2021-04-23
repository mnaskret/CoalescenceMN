#ifndef _CoalescerMN_CoalescerMN_h_
#define _CoalescerMN_CoalescerMN_h_

/**
 * \file
 * \author Michal Naskret
 * \date 10 Feb 2021
 */

#include <fwk/VModule.h>
#include <boost/utility.hpp>
#include <utl/Vector.h>
#include <evt/Event.h>
#include <utl/ShineUnits.h>
#include <utl/DatabasePDG.h>

#include <string>
#include <TLorentzVector.h>
#include <TVector3.h>

namespace CoalescerMN {

  /**
  \class CoalescerMN

  \brief A class used to coalesce protons and neutrons into deuterons or antiprotons and antineutrons into antideuterons.

  Now here a longer description in doxygen. You may use HTML.


  \author Michal Naskret
  \date 10 Feb 2021
  \version 
  \ingroup TODO: Put a group here.
  */

  class CoalescerMN : 
    public boost::noncopyable,
    public fwk::VModule {
  public:
    CoalescerMN() {}
    ~CoalescerMN() {}
    fwk::VModule::EResultFlag Init();
    fwk::VModule::EResultFlag Process(evt::Event& event, const utl::AttributeMap& attr);
    fwk::VModule::EResultFlag Finish();
  private:

    typedef std::list< evt::sim::VertexTrack >::const_iterator VtxTrackIterator;
    typedef std::vector<VtxTrackIterator> nucleonVector;
    typedef std::vector<TLorentzVector> combinedFourMomenta;

    utl::DatabasePDG& dPDG = utl::DatabasePDG::GetInstance();

    int listParticles(nucleonVector&, std::string);
    int checkIfCloseInMomentumSpace(TLorentzVector, TLorentzVector, int);
    int createParticle(evt::SimEvent&, combinedFourMomenta, int);
    int createSimultanousHelium4(evt::SimEvent&, nucleonVector&, nucleonVector&, int);
    int createSimultanousHelium3(evt::SimEvent&, nucleonVector&, nucleonVector&, int);
    int createSimultanousTriton(evt::SimEvent&, nucleonVector&, nucleonVector&, int);
    int createSimultanousDeuteron(evt::SimEvent&, nucleonVector&, nucleonVector&, int);

    utl::Vector fBeamMomentum;
    double fP0;
    double fAntiP0;
    int fProtonInitialNumber;
    int fNeutronInitialNumber;
    int fAntiProtonInitialNumber;
    int fAntiNeutronInitialNumber;

    int fProtonNumber;
    int fNeutronNumber;
    int fDeuteronNumber;
    int fHe4Number;
    int fHe3Number;
    int fTritonNumber;
    int fAntiProtonNumber;
    int fAntiNeutronNumber;
    int fAntiDeuteronNumber;
    int fAntiHe4Number;
    int fAntiHe3Number;
    int fAntiTritonNumber;
    TVector3 fBoostVectorToCOM;

    REGISTER_MODULE("CoalescerMN", CoalescerMN, "$Id$");
  };
}

#endif // _CoalescerMN_CoalescerMN_h_
