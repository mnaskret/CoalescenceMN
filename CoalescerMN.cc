#include "CoalescerMN.h"

#include <fwk/CentralConfig.h>

#include <evt/Event.h>
#include <evt/IndexedObjectLinker.h>

#include <utl/Branch.h>
#include <utl/ErrorLogger.h>
#include <utl/ShineUnits.h>
#include <utl/PhysicalConst.h>
#include <utl/DatabasePDG.h>
#include <utl/PDGParticleIds.h>

#include <string>
#include <sstream>
#include <TLorentzVector.h>
#include <TVector3.h>

using namespace fwk;
using namespace utl;
using namespace std;

namespace CoalescerMN {

  // Define the methods preferibly in the same order as listed in the header file.

  VModule::EResultFlag 
  CoalescerMN::Init()
  {
    // Initialize your module here. This method
    // is called once at the beginning of the run.
    // The eSuccess flag indicates the method ended
    // successfully.  For other possible return types,
    // see the VModule documentation.

    CentralConfig& cc = CentralConfig::GetInstance();

    Branch topBranch = cc.GetTopBranch("CoalescerMN");
    fP0 = topBranch.GetChild("p0").Get<double>()*utl::MeV;
    fAntiP0 = topBranch.GetChild("antiP0").Get<double>()*utl::MeV;
    fBeamMomentum.Set(0, 0, topBranch.GetChild("beamMomentum").Get<double>()*utl::GeV);
    int targetAtomicNumber = topBranch.GetChild("targetAtomicNumber").Get<int>();
    int targetMassNumber = topBranch.GetChild("targetMassNumber").Get<int>();
    int beamAtomicNumber = topBranch.GetChild("beamAtomicNumber").Get<int>();
    int beamMassNumber = topBranch.GetChild("beamMassNumber").Get<int>();
//    InitVerbosity(topBranch);

    std::ostringstream msg;
    msg << "Running the CoalescerMN module for beam momentum "<<fBeamMomentum.GetZ()<<" with p0="<<fP0<<"\n";
    INFO(msg);

    fProtonInitialNumber = 0;
    fNeutronInitialNumber = 0;
    fProtonNumber = 0;
    fNeutronNumber = 0;
    fDeuteronNumber = 0;
    fHe4Number = 0;
    fHe3Number = 0;
    fTritonNumber = 0;
    fAntiDeuteronNumber = 0;

    TLorentzVector targetFourMomentum;
    TLorentzVector beamFourMomentum;
    TLorentzVector COMFourMomentum;

    // Calculate boost vector to COM system
    double targetNucleonMass = targetAtomicNumber*utl::kProtonMass + (targetMassNumber-targetAtomicNumber)*utl::kNeutronMass;
    double beamNucleonMass = beamAtomicNumber*utl::kProtonMass + (beamMassNumber-beamAtomicNumber)*utl::kNeutronMass;
    targetFourMomentum.SetPxPyPzE(0, 0, 0, targetNucleonMass);
    beamFourMomentum.SetPxPyPzE(0, 0, fBeamMomentum.GetZ(), sqrt(fBeamMomentum.GetZ()*fBeamMomentum.GetZ()+beamNucleonMass*beamNucleonMass));

    COMFourMomentum = targetFourMomentum + beamFourMomentum;
    fBoostVectorToCOM = COMFourMomentum.BoostVector();

    msg.clear();
    msg.str("");
    msg<<"COMFourMomentum and boost: \n";
    msg<<COMFourMomentum[0]<<"\t"<<COMFourMomentum[1]<<"\t"<<COMFourMomentum[2]<<"\t"<<COMFourMomentum[3]<<"\n";
    msg<<fBoostVectorToCOM[0]<<"\t"<<fBoostVectorToCOM[1]<<"\t"<<fBoostVectorToCOM[2]<<"\n";
    INFO(msg);

    return eSuccess;
  }

  VModule::EResultFlag
  CoalescerMN::Process(evt::Event& event, const AttributeMap& /*attr*/)
  {

    std::ostringstream msg;
    msg << "CoalescerMN processing event\n";
    INFO(msg);

    //Create vectors for p, n , pbar and nbar
    nucleonVector neutrons;
    nucleonVector antineutrons;
    nucleonVector protons;
    nucleonVector antiprotons;

    evt::SimEvent& simEvent = event.GetSimEvent();

    //Iterate through Vtx tracks to find p, n, pbar and nbar
    std::set< evt::Index<evt::sim::VertexTrack> > indicesToSplit;
    for (VtxTrackIterator tIter  = simEvent.Begin<evt::sim::VertexTrack>(),
			  endIter  = simEvent.End<evt::sim::VertexTrack>();
			  tIter != endIter;
			  ++tIter)
    {
      const evt::sim::VertexTrack& vtxTrack = *tIter;

      // Check if particle is final generator
      if (vtxTrack.GetType() != evt::sim::VertexTrackConst::eGeneratorFinal)
        continue;

      const int pid = vtxTrack.GetParticleId();
      // Check if proton or neutron
      if (pid != 2212 && pid != 2112 && pid != -2212 && pid != -2112)
        continue;

      // If proton add to proton list
      if (pid == 2212)
      {
        protons.push_back(tIter);
      }

      // If neutron add to neutron list
      if (pid == 2112)
      {
        neutrons.push_back(tIter);
      }

      // If antiproton add to antiproton list
      if (pid == -2212)
      {
        antiprotons.push_back(tIter);
      }

      // If antineutron add to antineutron list
      if (pid == -2112)
      {
        antineutrons.push_back(tIter);
      }
    }

    // Print initial number of particles
    listParticles(protons, "protons");
    listParticles(neutrons, "neutrons");
    listParticles(antiprotons, "antiprotons");
    listParticles(antineutrons, "anineutrons");

    // Save initial number of particles
    fProtonInitialNumber += protons.size();
    fNeutronInitialNumber += neutrons.size();
    fAntiProtonInitialNumber += antiprotons.size();
    fAntiNeutronInitialNumber += antineutrons.size();

    msg.clear();
    msg.str("");
    msg << "Checking for He4 and anti-He4 simultanous\n";
    INFO(msg);

    // Create He4 particles first, then antiparticles
    createSimultanousHelium4(simEvent, protons, neutrons, 0);
    createSimultanousHelium4(simEvent, protons, neutrons, 1);

    msg.clear();
    msg.str("");
    msg << "Checking for He3 and anti-He3 simultanous\n";
    INFO(msg);

    // Create He3 particles first, then antiparticles
    createSimultanousHelium3(simEvent, protons, neutrons, 0);
    createSimultanousHelium3(simEvent, protons, neutrons, 1);

    msg.clear();
    msg.str("");
    msg << "Checking for triton and antitriton simultanous\n";
    INFO(msg);

    // Create tritons first, then antitritons
    createSimultanousTriton(simEvent, protons, neutrons, 0);
    createSimultanousTriton(simEvent, protons, neutrons, 1);

    msg.clear();
    msg.str("");
    msg << "Checking for deuterons and antideuterons simultanous\n";
    INFO(msg);

    // Create deuterons first, then antideuterons
    createSimultanousDeuteron(simEvent, protons, neutrons, 0);
    createSimultanousDeuteron(simEvent, protons, neutrons, 1);

    fProtonNumber += protons.size();
    fNeutronNumber += neutrons.size();

    return eSuccess;
  }

  VModule::EResultFlag
  CoalescerMN::Finish() 
  {
    std::ostringstream msg;
    msg << "There are "<<fProtonInitialNumber
        <<" protons and "<<fNeutronInitialNumber
        <<" neutrons in the beginning.\n";
    msg << "CoalescerMN found "<<fDeuteronNumber
        <<" deuterons, "<<fHe3Number
        <<" He3, "<<fTritonNumber
        <<" tritons and "<<fHe4Number<<" He4\n";
    msg << "There are remaining "
        <<fProtonNumber<<" protons, "
        <<fNeutronNumber<<" neutrons.\n";
    INFO(msg);

    return eSuccess;
  }

  int CoalescerMN::listParticles(nucleonVector &particlesVector, string particlesName)
  {
    std::ostringstream msg;
    msg<<"List of "<<particlesName<<" in the event: \n";
    for(nucleonVector::iterator particleIt = std::begin(particlesVector);
        particleIt != std::end(particlesVector);
        ++particleIt)
    {
      const evt::sim::VertexTrack& particleVtxTrack = *(*particleIt);
      utl::Vector particleMomentum = particleVtxTrack.GetMomentum();
      msg<<particleVtxTrack.GetIndex()
         <<", momentum "<<particleMomentum.GetX()
         <<"\t"<<particleMomentum.GetY()
         <<"\t"<<particleMomentum.GetZ()<<"\n";
    }
    INFO(msg);
  }

  int CoalescerMN::createParticle(evt::SimEvent& simEvent, combinedFourMomenta combinedMomenta, int pdgId)
  {
    //Sum component particles momenta boosted to COM
    TLorentzVector particleFourMomentum;
    for(std::size_t i = 0; i < combinedMomenta.size(); ++i) {
      combinedMomenta[i].Boost(-fBoostVectorToCOM);
      particleFourMomentum += combinedMomenta[i];
    }

    //Grab correct particle mass
    double particleMass = dPDG.GetParticle(pdgId).Mass()*utl::GeV;

    //Calculate correct energy and boost back to LAB frame
    particleFourMomentum.SetE(sqrt(pow(particleFourMomentum.Px(),2)+pow(particleFourMomentum.Py(),2)+pow(particleFourMomentum.Pz(),2)+pow(particleMass,2)));
    particleFourMomentum.Boost(fBoostVectorToCOM);

    //Add a new Vtx track to Sim event and set its type
    utl::Vector momentum(particleFourMomentum.Px(), particleFourMomentum.Py(), particleFourMomentum.Pz());
    evt::sim::VertexTrack& particleVtxTrack = simEvent.Make<evt::sim::VertexTrack>();
    particleVtxTrack.SetParticleId(pdgId);
    particleVtxTrack.SetMomentum(momentum);
    particleVtxTrack.SetType(evt::sim::VertexTrackConst::eGeneratorCoalesced);

    //Link created Vtx track to main vertex
    evt::sim::Vertex& mainVtx = simEvent.Get(simEvent.GetMainVertex().GetIndex());
    evt::IndexedObjectLinker::LinkDaughterToParent(particleVtxTrack, mainVtx);

    return 1;
  }

  int CoalescerMN::checkIfCloseInMomentumSpace(TLorentzVector firstFourMomentum, TLorentzVector secondFourMomentum, int isAntiParticle)
  {
    //Read correct p0 value
    double p0;
    if (isAntiParticle) p0 = fAntiP0;
    else p0 = fP0;

    //Create and initialize Lorentz Vectors
    TLorentzVector commonPNFourMomentum;

    //Boost to COM frame
    firstFourMomentum.Boost(-fBoostVectorToCOM);
    secondFourMomentum.Boost(-fBoostVectorToCOM);

    //Find and boost to common coalesced particles frame
    commonPNFourMomentum = firstFourMomentum + secondFourMomentum;
    TVector3 boostVectorToCommonPN = commonPNFourMomentum.BoostVector();

    firstFourMomentum.Boost(-boostVectorToCommonPN);
    secondFourMomentum.Boost(-boostVectorToCommonPN);

    //Create vectors of momenta of particles boosted to common frame
    utl::Vector firstMomentumCommon = utl::Vector(firstFourMomentum[0], firstFourMomentum[1], firstFourMomentum[2]);
    utl::Vector secondMomentumCommon = utl::Vector(secondFourMomentum[0], secondFourMomentum[1], secondFourMomentum[2]);

    //Check coalescence condition
    if((firstMomentumCommon-secondMomentumCommon).GetMag() < 2*fP0) return 1;

    return 0;
  }

  int CoalescerMN::createSimultanousHelium4(evt::SimEvent& simEvent, nucleonVector& protons, nucleonVector& neutrons, int isAntiParticle)
  {
    //Iterate through vector of protons
    for(nucleonVector::iterator protonIt1 = std::begin(protons); protonIt1 != std::end(protons); ++protonIt1)
    {
      int createdHelium4 = 0;

      //Grab first proton vertex track, momentum and energy
      const evt::sim::VertexTrack& protonVtxTrack1 = *(*protonIt1);
      utl::Vector protonMomentum1 = protonVtxTrack1.GetMomentum();
      double protonEnergy1 = sqrt(utl::kProtonMass*utl::kProtonMass + protonMomentum1.GetMag2());

      //Iterate through vector of neutrons
      for(nucleonVector::iterator neutronIt1 = std::begin(neutrons); neutronIt1 != std::end(neutrons); ++neutronIt1)
      {
        if(createdHelium4) break;

        //Grab first neutron vertex track, momentum and energy
        const evt::sim::VertexTrack& neutronVtxTrack1 = *(*neutronIt1);
        utl::Vector neutronMomentum1 = neutronVtxTrack1.GetMomentum();
        double neutronEnergy1 = sqrt(utl::kNeutronMass*utl::kNeutronMass + neutronMomentum1.GetMag2());

        //Iterate through vector of protons starting from the first proton after protonIt1
        for(nucleonVector::iterator protonIt2 = protonIt1+1; protonIt2 != std::end(protons); ++protonIt2)
        {
          if(createdHelium4) break;

          //Grab second proton vertex track, momentum and energy
          const evt::sim::VertexTrack& protonVtxTrack2 = *(*protonIt2);
          utl::Vector protonMomentum2 = protonVtxTrack2.GetMomentum();
          double protonEnergy2 = sqrt(utl::kProtonMass*utl::kProtonMass + protonMomentum2.GetMag2());

          //Iterate through vector of neutrons starting from the second neutron after neutronIt1
          for(nucleonVector::iterator neutronIt2 = neutronIt1+1; neutronIt2 != std::end(neutrons); ++neutronIt2)
          {

            //Grab second neutron vertex track, momentum and energy
            const evt::sim::VertexTrack& neutronVtxTrack2 = *(*neutronIt2);
            utl::Vector neutronMomentum2 = neutronVtxTrack2.GetMomentum();
            double neutronEnergy2 = sqrt(utl::kNeutronMass*utl::kNeutronMass + neutronMomentum2.GetMag2());

            // Create particles four momenta
            TLorentzVector protonFourMomentum1;
            TLorentzVector protonFourMomentum2;
            TLorentzVector neutronFourMomentum1;
            TLorentzVector neutronFourMomentum2;

            protonFourMomentum1.SetPxPyPzE(protonMomentum1.GetX(), protonMomentum1.GetY(), protonMomentum1.GetZ(), protonEnergy1);
            protonFourMomentum2.SetPxPyPzE(protonMomentum2.GetX(), protonMomentum2.GetY(), protonMomentum2.GetZ(), protonEnergy2);
            neutronFourMomentum1.SetPxPyPzE(neutronMomentum1.GetX(), neutronMomentum1.GetY(), neutronMomentum1.GetZ(), neutronEnergy1);
            neutronFourMomentum2.SetPxPyPzE(neutronMomentum2.GetX(), neutronMomentum2.GetY(), neutronMomentum2.GetZ(), neutronEnergy2);

            //Check if all 4 particles are close in the momentum space
            if(checkIfCloseInMomentumSpace(protonFourMomentum1, neutronFourMomentum1, isAntiParticle) && 
               checkIfCloseInMomentumSpace(protonFourMomentum1, neutronFourMomentum2, isAntiParticle) && 
               checkIfCloseInMomentumSpace(protonFourMomentum2, neutronFourMomentum1, isAntiParticle) && 
               checkIfCloseInMomentumSpace(protonFourMomentum2, neutronFourMomentum2, isAntiParticle) &&
               checkIfCloseInMomentumSpace(protonFourMomentum1, protonFourMomentum2, isAntiParticle) && 
               checkIfCloseInMomentumSpace(neutronFourMomentum1, neutronFourMomentum2, isAntiParticle))
            {
              createdHelium4 = 1;
              fHe4Number++;

              //Remove vtx tracks from event
              simEvent.Erase<evt::sim::VertexTrack>(*protonIt1);
              simEvent.Erase<evt::sim::VertexTrack>(*protonIt2);
              simEvent.Erase<evt::sim::VertexTrack>(*neutronIt1);
              simEvent.Erase<evt::sim::VertexTrack>(*neutronIt2);

              //Remove particles from vectors
              protons.erase(protonIt1--);
              protonIt2--;
              protons.erase(protonIt2--);

              neutrons.erase(neutronIt1--);
              neutronIt2--;
              neutrons.erase(neutronIt2--);

              //Save momenta of particles and combine into a vector
              combinedFourMomenta momentumHe4Combined;

              momentumHe4Combined.push_back(protonFourMomentum1);
              momentumHe4Combined.push_back(protonFourMomentum2);
              momentumHe4Combined.push_back(neutronFourMomentum1);
              momentumHe4Combined.push_back(neutronFourMomentum2);

              //Create He4
              if (isAntiParticle) createParticle(simEvent, momentumHe4Combined, -1000020040);
              else createParticle(simEvent, momentumHe4Combined, 1000020040);

              break;
            }
          }
        }
      }
    }
    return 1;
  }

  int CoalescerMN::createSimultanousHelium3(evt::SimEvent& simEvent, nucleonVector& protons, nucleonVector& neutrons, int isAntiParticle)
  {
    //Iterate through vector of protons
    for(nucleonVector::iterator protonIt1 = std::begin(protons); protonIt1 != std::end(protons); ++protonIt1)
    {
      int createdHelium3 = 0;

      //Grab first proton vertex track, momentum and energy
      const evt::sim::VertexTrack& protonVtxTrack1 = *(*protonIt1);
      utl::Vector protonMomentum1 = protonVtxTrack1.GetMomentum();
      double protonEnergy1 = sqrt(utl::kProtonMass*utl::kProtonMass + protonMomentum1.GetMag2());

      //Iterate through vector of neutrons
      for(nucleonVector::iterator neutronIt1 = std::begin(neutrons); neutronIt1 != std::end(neutrons); ++neutronIt1)
      {
        if(createdHelium3) break;

        //Grab first neutron vertex track, momentum and energy
        const evt::sim::VertexTrack& neutronVtxTrack1 = *(*neutronIt1);
        utl::Vector neutronMomentum1 = neutronVtxTrack1.GetMomentum();
        double neutronEnergy1 = sqrt(utl::kNeutronMass*utl::kNeutronMass + neutronMomentum1.GetMag2());

        //Iterate through vector of protons starting from the first proton after protonIt1
        for(nucleonVector::iterator protonIt2 = protonIt1+1; protonIt2 != std::end(protons); ++protonIt2)
        {

          //Grab second neutron vertex track, momentum and energy
          const evt::sim::VertexTrack& protonVtxTrack2 = *(*protonIt2);
          utl::Vector protonMomentum2 = protonVtxTrack2.GetMomentum();
          double protonEnergy2 = sqrt(utl::kProtonMass*utl::kProtonMass + protonMomentum2.GetMag2());

          TLorentzVector protonFourMomentum1;
          TLorentzVector protonFourMomentum2;
          TLorentzVector neutronFourMomentum1;
          protonFourMomentum1.SetPxPyPzE(protonMomentum1.GetX(), protonMomentum1.GetY(), protonMomentum1.GetZ(), protonEnergy1);
          protonFourMomentum2.SetPxPyPzE(protonMomentum2.GetX(), protonMomentum2.GetY(), protonMomentum2.GetZ(), protonEnergy2);
          neutronFourMomentum1.SetPxPyPzE(neutronMomentum1.GetX(), neutronMomentum1.GetY(), neutronMomentum1.GetZ(), neutronEnergy1);

          //Check if all 3 particles are close in the momentum space
          if(checkIfCloseInMomentumSpace(protonFourMomentum1, neutronFourMomentum1, isAntiParticle) && 
             checkIfCloseInMomentumSpace(protonFourMomentum2, neutronFourMomentum1, isAntiParticle) &&
             checkIfCloseInMomentumSpace(protonFourMomentum1, protonFourMomentum2, isAntiParticle))
          {
            createdHelium3 = 1;
            fHe3Number++;

            //Remove vtx tracks from event
            simEvent.Erase<evt::sim::VertexTrack>(*protonIt1);
            simEvent.Erase<evt::sim::VertexTrack>(*protonIt2);
            simEvent.Erase<evt::sim::VertexTrack>(*neutronIt1);

            //Remove particles from vectors
            protons.erase(protonIt1--);
            protonIt2--;
            protons.erase(protonIt2--);

            neutrons.erase(neutronIt1--);

            //Save momenta of particles and combine into a vector
            combinedFourMomenta momentumHe3Combined;

            momentumHe3Combined.push_back(protonFourMomentum1);
            momentumHe3Combined.push_back(protonFourMomentum2);
            momentumHe3Combined.push_back(neutronFourMomentum1);

            //Create He3
            if (isAntiParticle) createParticle(simEvent, momentumHe3Combined, -1000020030);
            else createParticle(simEvent, momentumHe3Combined, 1000020030);

            break;
          }
        }
      }
    }
    return 1;
  }

  int CoalescerMN::createSimultanousTriton(evt::SimEvent& simEvent, nucleonVector& protons, nucleonVector& neutrons, int isAntiParticle)
  {
    //Iterate through vector of protons
    for(nucleonVector::iterator protonIt1 = std::begin(protons); protonIt1 != std::end(protons); ++protonIt1)
    {
      int createdTriton = 0;

      //Grab first proton vertex track, momentum and energy
      const evt::sim::VertexTrack& protonVtxTrack1 = *(*protonIt1);
      utl::Vector protonMomentum1 = protonVtxTrack1.GetMomentum();
      double protonEnergy1 = sqrt(utl::kProtonMass*utl::kProtonMass + protonMomentum1.GetMag2());

      //Iterate through vector of neutrons
      for(nucleonVector::iterator neutronIt1 = std::begin(neutrons); neutronIt1 != std::end(neutrons); ++neutronIt1)
      {
        if(createdTriton) break;

        //Grab first neutron vertex track, momentum and energy
        const evt::sim::VertexTrack& neutronVtxTrack1 = *(*neutronIt1);
        utl::Vector neutronMomentum1 = neutronVtxTrack1.GetMomentum();
        double neutronEnergy1 = sqrt(utl::kNeutronMass*utl::kNeutronMass + neutronMomentum1.GetMag2());

        //Iterate through vector of neutrons starting from the first neutron after neutronIt1
        for(nucleonVector::iterator neutronIt2 = neutronIt1+1; neutronIt2 != std::end(neutrons); ++neutronIt2)
        {

          //Grab second neutron vertex track, momentum and energy
          const evt::sim::VertexTrack& neutronVtxTrack2 = *(*neutronIt2);
          utl::Vector neutronMomentum2 = neutronVtxTrack2.GetMomentum();
          double neutronEnergy2 = sqrt(utl::kNeutronMass*utl::kNeutronMass + neutronMomentum2.GetMag2());

          TLorentzVector protonFourMomentum1;
          TLorentzVector neutronFourMomentum1;
          TLorentzVector neutronFourMomentum2;
          protonFourMomentum1.SetPxPyPzE(protonMomentum1.GetX(), protonMomentum1.GetY(), protonMomentum1.GetZ(), protonEnergy1);
          neutronFourMomentum1.SetPxPyPzE(neutronMomentum1.GetX(), neutronMomentum1.GetY(), neutronMomentum1.GetZ(), neutronEnergy1);
          neutronFourMomentum2.SetPxPyPzE(neutronMomentum2.GetX(), neutronMomentum2.GetY(), neutronMomentum2.GetZ(), neutronEnergy2);

          //Check if all 3 particles are close in the momentum space
          if(checkIfCloseInMomentumSpace(neutronFourMomentum1, protonFourMomentum1, isAntiParticle) && 
             checkIfCloseInMomentumSpace(neutronFourMomentum2, protonFourMomentum1, isAntiParticle) &&
             checkIfCloseInMomentumSpace(neutronFourMomentum1, neutronFourMomentum2, isAntiParticle))
          {
            createdTriton = 1;
            fTritonNumber++;

            //Remove vtx tracks from event
            simEvent.Erase<evt::sim::VertexTrack>(*protonIt1);
            simEvent.Erase<evt::sim::VertexTrack>(*neutronIt1);
            simEvent.Erase<evt::sim::VertexTrack>(*neutronIt2);

            //Remove particles from vectors
            neutrons.erase(neutronIt1--);
            neutronIt2--;
            neutrons.erase(neutronIt2--);

            protons.erase(protonIt1--);

            //Save momenta of particles and combine into a vector
            combinedFourMomenta momentumTritonCombined;

            momentumTritonCombined.push_back(protonFourMomentum1);
            momentumTritonCombined.push_back(neutronFourMomentum1);
            momentumTritonCombined.push_back(neutronFourMomentum2);

            //Create triton
            if (isAntiParticle) createParticle(simEvent, momentumTritonCombined, -1000010030);
            else createParticle(simEvent, momentumTritonCombined, 1000010030);

            break;
          }
        }
      }
    }
    return 1;
  }

  int CoalescerMN::createSimultanousDeuteron(evt::SimEvent& simEvent, nucleonVector& protons, nucleonVector& neutrons, int isAntiParticle)
  {
    //Iterate through vector of protons
    for(nucleonVector::iterator protonIt1 = std::begin(protons); protonIt1 != std::end(protons); ++protonIt1)
    {
      int createdDeuteron = 0;

      //Grab first proton vertex track, momentum and energy
      const evt::sim::VertexTrack& protonVtxTrack1 = *(*protonIt1);
      utl::Vector protonMomentum1 = protonVtxTrack1.GetMomentum();
      double protonEnergy1 = sqrt(utl::kProtonMass*utl::kProtonMass + protonMomentum1.GetMag2());

      //Iterate through vector of protons
      for(nucleonVector::iterator neutronIt1 = std::begin(neutrons); neutronIt1 != std::end(neutrons); ++neutronIt1)
      {

        //Grab first neutron vertex track, momentum and energy
        const evt::sim::VertexTrack& neutronVtxTrack1 = *(*neutronIt1);
        utl::Vector neutronMomentum1 = neutronVtxTrack1.GetMomentum();
        double neutronEnergy1 = sqrt(utl::kNeutronMass*utl::kNeutronMass + neutronMomentum1.GetMag2());

        TLorentzVector protonFourMomentum1;
        TLorentzVector neutronFourMomentum1;
        protonFourMomentum1.SetPxPyPzE(protonMomentum1.GetX(), protonMomentum1.GetY(), protonMomentum1.GetZ(), protonEnergy1);
        neutronFourMomentum1.SetPxPyPzE(neutronMomentum1.GetX(), neutronMomentum1.GetY(), neutronMomentum1.GetZ(), neutronEnergy1);

        //Check if both particles are close in the momentum space
        if(checkIfCloseInMomentumSpace(protonFourMomentum1, neutronFourMomentum1, isAntiParticle))
        {
          createdDeuteron = 1;
          fDeuteronNumber++;

          //Remove vtx tracks from event
          simEvent.Erase<evt::sim::VertexTrack>(*protonIt1);
          simEvent.Erase<evt::sim::VertexTrack>(*neutronIt1);

          //Remove particles from vectors
          neutrons.erase(neutronIt1--);
          protons.erase(protonIt1--);

          //Save momenta of particles and combine into a vector
          combinedFourMomenta momentumDeuteronCombined;

          momentumDeuteronCombined.push_back(protonFourMomentum1);
          momentumDeuteronCombined.push_back(neutronFourMomentum1);

          //Create deuteron
          if (isAntiParticle) createParticle(simEvent, momentumDeuteronCombined, -1000010020);
          else createParticle(simEvent, momentumDeuteronCombined, 1000010020);

          break;
        }
      }
    }
    return 1;
  }
}
