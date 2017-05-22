// -*- C++ -*-
//
// Package:    QCDAnalyzer
// Class:      QCDAnalyzer
// 
/**\class QCDAnalyzer QCDAnalyzer.cc Appeltel/QCDAnalyzer/src/QCDAnalyzer.cc

 Description: Gen-Level analyzer for QCD processes (i.e. jets, charged particles, parton flavor)

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Eric Appelt
//         Created:  Tue Feb 25 14:25:23 CST 2014
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"
#include "HepPDT/ParticleID.hh"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"


#include <TH1.h>
#include <TH2.h>


//
// class declaration
//

class QCDAnalyzer : public edm::EDAnalyzer {
   public:
      explicit QCDAnalyzer(const edm::ParameterSet&);
      ~QCDAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      bool isDSEvent( const edm::Event&, const edm::EventSetup& );
      bool isInEtaRange( const reco::Candidate&, double, double ); 
      bool isInFlavor( const reco::MatchedPartons& ); 
      bool isInSpecies( const reco::GenParticle& ); 

      // ----------member data ---------------------------

      edm::InputTag genJetSrc_;
      edm::InputTag genParticleSrc_;
      bool doFlavor_;
      bool doSpecies_;
      bool onlyDS_;
      edm::InputTag flavorSrc_;
      std::vector<int> flavorId_;
      std::vector<int> speciesId_;
      bool useRapidity_;
      double jetEtaMin_;
      double jetEtaMax_;
      double hEtaMin_;
      double hEtaMax_;
      double jetRadius_;
      double pthatMin_;
      double pthatMax_;
      std::vector<double> jetPtBins_;
      std::vector<double> hPtBins_;
      std::vector<double> qScalePtBins_;
      std::vector<double> etaBins_;

      std::map<std::string,TH1F*> hist_;
      std::map<std::string,TH2F*> hist2D_;
};


//
// constructors and destructor
//

QCDAnalyzer::QCDAnalyzer(const edm::ParameterSet& iConfig):
genJetSrc_(iConfig.getParameter<edm::InputTag>("genJetSrc")),
genParticleSrc_(iConfig.getParameter<edm::InputTag>("genParticleSrc")),
doFlavor_(iConfig.getParameter<bool>("doFlavor")),
doSpecies_(iConfig.getParameter<bool>("doSpecies")),
onlyDS_(iConfig.getParameter<bool>("onlyDS")),
flavorSrc_(iConfig.getParameter<edm::InputTag>("flavorSrc")),
flavorId_(iConfig.getParameter<std::vector<int> >("flavorId")),
speciesId_(iConfig.getParameter<std::vector<int> >("speciesId")),
useRapidity_(iConfig.getParameter<bool>("useRapidity")),
jetEtaMin_(iConfig.getParameter<double>("jetEtaMin")),
jetEtaMax_(iConfig.getParameter<double>("jetEtaMax")),
hEtaMin_(iConfig.getParameter<double>("hEtaMin")),
hEtaMax_(iConfig.getParameter<double>("hEtaMax")),
jetRadius_(iConfig.getParameter<double>("jetRadius")),
pthatMin_(iConfig.getParameter<double>("pthatMin")),
pthatMax_(iConfig.getParameter<double>("pthatMax")),
jetPtBins_(iConfig.getParameter<std::vector<double> >("jetPtBins")),
hPtBins_(iConfig.getParameter<std::vector<double> >("hPtBins")),
qScalePtBins_(iConfig.getParameter<std::vector<double> >("qScalePtBins")),
etaBins_(iConfig.getParameter<std::vector<double> >("etaBins"))
{

   edm::Service<TFileService> fs;

   // Spectra histograms for eta (or y) within specified range

   // jet spectrum
   hist_["jetspectrum"] = fs->make<TH1F>("jetspectrum",";p_{T};counts",
                           jetPtBins_.size()-1, &jetPtBins_[0]);
   // hadron spectrum
   hist_["hspectrum"] = fs->make<TH1F>("hspectrum",";p_{T};counts",
                           hPtBins_.size()-1, &hPtBins_[0]);
   // charged hadron spectrum
   hist_["chspectrum"] = fs->make<TH1F>("chspectrum",";p_{T};counts",
                           hPtBins_.size()-1, &hPtBins_[0]);

   // positively charged hadron spectrum
   hist_["pchspectrum"] = fs->make<TH1F>("pchspectrum",";p_{T};counts",
                           hPtBins_.size()-1, &hPtBins_[0]);

   // positively charged hadron spectrum from u quark
   hist_["pchspectrumU"] = fs->make<TH1F>("pchspectrumU",";p_{T};counts",
                           hPtBins_.size()-1, &hPtBins_[0]);
   // negatively charged hadron spectrum from u quark
   hist_["nchspectrumU"] = fs->make<TH1F>("nchspectrumU",";p_{T};counts",
                           hPtBins_.size()-1, &hPtBins_[0]);
   // positively charged hadron spectrum from d quark
   hist_["pchspectrumD"] = fs->make<TH1F>("pchspectrumD",";p_{T};counts",
                           hPtBins_.size()-1, &hPtBins_[0]);
   // negatively charged hadron spectrum from d quark
   hist_["nchspectrumD"] = fs->make<TH1F>("nchspectrumD",";p_{T};counts",
                           hPtBins_.size()-1, &hPtBins_[0]);
   // positively charged hadron spectrum from gluons
   hist_["pchspectrumG"] = fs->make<TH1F>("pchspectrumG",";p_{T};counts",
                           hPtBins_.size()-1, &hPtBins_[0]);
   // negatively charged hadron spectrum from gluons
   hist_["nchspectrumG"] = fs->make<TH1F>("nchspectrumG",";p_{T};counts",
                           hPtBins_.size()-1, &hPtBins_[0]);

   // 2D histograms of (eta[y],pT) as above
   hist2D_["jetspectrum2D"] = fs->make<TH2F>("jetspectrum2D",";p_{T};counts",
                           etaBins_.size()-1, &etaBins_[0],
                           jetPtBins_.size()-1, &jetPtBins_[0]);
   hist2D_["hspectrum2D"] = fs->make<TH2F>("hspectrum2D",";p_{T};counts",
                           etaBins_.size()-1, &etaBins_[0],
                           hPtBins_.size()-1, &hPtBins_[0]);
   hist2D_["chspectrum2D"] = fs->make<TH2F>("chspectrum2D",";p_{T};counts",
                           etaBins_.size()-1, &etaBins_[0],
                           hPtBins_.size()-1, &hPtBins_[0]);
   hist2D_["pchspectrum2D"] = fs->make<TH2F>("pchspectrum2D",";p_{T};counts",
                           etaBins_.size()-1, &etaBins_[0],
                           hPtBins_.size()-1, &hPtBins_[0]);

   // pt-hat or momentum transfer scale of PYTHIA process
   // this will be 0 for diffractive events
   hist_["qscale"] = fs->make<TH1F>("qscale",";p_{T}-hat;counts",
                           qScalePtBins_.size()-1, &qScalePtBins_[0]);
   hist2D_["ch_qscale2D"] = fs->make<TH2F>("ch_qscale2D",";p_{T}-hat;p_{T} h^{#pm}",
                           qScalePtBins_.size()-1, &qScalePtBins_[0],
                           hPtBins_.size()-1, &hPtBins_[0]);
   hist2D_["jet_qscale2D"] = fs->make<TH2F>("jet_qscale2D",";p_{T}-hat;p_{T} h^{#pm}",
                           qScalePtBins_.size()-1, &qScalePtBins_[0],
                           jetPtBins_.size()-1, &jetPtBins_[0]);

   // leading charged hadron in approximate tracker acceptance (|eta|<2.5)
   hist_["lead_track"] = fs->make<TH1F>("lead_track",";p_{T};counts",
                           hPtBins_.size()-1, &hPtBins_[0]);

   // number of recorded events
   hist_["events"] = fs->make<TH1F>("events",";;events",1,0.,2.);

   // fragmentation function matrix - i.e. 2D histogram of hadrons vs jets
   // also for just charged and positively charged hadrons
   hist2D_["ffmatrix"] = fs->make<TH2F>("ffmatrix",";p_{T} Jet;p_{T} h",
                           hPtBins_.size()-1, &hPtBins_[0],
                           jetPtBins_.size()-1, &jetPtBins_[0]);

   hist2D_["cffmatrix"] = fs->make<TH2F>("cffmatrix",";p_{T} Jet;p_{T} h^{#pm}",
                           hPtBins_.size()-1, &hPtBins_[0],
                           jetPtBins_.size()-1, &jetPtBins_[0]);

   hist2D_["pcffmatrix"] = fs->make<TH2F>("pcffmatrix",";p_{T} Jet;p_{T} h^{+}",
                           hPtBins_.size()-1, &hPtBins_[0],
                           jetPtBins_.size()-1, &jetPtBins_[0]);

}


QCDAnalyzer::~QCDAnalyzer()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void
QCDAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // Pythia does not respect max pt-hat for MB process  
   // where the minimum is 0. In this case we need to 
   // remove events over the pt-hat maximum by hand.
   // We skip and do not count such events.
   // Here it is only coded for the MB_0_to_20 process
   edm::Handle<GenEventInfoProduct> genEvtInfo;
   iEvent.getByLabel("generator", genEvtInfo);
   
   if( pthatMin_ < 1. )
   {
     if( genEvtInfo->qScale() > pthatMax_ ) return;
   }

   // Check if the event is DS and skip if configured
   if( ! isDSEvent(iEvent,iSetup) && onlyDS_ ) return;
  
   Handle<reco::GenParticleCollection> pcol;
   iEvent.getByLabel(genParticleSrc_ ,pcol);

   Handle<reco::GenJetCollection> gcol;
   if( !doFlavor_ ) iEvent.getByLabel(genJetSrc_ ,gcol);

   Handle<reco::JetMatchedPartonsCollection> fcol;
   if( doFlavor_) iEvent.getByLabel(flavorSrc_ ,fcol);

   hist_["events"]->Fill(1);
   hist_["qscale"]->Fill(genEvtInfo->qScale());

   // genjet spectra
   if ( doFlavor_ )
   {
     for( const auto & mjp : *fcol )
     {
        if( ! isInFlavor( mjp.second ) ) continue;
        const reco::Jet *aJet = mjp.first.get();
        if( isInEtaRange( *aJet, jetEtaMin_, jetEtaMax_ ) ) 
        {
           hist_["jetspectrum"]->Fill(aJet->pt());
           hist2D_["jet_qscale2D"]->Fill(genEvtInfo->qScale(), aJet->pt());
        }
        hist2D_["jetspectrum2D"]->Fill(aJet->eta(), aJet->pt());
     }
   }
   else
   {
     for( const auto & jet : *gcol )
     {
       if( isInEtaRange( jet, jetEtaMin_, jetEtaMax_ ) ) 
       {
         hist_["jetspectrum"]->Fill(jet.pt());
         hist2D_["jet_qscale2D"]->Fill(genEvtInfo->qScale(),jet.pt());
       }
       hist2D_["jetspectrum2D"]->Fill(jet.eta(), jet.pt());
     }
   }


//   try mother id
//   Handle<reco::GenParticleCollection> pcol;
//   iEvent.getByLabel(genParticleSrc_ ,pcol);
//  for (unsigned i=0; i<iEvent.size(); ++i) {
     //const Particle &p = iEvent[i];
//     std::cout<<"particle id: " << iEvent[i].id() << std::endl;     
//     std::cout<<"mother id: " << iEvent[i].mother1() << std::endl;
     //std::cout<<"particle id: " << p.id() << std::endl;
     //std::cout<<"mother id: " << (p.mother1()).id() << std::endl;
//  }

   // charged particle spectra and ff matrix
   double lead_track_pt = 0.0;
   for( const auto & h : *pcol )
   {
//	   std::cout << "particle id: " << h.pdgId() << std::endl;
     // skip decayed  particles

     //std::cout<<"particle id: " << h.pdgId() << std::endl;
     //std::cout<<"mother id: " << (h.mother1()).pdgId() << std::endl;
/*
   if(h.pdgId()==-2 && h.numberOfDaughters()>1){
     std::cout<<"start------" << std::endl;
     std::cout<<"particle id: " << h.pdgId() <<"  , N_daughter = "<< h.numberOfDaughters() << std::endl;
     std::cout<<"particle chg " << h.charge() << std::endl;
     for (unsigned int idx=0; idx<h.numberOfDaughters(); idx++){
       std::cout<<"daughter id: " << h.daughter(idx)->pdgId() << std::endl;
       std::cout<<"daughter chg: " << h.daughter(idx)->charge() << std::endl;
     }
     std::cout<<"------end" << std::endl;
   }
*/
/*
   if(h.pdgId()==211 && h.numberOfMothers()>0){
     std::cout<<"start------" << std::endl;
     std::cout<<"particle id: " << h.pdgId() << "  , N_mother = "<< h.numberOfMothers() << std::endl;
     for (unsigned int idx=0; idx<h.numberOfMothers(); idx++){
       std::cout<<"mother id: " << h.mother(idx)->pdgId() << std::endl;
     }
     std::cout<<"------end" << std::endl;
   }
*/

   if(h.pdgId()==2 && h.numberOfDaughters()>1){ // u quark
     for (unsigned int idx=0; idx<h.numberOfDaughters(); idx++){
       if( h.daughter(idx)->eta()<=hEtaMax_ && h.daughter(idx)->eta()>= hEtaMin_ ){
         if(h.daughter(idx)->charge() > 0) hist_["pchspectrumU"]->Fill(h.daughter(idx)->pt());
         if(h.daughter(idx)->charge() < 0) hist_["nchspectrumU"]->Fill(h.daughter(idx)->pt());
       }
     }
   }
   if(h.pdgId()==1 && h.numberOfDaughters()>1){ //d quark
     for (unsigned int idx=0; idx<h.numberOfDaughters(); idx++){
       if( h.daughter(idx)->eta()<=hEtaMax_ && h.daughter(idx)->eta()>= hEtaMin_ ){
         if(h.daughter(idx)->charge() > 0) hist_["pchspectrumD"]->Fill(h.daughter(idx)->pt());
         if(h.daughter(idx)->charge() < 0) hist_["nchspectrumD"]->Fill(h.daughter(idx)->pt());
       }
     }
   }
   if(h.pdgId()==21 && h.numberOfDaughters()>1){ //gluons
     for (unsigned int idx=0; idx<h.numberOfDaughters(); idx++){
       if( h.daughter(idx)->eta()<=hEtaMax_ && h.daughter(idx)->eta()>= hEtaMin_ ){
         if(h.daughter(idx)->charge() > 0) hist_["pchspectrumG"]->Fill(h.daughter(idx)->pt());
         if(h.daughter(idx)->charge() < 0) hist_["nchspectrumG"]->Fill(h.daughter(idx)->pt());
       }
     }
   }


     if( h.status() != 1  ) continue;

     // update leading track
     if( h.charge() != 0 && fabs(h.eta()) < 2.5 && lead_track_pt < h.pt() )
       lead_track_pt = h.pt(); 

     if( doSpecies_ && ! isInSpecies(h) ) continue;
       
     hist2D_["hspectrum2D"]->Fill(h.eta(),h.pt());
     if( h.charge() != 0 ) hist2D_["chspectrum2D"]->Fill(h.eta(),h.pt());
     if( h.charge() > 0 ) hist2D_["pchspectrum2D"]->Fill(h.eta(),h.pt());

     if( isInEtaRange(h, hEtaMin_, hEtaMax_) )
     {
       hist_["hspectrum"]->Fill(h.pt());
       if( h.charge() != 0 ) hist_["chspectrum"]->Fill(h.pt());
       if( h.charge() != 0 ) hist2D_["ch_qscale2D"]->Fill(genEvtInfo->qScale(),h.pt());
       if( h.charge() > 0 ) hist_["pchspectrum"]->Fill(h.pt());

       // associate with the highest-pt jet for which 
       // the track is found in the jet cone
     
       double maxPtJet = 0.;
       if( doFlavor_)
       {
         for( const auto & mjp : *fcol )
         {
           const reco::Jet *aJet = mjp.first.get();
           if ( !isInEtaRange(*aJet,jetEtaMin_,jetEtaMax_) ) continue;
           if( !isInFlavor( mjp.second ) ) continue;
           double dr = deltaR(*aJet,h);
           if( dr < jetRadius_ && aJet->pt() > maxPtJet)
             maxPtJet = aJet->pt();
         }
       }
       else
       {
         for( const auto & jet : *gcol )
         {
           if( ! isInEtaRange( jet, jetEtaMin_, jetEtaMax_ ) ) continue; 
           double dr = deltaR(jet,h);
           if( dr < jetRadius_ && jet.pt() > maxPtJet)
             maxPtJet = jet.pt();
         }
       }
       hist2D_["ffmatrix"]->Fill( maxPtJet, h.pt());
       if( h.charge() != 0 ) hist2D_["cffmatrix"]->Fill( maxPtJet, h.pt());
       if( h.charge() > 0 ) hist2D_["pcffmatrix"]->Fill( maxPtJet, h.pt());
     }
   }

   hist_["lead_track"]->Fill(lead_track_pt);
}


bool 
QCDAnalyzer::isInEtaRange( const reco::Candidate& c, double etaMin, double etaMax )
{
     if( useRapidity_ == false &&  c.eta() <= etaMax && c.eta() >= etaMin )
       return true;
     if( useRapidity_ == true &&  c.y() <= etaMax && c.y() >= etaMin )
       return true;
     return false;
}

bool
QCDAnalyzer::isInFlavor( const reco::MatchedPartons & aMatch )
{
  int flavor = 0;
  if( aMatch.heaviest().isNonnull() )
  {
    flavor = aMatch.heaviest().get()->pdgId();
  }
  for( const auto inFlavor : flavorId_ )
  { if( flavor == inFlavor ) return true; }
 
  return false;
}

bool 
QCDAnalyzer::isInSpecies( const reco::GenParticle & h )
{
  for( const auto inSpecies : speciesId_ )
  { if( h.pdgId() == inSpecies ) return true; }

  return false;  
}

bool 
QCDAnalyzer::isDSEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

     using namespace edm;

     bool posDS = false; bool negDS = false;

     edm::ESHandle<ParticleDataTable> particleDataTable_;
     iSetup.getData(particleDataTable_);

     Handle<reco::GenParticleCollection> gcol;
     iEvent.getByLabel(genParticleSrc_, gcol);
     for( const auto & gen : * gcol )
     {
       // see if genpartice counts for DS
       HepPDT::ParticleID particleID(gen.pdgId());
       if (particleID.isValid())
       {
         ParticleData const * particleData = particleDataTable_->particle(particleID);
         if (particleData)
         { 
           double tau =  particleDataTable_->particle(particleID)->lifetime();
           if ( tau  > 1e-18 || tau == 0.0 )
           {
             if( gen.energy() > 3.0 && gen.eta() > 3.0 && gen.eta() < 5.0 ) posDS = true;
             if( gen.energy() > 3.0 && gen.eta() < -3.0 && gen.eta() > -5.0 ) negDS = true;
           }
         }
       }
     }

     if( posDS && negDS ) return true;
     else return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
QCDAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
QCDAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
QCDAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
QCDAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
QCDAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
QCDAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
QCDAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(QCDAnalyzer);
