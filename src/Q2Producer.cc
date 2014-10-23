// -*- C++ -*-
//
// Package:    Q2Producer
// Class:      Q2Producer
// 
/**\class Q2Producer Q2Producer.cc QWAna/Q2Producer/src/Q2Producer.cc

Description: HeavyIon q2

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Quan Wang
//         Created:  Thu Oct 23 11:01:36 CEST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

//
// class declaration
//

class Q2Producer : public edm::EDProducer {
	public:
		explicit Q2Producer(const edm::ParameterSet&);
		~Q2Producer();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() ;
		virtual void produce(edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		virtual void beginRun(edm::Run&, edm::EventSetup const&);
		virtual void endRun(edm::Run&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

		// ----------member data ---------------------------
		edm::InputTag vtxCollection_;
		edm::InputTag caloCollection_;
		edm::InputTag trackCollection_;
		bool	useweight_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
Q2Producer::Q2Producer(const edm::ParameterSet& iConfig)
{
	vtxCollection_  = iConfig.getParameter<edm::InputTag>("vtxCollection_");
	caloCollection_  = iConfig.getParameter<edm::InputTag>("caloCollection_");
	trackCollection_  = iConfig.getParameter<edm::InputTag>("trackCollection_");
	useweight_ = iConfig.getUntrackedParameter<bool>("useweight_",true);
	produces<std::vector<double>>("recolevel");
	//register your products
	/* Examples
	   produces<ExampleData2>();

	//if do put with a label
	produces<ExampleData2>("label");

	//if you want to put into the Run
	produces<ExampleData2,InRun>();
	*/
	//now do what ever other initialization is needed

}


Q2Producer::~Q2Producer()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
	void
Q2Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	edm::Handle<reco::VertexCollection> vertexCollection;
	iEvent.getByLabel(vtxCollection_,vertexCollection);
	
	Handle<CaloTowerCollection> calotower;
	iEvent.getByLabel(caloCollection_,calotower);

	double q2trklowpt = 0;
	double q2trkhighpt = 0;
	double q2trk = 0;
	double q2HFp = 0;
	double q2HFm = 0;
	double q2HF = 0;

	double q3trklowpt = 0;
	double q3trkhighpt = 0;
	double q3trk = 0;
	double q3HFp = 0;
	double q3HFm = 0;
	double q3HF = 0;

	// calo
	if ( calotower.isValid() ) {
		double s2p = 0;
		double s3p = 0;
		double s2m = 0;
		double s3m = 0;
		double s2 = 0;
		double s3 = 0;

		double c2p = 0;
		double c3p = 0;
		double c2m = 0;
		double c3m = 0;
		double c2 = 0;
		double c3 = 0;

		double wp = 0;
		double wm = 0;
		double w = 0;
		for ( const auto calo : *calotower) {
			double eta = calo.eta();
			double phi = calo.phi();
			double weight = calo.emEt() + calo.hadEt();

			double lsin2 = weight*sin(2*phi);
			double lsin3 = weight*sin(3*phi);
			double lcos2 = weight*cos(2*phi);
			double lcos3 = weight*cos(3*phi);

			if ( eta > 3. and eta < 5. ) {
				s2p += lsin2;
				c2p += lcos2;
				s2 += lsin2;
				c2 += lcos2;

				s3p += lsin3;
				c3p += lcos3;
				s3 += lsin3;
				c3 += lcos3;

				wp += weight;
				w += weight;
			}
			if ( eta < -3. and eta > -5. ) {
				s2m += lsin2;
				c2m += lcos2;
				s2 += lsin2;
				c2 += lcos2;

				s3m += lsin3;
				c3m += lcos3;
				s3 += lsin3;
				c3 += lcos3;

				wm += weight;
				w += weight;
			}
		}
		q2HFp = sqrt( s2p*s2p + c2p*c2p ) / wp;
		q2HFm = sqrt( s2m*s2m + c2m*c2m ) / wm;
		q2HF = sqrt( s2*s2 + c2*c2 ) / w;

		q3HFp = sqrt( s3p*s3p + c3p*c3p ) / wp;
		q3HFm = sqrt( s3m*s3m + c3m*c3m ) / wm;
		q3HF = sqrt(  s3 *s3  + c3 *c3  ) / w;
	}

	// track

	Handle<reco::TrackCollection> tracks;
	iEvent.getByLabel(trackCollection_, tracks);
	if ( tracks.isValid() ) {
		double s2p = 0;
		double s3p = 0;
		double s2m = 0;
		double s3m = 0;
		double s2 = 0;
		double s3 = 0;

		double c2p = 0;
		double c3p = 0;
		double c2m = 0;
		double c3m = 0;
		double c2 = 0;
		double c3 = 0;

		double wp = 0;
		double wm = 0;
		double w = 0;
	
		for ( const auto trk: *tracks) {
			double eta = trk.eta();
			double phi = trk.phi();
			double weight = trk.pt();
			if ( fabs(eta) > 2.4 ) continue;
			if ( weight > 2.5 ) weight = 2.;

			double lsin2 = weight*sin(2*phi);
			double lsin3 = weight*sin(3*phi);
			double lcos2 = weight*cos(2*phi);
			double lcos3 = weight*cos(3*phi);

			if ( weight > 1. ) {
				s2p += lsin2;
				c2p += lcos2;
				s2 += lsin2;
				c2 += lcos2;

				s3p += lsin3;
				c3p += lcos3;
				s3 += lsin3;
				c3 += lcos3;

				wp += weight;
				w += weight;
			} else {
				s2m += lsin2;
				c2m += lcos2;
				s2 += lsin2;
				c2 += lcos2;

				s3m += lsin3;
				c3m += lcos3;
				s3 += lsin3;
				c3 += lcos3;

				wm += weight;
				w += weight;
			}
		}
		q2trklowpt = sqrt( s2m*s2m + c2m*c2m ) /wm;
		q2trkhighpt = sqrt( s2p*s2p + c2p*c2p ) /wp;
		q2trk = sqrt( s2*s2 + c2*c2 ) /w;

		q3trklowpt = sqrt( s3m*s3m + c3m*c3m ) /wm;
		q3trkhighpt = sqrt( s3p*s3p + c3p*c3p ) /wp;
		q3trk = sqrt( s3*s3 + c3*c3 ) /w;
	}

	std::auto_ptr<std::vector<double>> ptr(new std::vector<double>);
	ptr->push_back(q2trklowpt);
	ptr->push_back(q2trkhighpt);
	ptr->push_back(q2trk);
	ptr->push_back(q2HFp);
	ptr->push_back(q2HFm);
	ptr->push_back(q2HF);
	ptr->push_back(q3trklowpt);
	ptr->push_back(q3trkhighpt);
	ptr->push_back(q3trk);
	ptr->push_back(q3HFp);
	ptr->push_back(q3HFm);
	ptr->push_back(q3HF);

	iEvent.put(ptr, "recolevel");

}

// ------------ method called once each job just before starting event loop  ------------
	void 
Q2Producer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Q2Producer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
	void 
Q2Producer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
	void 
Q2Producer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
Q2Producer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
Q2Producer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Q2Producer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Q2Producer);
