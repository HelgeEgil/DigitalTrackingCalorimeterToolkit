#ifndef Track_h
#define Track_h

#include <TClonesArray.h>
#include <vector>

#include "Classes/Cluster/Cluster.h"
// #include "Classes/Track/conversionFunctions.h"

class TGraph;
class TGraphErrors;

namespace DTC {
class Clusters;

const Int_t MaxTrackLength = 250;

class Track : public TObject {

   private:
      TClonesArray track_;
      Float_t  fitEnergy_;
      Float_t  fitRange_;
      Float_t  fitScale_;
      Float_t  fitError_;
      Float_t  fitChi2_;
      Float_t  fitSigma_;
      Float_t  precomputeFactor;
      Bool_t   isIncomplete_;

   public:
      Track();
      Track(Cluster *cluster);

      virtual ~Track(); 


      // ROOT & I/O
      virtual Cluster * At(Int_t i)                      { return ((Cluster*) track_.At(i)); }
      virtual Cluster * Last()                           { return ((Cluster*) track_.At(GetEntriesFast()-1)); }
      virtual Int_t     GetEntriesFast()                 { return track_.GetEntriesFast(); }
      virtual Int_t     GetEntries()                     { return track_.GetEntries(); }
      virtual void      Compress()                       { track_.Compress(); }
      virtual void      Clear(Option_t * option = "")    { track_.Clear(option); }
      void              sortTrack()                      { track_.Sort(); }

      // Sorting
      Int_t             Compare(const TObject* obj) const;
      Bool_t            isSortable() const { return kTRUE; }

      // Add and remove clusters
      virtual TObject * removeClusterAt(Int_t i)         { return track_.RemoveAt(i); }
      virtual void      removeCluster(Cluster *c)        { track_.Remove((TObject*) c); }
      virtual void      setTrack(Track *copyTrack, Int_t startOffset = 0); // copy whole track
      virtual void      appendCluster(Cluster *copyCluster, Int_t startOffset = 0); // copy cluster
      virtual void      appendPoint(Float_t x, Float_t y, Int_t layer, Int_t size = -1, Int_t eventID = -1, Bool_t isSecondary = false, Int_t PDG = 0);
      void              removeNANs();

      // Getters and setters
      virtual Float_t   getX(Int_t i)                    { return At(i)->getX(); }
      virtual Float_t   getY(Int_t i)                    { return At(i)->getY(); }
      virtual Int_t     getLayer(Int_t i)                { return At(i)->getLayer(); }
      Float_t           getXmm(Int_t i)                  { return At(i)->getXmm(); }
      Float_t           getYmm(Int_t i)                  { return At(i)->getYmm(); }
      Float_t           getLayermm(Int_t i)              { return At(i)->getLayermm(); }
      virtual Int_t     getSize(Int_t i)                 { return At(i)->getSize(); }
      virtual Int_t     getEventID(Int_t i)              { return At(i)->getEventID(); }
      virtual Int_t     getClusterSizeError(Int_t i)     { return At(i)->getError(); } 
      virtual Float_t   getDepositedEnergy(Int_t i, Bool_t c = false)      { return At(i)->getDepositedEnergy(c); }
      virtual Float_t   getDepositedEnergyError(Int_t i, Bool_t c = false) { return At(i)->getDepositedEnergyError(c); }
      Bool_t            isUsed(Int_t i)                  { return At(i)->isUsed(); }
      Bool_t            isSecondary(Int_t i)             { return At(i)->isSecondary(); }
      void              setIncomplete(Bool_t inc)        { isIncomplete_ = inc; }
      Bool_t            isIncomplete()                   { return isIncomplete_; }
      
      // TRACK PROPERTIES - Event IDs
      Int_t             getModeEventID(); // (vs. median)
      Bool_t            isOneEventID();
      Bool_t            isFirstAndLastEventIDEqual();

      // TRACK PROPERTIES - General track properties
      // trackProperties.C
      Int_t             getFirstLayer();
      Int_t             getLastLayer();
      Bool_t            hasLayer(Int_t layer);
      Int_t             getIdxFromLayer(Int_t i);
      Bool_t            doesTrackEndAbruptly();
      Float_t           getRiseFactor(Int_t n = 5);
      Int_t             getNMissingLayers();
      void              propagateSecondaryStatus();
      Float_t           getAverageDepositedEnergy(Int_t fromIdx = 0, Int_t toIdxExclusive = -1);

      // TRACK PROPERTIES - Ranges and energies
      // trackRangeCalculations.C
      Float_t           getTrackLengthmm();
      Float_t           getTrackLengthmmAt(Int_t i);  // TL between i-1 and i in mm
      Float_t           getTrackLengthWEPLmmAt(Int_t i);
      Float_t           getRangemm();
      Float_t           getRangemmAt(Int_t i);
      Float_t           getRangeWEPLAt(Int_t i);
      Float_t           getWEPL();
      Float_t           getEnergy();
      Float_t           getEnergyStraggling();

      // TRACK PROPERTIES - Angles
      // trackAngleCalculations.C
      Float_t           getSlopeAngleAtLayer(Int_t i);
      Float_t           getSlopeAngleBetweenLayers(Int_t i);
      Float_t           getSlopeAngle();
      Float_t           getSlopeAngleChangeBetweenLayers(Int_t i);
      Float_t           getSlopeAngleDifferenceSum();
      Float_t           getSlopeAngleDifferenceSumInTheta0();
      Float_t           getSinuosity();
      Float_t           getProjectedRatio();
      Float_t           getMaximumSlopeAngleChange();
      Float_t           getAbsorberLength(Int_t i);

      // TRACK PROPERTIES - Cluster properties
      // trackClusterProperties.C
      Int_t             getClusterIdx(Cluster * cluster);
      Int_t             getClusterIdx(Float_t x, Float_t y, Int_t layer);
      Bool_t            isClusterInTrack(Cluster * cluster);
      Int_t             getClusterFromLayer(Int_t layer);
      Clusters        * getConflictClusters();
      Int_t             getNumberOfConflictClusters();
      Bool_t            isUsedClustersInTrack();
      Float_t           getAverageCS();
      Float_t           getAverageCSLastN(Int_t i);
      Float_t           getMeanSizeToIdx(Int_t i);
      Float_t           getStdSizeToIdx(Int_t toIdx);

      // TRACK PROPERTIES - Extrapolations
      // trackExtrapolations.C
      Cluster         * getInterpolatedClusterAt(Int_t layer);
      Cluster         * getExtrapolatedClusterAt(Float_t mmBeforeDetector);
      vector<Float_t>   getLateralDeflectionFromExtrapolatedPosition(Int_t layer);
      void              extrapolateToLayer0();

      // TRACK PROPERTIES - Track fitting, scoring
      // trackFitting.C
      TGraphErrors    * doTrackFit(Bool_t isScaleVariable = false, Bool_t useTrackLength = kUseCSDA);
//      TGraphErrors    * doRangeFit(Bool_t isScaleVariable = false);
      Float_t           getFitParameterRange();
      Float_t           getFitParameterScale();
      Float_t           getFitParameterError();
      Float_t           getFitParameterChiSquare();
      Float_t           getFitParameterSigma();
      Float_t           getTrackScore();

      ClassDef(Track,5)
};
}

ostream& operator<<(ostream &os, DTC::Track& t);

#endif
