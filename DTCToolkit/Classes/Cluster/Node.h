#ifndef Node_h
#define Node_h

#include <iostream>
#include <vector>

#include <TObject.h>
#include <math.h>

#include "GlobalConstants/Constants.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Track/Track.h"

using namespace std;

namespace DTC {
class Node : public TObject {
private:
   Node    * parent_;
   Node*     children_[nChildrenInNode]; // max 5 children per node
   Int_t     lastChild_;
   Cluster * connectedCluster_;
   Float_t   score_;
   Bool_t    isExplored_;

public:
   Node(Node *parent, Cluster *connectedCluster, Float_t score);
   virtual ~Node();
   
   // inline getters and setters
   Node         * getParent()          { return parent_; }
   Node         * getChild(Int_t i)    { return (Node*) children_[i]; }
   Cluster      * getCluster()         { return connectedCluster_; }
   Float_t        getScore()           { return score_; }
   void           markExplored()       { isExplored_ = true; }
   Bool_t         isExplored()         { return isExplored_; }
   Int_t          getFirstEmptyChild();
   Int_t          getWorstChild();
   Int_t          getBestChild();
   void           addChild(Node *child);
   void           removeChild(Int_t i);
   void           deleteNodeTree();

   Cluster      * getParentCluster();
   Cluster      * getExtrapolatedCluster();
   Cluster      * getNextUnexploredCluster();
   Int_t          getNChildren();
   Int_t          getClusterCount();
   void           getEndNodes(vector<Node*> *endNodes = nullptr);
   void           getUnexploredEndNodes(vector<Node*> *endNodes = nullptr);
   void           getAllNodes(vector<Node*> *allNodes = nullptr);

   Float_t        getNodeAngle(Cluster *c);
   Float_t        getNextScore(Cluster *c);

   Node         * getBestNode();
   Track        * getTrackFromBestNode(Node * bestNode);
   Track        * getBestTrack();

   ClassDef(Node,1)
};
}

#endif
