#include <iostream>
#include <vector>
#include <math.h>

#include "Classes/Cluster/Node.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Track/Track.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

using namespace DTC;

Node::Node(Node *parent, Cluster *connectedCluster, Float_t score) {
   parent_ = parent;
   connectedCluster_ = connectedCluster;
   score_ = score;
   isExplored_ = false;
   lastChild_ = 0;
   for (Int_t i=0; i<nChildrenInNode; i++) children_[i] = nullptr;
}

Node::~Node() {
   if (getNChildren() > 0) {
      for (Int_t i=0; i<nChildrenInNode; i++) { // max 5 children per node
         if (getChild(i)) {
            delete children_[i]; // illegal in 45 states
         }
      }
   }
}

Cluster * Node::getParentCluster() {
   Node * parent = getParent();
   if (!parent) return nullptr;
   else return parent->getCluster();
}


void Node::addChild(Node *node) {
   Int_t i = getFirstEmptyChild();
   if (i<0) {
      Int_t worstChild = getWorstChild();
      Float_t worstScore = getChild(worstChild)->getScore();
      if (worstScore > node->getScore()) {
         children_[worstChild] = node; 
      }
   }
   else {
      children_[i] = node;
   }
}

Int_t Node::getFirstEmptyChild() {
   for (Int_t i=0; i<nChildrenInNode; i++) {
      if (!getChild(i)) return i;
   }
   return -1;
}

Int_t Node::getWorstChild() {
   Int_t worstChild = -1;
   Float_t worstScore = -1;

   for (Int_t i=0; i<getNChildren(); i++) {
      if (getChild(i)->getScore() > worstScore) {
         worstChild = i;
         worstScore = getChild(i)->getScore();
      }
   }
   return worstChild;
}

Cluster * Node::getExtrapolatedCluster() {
   if (!getParent()) return nullptr;
   Cluster *p2 = getCluster();
   Cluster *p1 = getParentCluster();

   return getTrackExtrapolationCluster(p1, p2);
}

Cluster * Node::getNextUnexploredCluster() {
   for (Int_t i=lastChild_; i<nChildrenInNode; i++) {
      if (!getChild(i)) continue;

      if (!getChild(i)->isExplored()) {
         lastChild_ = i;
         return getChild(i)->getCluster();
      }
   }
}

Int_t Node::getNChildren() {
   Int_t nClusters = 0;
   for (Int_t i=0; i<nChildrenInNode; i++) {
      if (getChild(i)) nClusters++;
   }
   return nClusters;
}

void Node::getEndNodes(vector<Node*> *endNodes) {
   if (!endNodes) endNodes = new vector<Node*>;
   
   Int_t nChildren = getNChildren();

   if (nChildren == 0) {
      endNodes->push_back(this);
   }

   else {
      for (Int_t i=0; i<nChildren; i++) {
         if (getChild(i)) {
            getChild(i)->getEndNodes(endNodes);
         }
      }
   }
}

void Node::getUnexploredEndNodes(vector<Node*> *endNodes) {
   if (!endNodes) endNodes = new vector<Node*>;

   Int_t nChildren = getNChildren(); 

   if (nChildren == 0 && !isExplored()) {
      endNodes->push_back(this);
   }
   else if (nChildren > 0) {
      for (Int_t i=0; i<nChildren; i++) {
         if (getChild(i)) {
            getChild(i)->getUnexploredEndNodes(endNodes);
         }
      }
   }
}

void  Node::getAllNodes(vector<Node*> *allNodes) {
   if (!allNodes) {
      allNodes = new vector<Node*>;
   }
      
   allNodes->push_back(this);
   
   Int_t nChildren = getNChildren();

   if (nChildren > 0) {
      for (Int_t i=0; i<nChildren; i++) {
         if (getChild(i)) {
            getChild(i)->getAllNodes(allNodes);
         }
      }
   }
}

Int_t Node::getClusterCount() {
   vector<Node*> * allNodes = new vector<Node*>;
   getAllNodes(allNodes);
   return allNodes->size();
}

Float_t Node::getNodeAngle(Cluster *p3) {
   Cluster *p1 = getParentCluster();
   Cluster *p2 = getCluster();

   if (!p1) return getDotProductAngle(p2, p2, p3);
   else     return getDotProductAngle(p1, p2, p3);
}

Node* Node::getBestNode() {
   vector<Node*> * endNodes = new vector<Node*>;
   Float_t         bestScore = 1e7;
   Float_t         score;
   Node          * bestNode = nullptr;
   Node          * thisNode = nullptr;

   getEndNodes(endNodes);

   for (UInt_t i=0; i<endNodes->size(); i++) {
      thisNode = endNodes->at(i);
      score = thisNode->getScore();
      if (score < bestScore) {
         bestScore = score;
         bestNode = thisNode;
      }
   }

   return bestNode;
}

Float_t Node::getNextScore(Cluster *c) {
   Float_t newAngle = getNodeAngle(c);
   return sqrt(pow(getScore(), 2) + pow(newAngle, 2));
}

Track* Node::getTrackFromBestNode(Node* bestNode) {
   vector<Node*>  reverseTrack;
   Node         * nextSegment = bestNode;
   Track        * bestTrack = new Track();


   Int_t n = 0;
   while (nextSegment && n<100) {
      reverseTrack.push_back(nextSegment);
      nextSegment = (Node*) nextSegment->getParent();
      n++;
   }
   
   while (reverseTrack.size() > 0) {
      bestTrack->appendCluster(reverseTrack.back()->getCluster());
      reverseTrack.pop_back();
   }

   return bestTrack;
}

Track* Node::getBestTrack() {
   Node *bestNode = getBestNode();
   return getTrackFromBestNode(bestNode);
}
