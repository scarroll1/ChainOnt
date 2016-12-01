/*
 * The contents of this file are subject to the Mozilla Public
 * License Version 1.1 (the "License"); you may not use this file
 * except in compliance with the License. You may obtain a copy of
 * the License at http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS
 * IS" basis, WITHOUT WARRANTY OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * rights and limitations under the License.
 *
 * The Original Code is ChainOnt Protein Classifier.
 *
 * The Initial Developer of the Original Code is Steven D. Carroll
 * Copyright (C) 2006.  All Rights Reserved.
 *
 * Contributor(s):
 *   Steven D. Carroll <carroll@genome.chop.edu>
 *   Vladimir Pavlovic <vladimir@cs.rutgers.edu>
 */

class Node{ //Graph node - adj list
public:
  //data members
  int nummsgs; //# of messages from one node to other
  long double *m; //array of messages - variable size
  long double s; //similarity
  int id; 
  bool gridedge;
  Node *next;

  //function members
  //default constructor
  Node(){
    m=NULL;
    nummsgs=1; //arbitrary value for default constructor only
    id=-1;
    next=NULL;
  }//default constructor

  //constructor
  Node(int nummsgs,long double *m,long double s,int id){
    this->nummsgs=nummsgs;
    this->m=new long double[nummsgs];
    for(int p=0; p<nummsgs; p++)
      this->m[p]=m[p];
    this->s=s;
    this->id=id;
    this->next=NULL;
    this->gridedge = false;
  }//constructor

  Node(int nummsgs,long double *m,long double s,int id,bool gridedge){
    this->nummsgs=nummsgs;
    this->m=new long double[nummsgs];
    for(int p=0; p<nummsgs; p++)
      this->m[p]=m[p];
    this->s=s;
    this->id=id;
    this->next=NULL;
    this->gridedge = gridedge;
  }//constructor

};//Node
   
