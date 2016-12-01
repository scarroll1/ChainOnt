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

class Qnode{
public:
  //data members
  int info;
  Qnode *next;
  //function members
  //default constructor
  Qnode(){
    info = -1;
    next = NULL;
  }//default const
  //constructor
  Qnode(int id,Qnode *next){
    info=id;
    this->next=next;
  }
};
