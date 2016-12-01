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

#include <stdlib.h>
#include <stdio.h>

class LinkedQueue{
public:
  //data members
  Qnode *head;
  Qnode *rear;

  //function members
  //constructor
  LinkedQueue(){
    head=NULL;
    rear=NULL;
  }//constructor

  int dequeue(){
    if(empty()){
      printf("tried to dequeue from an empty queue\n");
      exit(1);
    }
    int retval=head->info;
    Qnode *aux = head;
    head=head->next;
    delete aux;
    return retval;
  }//dequeue

  void enqueue(int x){
    if(full()){ //for now full is always false
      printf("tried to enqueue into a full queue\n");
      exit(1);
    }
    Qnode *aux = new Qnode(x,NULL);
    if(empty()){
      head=aux;
      rear=aux;
    }
    else{
      rear->next=aux;
      rear=aux;
    }
  }//enqueue

  int full(){
    return 0; //false
  }//full

  int empty(){
    return (head==NULL);
  }//empty
};
