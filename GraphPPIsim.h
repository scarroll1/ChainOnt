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
#include <math.h>
#include "hash_set.h"

class Graph{
  typedef hash_set<int> int_hash_set;
public:
  //data members
  int numnodes;
  Node *neighbors;
  LinkedQueue jci;
  int *sendorder;
  long double **phi; //2-D array - must be double for function nodes/cond probs
  long double **belief;
  long double maxdiff;  
  int *card; //array of cardinalities of nodes
  int *numpars; //array of number of parents - only used for function nodes, ie even #'s
  int **parents;
  int numterms; //number of GO nodes
  int numprots; //numprots - numnodes (above) is numterms*2*numprots
  long double *sdist; //shortest distance to any evidence protein
  long double threshold; //distance above which protein not included
  bool *visited;
  bool *jcibool;
  bool *evidpos;
  bool *evidneg;
  int numreachable;
  long double maxsim;
  long double **mult;
  //function members
  //constructor
  Graph(int numnodes,bool bayes){
    maxsim = 1; 
    this->numnodes=numnodes;
    visited = new bool[numnodes];
    jcibool = new bool[numnodes];
    evidpos = new bool[numnodes];
    evidneg = new bool[numnodes];
    numpars = new int[numnodes];
    parents = new int*[numnodes];
    neighbors = new Node[numnodes];
    if(neighbors==NULL){
      printf("Out of memory - Graph constructor, neighbors\n");
      exit(1);
    } 
    sendorder = new int[numnodes];
    if(sendorder==NULL){
      printf("Out of memory - Graph constructor, sendorder\n");
      exit(1);
    }
    long double *msgs = new long double[1];
    msgs[0] = -1.0;
    for(int i=0; i<numnodes; i++){
      Node *temp = new Node(1,msgs,-1.0,-1);
      neighbors[i] = *temp;
      sendorder[i]=-1;
      visited[i] = false;
      jcibool[i] = false;
      evidpos[i] = false;
      evidneg[i] = false;
    }//for i
    belief = new long double*[numnodes];
    if(belief==NULL){
      printf("Out of memory - Graph constructor, belief\n");
      exit(1);
    }
    phi = new long double*[numnodes];
    if(phi==NULL){   
      printf("Out of memory - Graph constructor, phi\n");   
      exit(1);
    }
    if(!bayes){
      for(int i=0; i<numnodes; i++){
        phi[i] = new long double[2];
        belief[i] = new long double[2];
      }//for i
    }//if(!bayes)
    else{
      for(int i=0; i<numnodes/2; i++){
        phi[2*i+1] = new long double[2];
        belief[2*i+1] = new long double[2];
        numpars[2*i] = 0;
      }//for i
    }//else
    int_hash_set evid;
    int_hash_set evidpos;
    int_hash_set evidneg;
    int_hash_set sentout;
    int_hash_set jcihash;
    LinkedQueue jci;
    card = new int[numnodes];
    printf("Graph const cp8 - numnodes = %d\n",numnodes);
  }//constructor

  void add(int prot,Node *n){
    if(neighbors[prot].id==-1){
      neighbors[prot]=*n;
    }//if
    else{
      Node *cur=&neighbors[prot];
      while(cur->next!=NULL){
        cur=cur->next;
      }//while
      cur->next=n;
    }//else
  }//add

  Node* findnode(int node1,int node2){
    Node *cur = &neighbors[node1];
    while(cur!=NULL){
      if((cur->id)==node2)
        break;
      cur=cur->next;
    }//while
    return cur; //returns NULL if not there
  }//findnode

  bool contains(int_hash_set ihs,int check){
    int_hash_set::const_iterator it = ihs.find(check);
    if(it==ihs.end())
      return false;
    return true;
  }//contains

  void definerestofphi(){ //for BN nodes, ie term nodes
    int limit = numnodes/2;
    for(int i=0;i<limit;i++){
      int check=2*i+1;
      if(evidpos[check]){
        phi[check][0]=0;
        phi[check][1]=1;
      }//if
      else if(evidneg[check]){
        phi[check][0]=1;
        phi[check][1]=0;
      }//else if
      else{
        phi[check][0]=1;
        phi[check][1]=1;
      }//else
    }//for i
  }//definerestofphi

  void sendmsgsout(int i,FILE *status){ //i in Yed's msg passing
    for(Node *cur=&neighbors[i];cur!=NULL;cur=cur->next){
      int j=cur->id; //j in Yed's msg passing
      if(!evidpos[j]&&!evidneg[j]){
        long double **mij = new long double*[card[i]];
        for(int ct=0; ct<card[i]; ct++)
          mij[ct] = new long double[card[j]];
        if((i%2==1)&&(j%2==1)){ //2 term nodes, same term, diff't orfs - find out which term
          int inode = i%(numterms*2);
          int term = inode/2; //integer division
          int jnode = j%(numterms*2);
          if(jnode!=inode){
            printf("jnode != inode");
            exit(EXIT_FAILURE);
          }
          long double simmult;
          long double gridmult;
          for(int x=0; x<card[i]; x++){
            for(int y=0; y<card[j]; y++){
              if(x==y)
                simmult=cur->s; //will be 0/1 for function nodes
              else
                simmult=maxsim-cur->s; //maxsim not always 1
              int ind = 2*x + y; //binary -> base 10
              if(cur->gridedge)
                gridmult = mult[term][ind];
              else
                gridmult = 1;
              mij[x][y]=(phi[i][x])*gridmult*simmult;
            }//for y
          }//for x
        }//if
        else{ //bayes is true and at least one of i and j is even, ie function node - must be exactly one
          long double mult;
          int funid;
          int termid;
          if(i%2==0){
            funid = i;
            termid = j;
          }
          else{ 
            funid = j;
            termid = i;
          }
          int parentpos = -1; //if still -1 after next loop, it must be child
          for(int p=0; p<numpars[funid]; p++)
            if(termid==parents[funid][p]){
              parentpos = p;
              break;
            }
          for(int x=0; x<card[i]; x++) //card[i] = 2^(numparents+1)
            for(int y=0; y<card[j]; y++){ //card[j] = 2, or vice-versa
              int limit; //relevant position in binary number
              if(parentpos==-1) //child, highest order binary digit
                limit = numpars[funid]; //child is highest order bin
                                        //digit, which corresponds 
                                        //to 2^numparents
              else
                limit = parentpos; 
              int exp = 1; //will be 2^a
              for(int ct=0; ct<limit; ct++)
                exp=exp*2;
              int counter;
              int bindig;
              if(funid==i){
                counter = x;
                bindig = y;
              }
              else{
                counter=y;
                bindig=x;
              }
              if((counter/exp)%2==bindig) //extracts the bin digit
                mult = 1;
              else
                mult = 0;
              mij[x][y]=phi[i][x]*mult;
            }//for y
        }//else - bayes is true
        for(Node *cur2=&neighbors[i];cur2!=NULL;cur2=cur2->next){
          if(cur2->id!=j){
            for(int x=0; x<card[i]; x++){
              for(int y=0; y<card[j]; y++){
                mij[x][y]=mij[x][y]*cur2->m[(x+card[cur2->id])];
              }//for y                  
            }//for x
          }//if
        }//for cur2

        long double *mij2 = new long double[card[j]];
        for(int y=0; y<card[j]; y++){
          mij2[y]=0;
          for(int x=0; x<card[i]; x++)
            mij2[y]+=mij[x][y];
        }//for y

        for(int ct=0; ct<card[i]; ct++)
          delete mij[ct];
        delete mij;

        long double denom = 0;
        for(int y=0; y<card[j]; y++){
          denom+=mij2[y];
        }
        if(denom==0){
          fprintf(status,"denom=0, i = %d, j = %d\n",i,j);
          fflush(status);
          printf("denom=0, i = %d, j = %d\n",i,j);
          exit(EXIT_FAILURE);
        }
        for(int y=0; y<card[j]; y++){
          mij2[y]=mij2[y]/denom;
        }//for y

        if((fabs(cur->m[0]-mij2[0])>0.0001)||(fabs(cur->m[1]-mij2[1])>0.0001)){
          if((fabs(cur->m[0]-mij2[0])>maxdiff)||(fabs(cur->m[1]-mij2[1])>maxdiff)){
            if((fabs(cur->m[0]-mij2[0]))>(fabs(cur->m[1]-mij2[1])))
              maxdiff = fabs(cur->m[0]-mij2[0]);
            else
              maxdiff = fabs(cur->m[1]-mij2[1]);
          }//if((fabs...maxdiff))
        }//if((fabs...0.0001)

        Node *cur2=&neighbors[j];
        while(cur2->id!=i)
          cur2=cur2->next;
        for(int y=0; y<card[j]; y++){
          cur->m[y] = mij2[y];
          if(cur->m[y]!=cur->m[y]){
            fprintf(status,"cur->m[y] is NaN\n");
            fprintf(status,"i = %d, j = %d, y = %d\n",i,j,y);
            printf("i = %d, j = %d, y = %d\n",i,j,y);
            exit(EXIT_FAILURE);
          }
          cur2->m[(y+card[i])] = mij2[y];
        }//for y
        delete mij2;
      }//if(!evidpos[j]&&!evidneg[j])
    }//for cur
  }//sendmsgsout

  void calcsendorder(){
    int numout=0; //number that have sent msgs out
    int nodesperprot = numterms*2;
    for(int i=0; i<numnodes; i++){
      int protind = i/nodesperprot; //integer division
      if(evidpos[i]&&(sdist[protind]<=threshold)){
        sendorder[numout] = i;
        for(Node *cur=&neighbors[i];cur!=NULL;cur=cur->next){
          int j=cur->id;
          if(j==-1){
	    printf("tried to enqueue -1 - make sure end of line after last line of goedges, other files\n");
	    exit(EXIT_FAILURE);
	  }//if(j==-1)
	  
          if(!visited[j]&&!jcibool[j]){
            jci.enqueue(j);
            jcibool[j] = true;
          }//if(!visited...)
        }//for cur
        visited[i] = true;
        numout++;
      }//if(evidpos[i])
    }//for i
    for(int i=0; i<numnodes; i++){
      int protind = i/nodesperprot; //integer division
      if(evidneg[i]&&(sdist[protind]<=threshold)){
        sendorder[numout] = i;
        for(Node *cur=&neighbors[i];cur!=NULL;cur=cur->next){
          int j=cur->id;
          if(j==-1){
            printf("tried to enqueue -1 - make sure end of line after last line of goedges, other files\n");
            exit(EXIT_FAILURE);
          }//if(j==-1)

          if(!visited[j]&&!jcibool[j]){
            jci.enqueue(j);
            jcibool[j] = true;
          }//if(!visited...)
        }//for cur
        visited[i] = true;
        numout++;    
      }//if(evidneg[i])
    }//for i
    while(!jci.empty()){
      int injci=jci.dequeue();
      if(!visited[injci]){
        sendorder[numout]=injci;
        for(Node *cur=&neighbors[injci];cur!=NULL;cur=cur->next){
          int j=cur->id;
          if(!visited[j]&&!jcibool[j]){
            jci.enqueue(j);
            jcibool[j] = true;
          }//if(!visited...)
        }//for cur
        visited[injci] = true;
        numout++;
      }//if(!visited[injci])
    }//while(!jci.empty())
    delete visited;
    delete jcibool;
    bool *senddefined = new bool[numnodes];
    for(int i=0; i<numnodes; i++)
      senddefined[i] = false;
    int i = 0;
    while(sendorder[i]>=0){
      senddefined[sendorder[i]] = true;
      i++;
    }//while
    numreachable = i;
    for(int i=0; i<numnodes; i++){
      if(!senddefined[i]){
	sendorder[numreachable] = i;
	senddefined[i] = true;
        numreachable++;
      }//if
    }//for i
    delete senddefined;
    printf("numnodes = %d\n",numnodes);
    printf("numreachable = %d\n",numreachable);
  }//calcsendorder

  void calcsendorder2(){
    int nodesperprot = numterms*2;
    int ct = 0;
    for(int i=0; i<numnodes; i++){
      int protind = i/nodesperprot; //integer division
      if(sdist[protind]<=threshold){
        sendorder[ct] = i; //will go to unreachables if no threshold-
                           //avoids doing a dfs or bfs
        ct++;
      }
    }
    numreachable = ct;
  }//calcsendorder2

  void propagate(FILE *status){
    fprintf(status,"numreachable = %d\n",numreachable);
    fflush(status);
    int t=0; //t=# iterations
    bool finished=false;
    while(!finished){
      fprintf(status,"iteration %d - starting down\n",t);
      fflush(status);
      maxdiff=-1;   
      for(int i=0;i<numreachable;i++){
        sendmsgsout(sendorder[i],status);//down
      }//for i
      //writemsgs(); //will write to same file each time, overwrites
      if(maxdiff==-1){ //ie if converged
        fprintf(status,"converged on down - iteration %d\n",t);
        fprintf(status,"maxdiff = %Lf\n",maxdiff);
        fflush(status);
        finished=true;
      }//if(maxdiff==-1)
      else{ //only go up if not converged on down
        fprintf(status,"maxdiff = %Lf\n",maxdiff);
        fprintf(status,"iteration %d - starting 2nd down\n",t);
        fflush(status);
        maxdiff=-1;
	//for(int i=numreachable-1;i>=0;i--) //up
	for(int i=0;i<numreachable;i++) //down again
          sendmsgsout(sendorder[i],status);
      
        //writemsgs();
    
        t++;
        if(maxdiff==-1)
          finished=true;
        fprintf(status,"maxdiff = %Lf\n",maxdiff);
        fflush(status);
      }//else
    }//while(!finished)
    fprintf(status,"num iterations (down and up) = %d\n",t);
    fflush(status);
  }//propagate

  void writesendorder(){
    printf("in writesendorder()\n");    
    FILE *send = fopen("sendorder.txt","w");
    for(int i=0; i<numnodes; i++)
      fprintf(send,"%d\n",sendorder[i]);
    fflush(send);
    fclose(send);
  }//writesendorder

  void writemsgs(){
    FILE *msgs = fopen("msgs.txt","w"); //will overwrite each time
    for(int i=0; i<numnodes; i++){
      for(Node *cur=&neighbors[i];cur!=NULL;cur=cur->next){
        int j = cur->id;
        for(int y=0; y<card[j]; y++)
          fprintf(msgs,"%Lf\n",cur->m[y]);  
      }//for cur
    }//for i
    fflush(msgs);
    fclose(msgs);
  }//writemsgs

  void readmsgs(FILE *status){ //for resume 
    FILE *msgs = fopen("msgs.txt","r"); //read in same order as wrote
    long double msgin;
    for(int i=0; i<numnodes; i++){
      for(Node *cur=&neighbors[i];cur!=NULL;cur=cur->next){
        int j = cur->id;
        for(int y=0; y<card[j]; y++){
          fscanf(msgs,"%Lf",&msgin);
          cur->m[y]=msgin;
        }//for y
      }//for cur
    }//for i
    //fflush(msgs);  
    fclose(msgs);
  }//readmsgs

  void readsendorder(){
    FILE *sendo = fopen("sendorder.txt","r"); //read in same order as wrote
    int so;
    for(int i=0; i<numnodes; i++){
      fscanf(sendo,"%d",&so);
      sendorder[i]=so;
    }//for i
    fclose(sendo);
  }//readsendorder

  void computebeliefs(FILE *status){
    for(int i=0;i<numnodes;i++){
      int nodesperprot = numterms*2;
      int protind = i/nodesperprot; //integer division
      if((sdist[protind]<=threshold)&&(i%2==1)&&(!evidpos[i])&&(!evidneg[i])){ //only compute if within threhsold, unannotated term node
        for(int j=0;j<2;j++){
          long double b=phi[i][j];
          for(Node *cur=&neighbors[i];cur!=NULL;cur=cur->next){
            Node *cur2=&neighbors[cur->id];
            while(cur2->id!=i)
              cur2=cur2->next;
            if(j==0)
              b=b*cur2->m[0];
            else //here j=1
              b=b*cur2->m[1];
          }//for cur
          belief[i][j]=b;
        }//for j
        double denom=belief[i][0]+belief[i][1];
        belief[i][0]=belief[i][0]/denom;
        belief[i][1]=belief[i][1]/denom;
      }//if
    }//for i
  }//computebeliefs()
        
  void printresults(char *outfn){
    FILE *out = fopen(outfn,"w");
    int nodesperprot = numterms*2;
    for(int i=0;i<numnodes;i++){
      int protind = i/nodesperprot; //integer division  
      if((sdist[protind]<=threshold)&&(i%2==1)&&(!evidpos[i])&&(!evidneg[i])){ //only prints for non-evidential term nodes in thresh
        int termind = (i%nodesperprot-1)/2;
        fprintf(out,"%d,%d,%Lf\n",protind,termind,belief[i][1]);
      }//if(!contains(evid,i))
    }//for i
    delete(neighbors);
    delete(sendorder);
    delete(phi);
    delete(belief);
    delete(card);
    delete(numpars);
    delete(parents);
    if(fclose(out)!=0){
      printf("error closing out\n");
      exit(EXIT_FAILURE);  
    }
    fclose(out);
  }//printresults()

}; //Graph
