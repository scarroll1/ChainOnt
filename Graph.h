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
  Node *neighbors; //adjacency list
  LinkedQueue jci;
  int *sendorder;
  long double **phi; //2-D array - must be double for function nodes/cond probs
  long double **belief;
  long double maxdiff;  
  int *card; //array of cardinalities of nodes
  int *numpars; //array of number of parents - only used for function nodes, ie even #'s
  int **parents; //will only be used for function nodes
  int numterms; //number of GO nodes, i.e. terms
  int numprots; //numprots - numnodes (above) is numterms*2*numprots
  long double *sdist; //shortest distance to any annotated protein
  long double threshold; //distance above which protein not included
  bool *visited;
  bool *jcibool;
  bool *evidpos;
  bool *evidneg;
  int numreachable;
  long double maxsim;
  
  //function members
  //constructor
  Graph(int numnodes,bool bayes){
    //printf("Graph const cp1 - numnodes = %d\n",numnodes);
    maxsim = 1; 
    this->numnodes=numnodes;
    visited = new bool[numnodes];
    jcibool = new bool[numnodes];
    evidpos = new bool[numnodes];
    evidneg = new bool[numnodes];
    numpars = new int[numnodes];
    parents = new int*[numnodes]; //won't all be used
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
        phi[2*i+1] = new long double[2]; //term nodes binary - function nodes' phi declared after constructor call
        belief[2*i+1] = new long double[2]; //belief not necessary for function nodes
        numpars[2*i] = 0; //start assuming no parents, fill in with parents later
      }//for i
    }//else

    LinkedQueue jci;
    card = new int[numnodes];

    //belief only needed for GO nodes, not "function nodes" introduced when 
    //converting from BN to MN - as a result, the entries of belief 
    //corresponding to function nodes will not be used/calculated
    
  }//constructor

  void add(int prot,Node *n){
    if(neighbors[prot].id==-1){
      neighbors[prot]=*n;
    }//if
    else{
      Node *cur=&neighbors[prot];
      while(cur->next!=NULL)
        cur=cur->next;
      cur->next=n;
    }//else
  }//add

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
      else{ //uniform
        phi[check][0]=1;
        phi[check][1]=1;
      }//else
    }//for i
  }//definerestofphi
  
  void sendmsgsout(int i,FILE *status){ //i in Yedidia's msg passing
    for(Node *cur=&neighbors[i];cur!=NULL;cur=cur->next){
      int j=cur->id; //j in Yedidia's msg passing
      if(!evidpos[j]&&!evidneg[j]){
        long double **mij = new long double*[card[i]];
        for(int ct=0; ct<card[i]; ct++)
          mij[ct] = new long double[card[j]];
        if((i%2==1)&&(j%2==1)){
          long double mult;
          for(int x=0; x<card[i]; x++){
            for(int y=0; y<card[j]; y++){
              if(x==y)
                mult=cur->s; //will be 0/1 for function nodes
              else
                mult=maxsim-cur->s;
              mij[x][y]=phi[i][x]*mult;
            }//for y
          }//for x
        }//if
        else{ //exactly one of i and j is even, ie function node
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
            } //if(j...)
          for(int x=0; x<card[i]; x++) //card[i] = 2^(numparents+1)
            for(int y=0; y<card[j]; y++){ //card[j] = 2
			                  //or vice-versa
              int limit; //relevant position in binary number
              if(parentpos==-1) //child, highest order binary digit
                limit = numpars[funid]; //child is highest order bin
                          //digit, which corresponds to 2^numparents
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
        }//else
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
          cur2=cur2->next; //have to change messages from i to j here now as well - will still save time

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

  void calcsendorder(){ //if no convergence, try interspersing b/t
                        //function and term nodes (odds and evens)
                        //calcsendorder2 does this in a very
                        //simple way
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
            printf("i = %d\n",i);
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
            printf("i = %d\n",i);
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
    
    for(int i=0; i<numnodes; i++){
      senddefined[i] = false;
    }

    int i = 0;
    while((i<numnodes)&&(sendorder[i]>=0)){
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
    
    //printf("about to call writesendorder()\n");    
    //writesendorder(); //for checkpointing
  }//calcsendorder

  void calcsendorder2(){
    int nodesperprot = numterms*2;
    int ct = 0;
    for(int i=0; i<numnodes; i++){
      int protind = i/nodesperprot; //integer division
      if(sdist[protind]<=threshold){ //only prints for non-evidential term nodes
        sendorder[ct] = i; //will go to unreachables if no theshold-
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

      //writemsgs(); //also for checkpointing
                     //will write to same file each time, overwrites
       
      if(maxdiff==-1){ //ie if converged
        fprintf(status,"converged on down - iteration %d\n",t);
        fprintf(status,"maxdiff = %Lf\n",maxdiff);
        fflush(status);
        finished=true;
      }//if(maxdiff==-1)
      else{ //only go up if not converged on down
        fprintf(status,"maxdiff = %Lf\n",maxdiff);
        //fprintf(status,"iteration %d - starting up\n",t);
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
    printf("num iterations (down and up) = %d\n",t);
    fprintf(status,"num iterations (down and up) = %d\n",t);
      //says down and up, may be down-down
    fflush(status);
  }//propagate

  void writesendorder(){
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

  void readmsgs(FILE *status){ //for resume if crash, etc.
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
      if((sdist[protind]<=threshold)&&(i%2==1)&&(!evidpos[i])&&(!evidneg[i])){ //only compute if within threhsold, term node, not evidential
        for(int j=0;j<2;j++){
          long double b=phi[i][j];
          for(Node *cur=&neighbors[i];cur!=NULL;cur=cur->next){
            Node *cur2=&neighbors[cur->id];
            while(cur2->id!=i)
              cur2=cur2->next;
            if(j==0)
              b=b*cur2->m[0];
            else //j==1
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
      if((sdist[protind]<=threshold)&&(i%2==1)&&(!evidpos[i])&&(!evidneg[i])){ //only prints for non-evidential term nodes within thresh
        int termind = (i%nodesperprot-1)/2;
        fprintf(out,"%d,%d,%Lf\n",protind,termind,belief[i][1]);
      }//if
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
  }//printresults()

  void computemarg(){ //computes all joints, then marginalizes
		      //can run this on a sufficiently small graph
		      //compare results with BP
    //assume 16 nodes
    //printf("numnodes = %d\n",numnodes);
    int numjoints=1;
    for(int i=0; i<numnodes; i++)
      numjoints = numjoints*card[i];
    //printf("numjoints = %d\n",numjoints);
    double *jointprob= new double[numjoints];
    int **dig = new int*[numjoints];
    for(int i=0; i<numjoints; i++)
      dig[i] = new int[numnodes];

    int jointind=0;
    for(int i15=0;i15<card[15];i15++)
    for(int i14=0;i14<card[14];i14++)
    for(int i13=0;i13<card[13];i13++)
    for(int i12=0;i12<card[12];i12++)
    for(int i11=0;i11<card[11];i11++)
    for(int i10=0;i10<card[10];i10++)
    for(int i9=0;i9<card[9];i9++)
    for(int i8=0;i8<card[8];i8++)
    for(int i7=0;i7<card[7];i7++)
    for(int i6=0;i6<card[6];i6++)
    for(int i5=0;i5<card[5];i5++)
    for(int i4=0;i4<card[4];i4++)
    for(int i3=0;i3<card[3];i3++)
    for(int i2=0;i2<card[2];i2++)
    for(int i1=0;i1<card[1];i1++)
    for(int i0=0;i0<card[0];i0++){
      //printf("jointind = %d\n",jointind);
      if(jointind==numjoints)
        printf("problem - jointind==numjoints\n");
      dig[jointind][0]=i0;
      dig[jointind][1]=i1;
      dig[jointind][2]=i2;
      dig[jointind][3]=i3;
      dig[jointind][4]=i4;
      dig[jointind][5]=i5;
      dig[jointind][6]=i6;
      dig[jointind][7]=i7;
      dig[jointind][8]=i8;
      dig[jointind][9]=i9;
      dig[jointind][10]=i10;
      dig[jointind][11]=i11;
      dig[jointind][12]=i12;
      dig[jointind][13]=i13;
      dig[jointind][14]=i14;
      dig[jointind][15]=i15;
      jointind++;     
    }//for i0
    if(numjoints!=jointind)
      printf("numjoints!=jointind\n");

    printf("dig[%d][%d] = %d\n",(numjoints-1),(numnodes-1),dig[numjoints-1][numnodes-1]);


    double *joint = new double[numjoints];
    double denom=0.0; //for normalization
    for(int jointind=0; jointind<numjoints; jointind++){
      //printf("jointind = %d\n",jointind);
      double psiprod = 1.0;
      bool breaknode=false;
      for(int node=0; node<numnodes; node++){
        for(Node *cur=&neighbors[node]; cur!=NULL; cur=cur->next){
          if(node<(cur->id)){ //ensures you only count each edge once
            int dig1 = dig[jointind][node]; //node 0 corresponds to 0th "digit"
            int dig2 = dig[jointind][cur->id]; //not binary digits
            bool cond2 = ((node%2==1)&&((cur->id)%2==1));
            if(cond2){
              if(dig1==dig2)
                psiprod = psiprod*cur->s;
              else
                psiprod = psiprod*(1-cur->s);
            }//if(cond2)
            else{ //here one is even, other odd
              int termid;
              int funid;
              if(node%2==1){
                termid=node;
                funid=cur->id;
              }
              else{
                termid=cur->id;
                funid=node;
              }
              int parentpos=-1;
              for(int i=0; i<numpars[funid]; i++){
                if(parents[funid][i]==termid){
                  parentpos=i;
                  break;
                }
              }//for i
              int limit; //relevant position in binary number
              if(parentpos==-1) //child, highest order binary digit
                limit = numpars[funid];
              else  
                limit = parentpos;
              int exp = 1; //will be 2^a
              for(int ct=0; ct<limit; ct++)
                exp=exp*2;
              int counter = dig[jointind][funid];//x 0..card[i]
              int bindig = dig[jointind][termid];//y 0..card[j]
                
              if((counter/exp)%2==bindig) //extracts the bin digit 
                psiprod = psiprod*1; 
              else{
                psiprod = 0; //really = psiprod*0 - not used
                             //just for readability
                breaknode=true;
                break; //breaks out of for cur loop
              }//else

            }//else - here one is even, other odd
              
          }//if(node<cur->id)
        }//for cur
        if(breaknode)
          break; //breaks out of for node loop
      }//for node
      double phiprod=1.0;
      if(breaknode)
        joint[jointind]=0;
      else{ //here we know psiprod!=0
        for(int node=0; node<numnodes; node++){
          if(evidpos[node]||evidneg[node]||(node%2)==0) //function nodes are even, evidential, not in evid
            phiprod=phiprod*phi[node][dig[jointind][node]];
            if(phiprod==0)
              break; //breaks out of node loop
        }//for node
      }//else - psiprod!=0
      joint[jointind]=psiprod*phiprod;
      denom+=joint[jointind];
    }//for jointind
    //normalize
    for(int jointind=0; jointind<numjoints; jointind++)
      joint[jointind]=joint[jointind]/denom;

    //NOW MARGINALIZE
    double **marg = new double*[numnodes];
    for(int node=0; node<numnodes; node++){
      marg[node] = new double[2]; //will only need for binary term nodes
      if(!evidpos[node]&&!evidneg[node]&&(node%2==1)){ //o/w no need to do anything - even nodes are functions, not printing
        for(int m=1; m<2; m++){ //just doing for +, - is just 1-that
          marg[node][m]=0;
          for(int k=0; k<numjoints; k++){
            if(dig[k][node]==m)
              marg[node][m]+=joint[k];
          }//for k
          printf("%d,%d,%lf\n",node,m,marg[node][m]);
        }//for m  
      }//if(!contains...)
    }//for node
  
    delete(dig);
    delete(marg);

    delete(neighbors);
    delete(sendorder);
    delete(phi);
    delete(belief);
    delete(card);
    delete(numpars);
    delete(parents);
    
    printf("end computemarg\n");
  }//computemarg

  void dijkstrapath(FILE *status,char *outfn,int source,bool edges){
    long double infin = std::numeric_limits<long double>::max();
    int_hash_set fringe;
    long double *dist = new long double[numnodes];
    bool *done = new bool[numnodes];
    int *parent = new int[numnodes];
    for(int i=0; i<numnodes; i++){ //step 1
      dist[i] = infin; //infinity
      done[i] = false;
      parent[i] = -1;
    }//for i
    dist[source]=0; //step 2
    done[source] = true;
    fprintf(status,"%d\n",source);
    fflush(status);
    for(Node *cur=&neighbors[source];cur!=NULL;cur=cur->next){ //step 3
      if(edges)
        dist[cur->id] = 1;
      else
        dist[cur->id] = 1-cur->s;
      if((cur->id)!=-1){	
        fringe.insert(cur->id);
        parent[cur->id] = source;
      }//if
    }//for cur   
    while(!fringe.empty()){ //step 4
      long double mindist;
      mindist = infin;
      int vm = -1; //nonsensical value
      int_hash_set::iterator it=fringe.begin();
      int_hash_set::iterator vmptr = it;
      while(it!=fringe.end()){
        int infringe=*it;
        if(dist[infringe]<mindist){  
          mindist = dist[infringe];
          vm = infringe;
          vmptr = it; //for quicker erase
        }//if
        it++;
      }//while(it!=evid.end())
      fringe.erase(vmptr); //erases vm from fringe - step 4a
      done[vm] = true;
      fprintf(status,"%d\n",vm);
      fflush(status);
      for(Node *cur=&neighbors[vm];cur!=NULL;cur=cur->next){ //step 4b
        if(edges){
          int w = cur->id;
          if(!done[w]){
            if(dist[w]==infin){ //ie if unseen
              dist[w] = dist[vm] + 1;
              fringe.insert(w);
              parent[w] = vm;
            }//if
            else if(contains(fringe,w)){
              long double newdist = dist[vm] + 1;
              if(newdist<dist[w]){
                dist[w] = newdist;
                parent[w] = vm;
              }//if
            }//else if
          }//if(!done)
        }//if(edges)
        else{ 
          int w = cur->id;
          if(!done[w]){
            if(dist[w]==infin){ //ie if unseen
              dist[w] = 1 - ((1-dist[vm])*cur->s);
              fringe.insert(w);
              parent[w] = vm;
            }//if
            else if(contains(fringe,w)){
              long double newdist = 1 - ((1-dist[vm])*cur->s);
              if(newdist<dist[w]){
                dist[w] = newdist;
                parent[w] = vm;  
              }//if
            }//else if  
          }//if !done
        }//else
      }//for cur
    }//while(!fringe.empty())
    //now write results to outfile
    FILE *dijk = fopen(outfn,"w");
      for(int i=0; i<numnodes; i++){
        fprintf(dijk,"%d %Lf",i,dist[i]);
        if((dist[i]<infin)&&(i!=source)){ //ie if reachable from source
          int onpath = parent[i];
          while(onpath!=source){
            fprintf(dijk," %d",onpath);
            onpath = parent[onpath];
          }//while
          fprintf(dijk," %d\n",source);
        }//if(dist[i]<2)
        else{
          fprintf(dijk,"\n");
        }//else
      }//for i
     
    if(fclose(dijk)!=0){
      printf("error closing dijk\n");
      exit(EXIT_FAILURE);
    } 
  }//dijkstrapath

  void conn_comp(FILE *status,char *outfn,FILE *ann){
    printf("inside conn_comp\n");
    FILE *outc = fopen(outfn,"w");
    int annprot;
    while(fscanf(ann,"%d",&annprot)==1)
      evidpos[annprot] = true; //might actually be a negatively ann.
                               //prot, but doesn't matter here
    for(int i=0; i<numnodes; i++){
      if(evidpos[i]&&!visited[i]){
        printf("about to call dfs on %d\n",i);
        fprintf(outc,"CONNECTED COMPONENT (ANNS):\n");
	fflush(outc);
        dfs(i,outc);
      }//if
    }//for i
  }//conn_comp 

  void dfs(int source, FILE *out){
    if(evidpos[source]) //again, from conn_comp, may be a negative
                        //prot, but doesn't matter for conn_comp
      //if you want to print out all proteins in each connected
      //component, regardless of whether annotated, comment out above
      //conditional
      fprintf(out,"%d\n",source);
    visited[source] = true;
    for(Node *cur=&neighbors[source];cur!=NULL;cur=cur->next){
      if((cur->id!=-1)&&!visited[cur->id])
        dfs(cur->id,out);
    }//for cur
  }//dfs  
   
}; //Graph
