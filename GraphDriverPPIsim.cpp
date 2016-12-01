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
#include <string>
#include <fstream>
#include "Qnode.h"
#include "LinkedQueue.h"
#include "Node.h"
#include "GraphPPIsim.h"

using namespace std;

Graph readGraph_ppisim(char *fn0,char *fn1,char *fn2,char *fn3,char
*fn4, char *fn5, long double thresh, FILE *status){
  FILE *f0 = fopen(fn0,"r"); //pgraph.txt
  FILE *f1 = fopen(fn1,"r"); //gridedges.txt
  FILE *f2 = fopen(fn2,"r"); //gograph.txt
  ifstream f3(fn3); //annotations.txt
  int numprts;
  fscanf(f0,"%d",&numprts);
  printf("numprts = %d\n",numprts);
  int numgonodes;
  fscanf(f2,"%d",&numgonodes);
  if(fclose(f2)!=0){ //will open stream later
    printf("error closing f2\n");
    exit(EXIT_FAILURE);
  }
  fprintf(status,"numprts = %d\n",numprts);
  fprintf(status,"numgonodes = %d\n",numgonodes);
  fflush(status);
  Graph g(2*numgonodes*numprts,true);
  g.numprots = numprts;
  g.numterms = numgonodes;  
  printf("numgonodes = %d\n",numgonodes);
  printf("numprts = %d\n",numprts);
  printf("g.numnodes = %d\n",g.numnodes);
  long double *sd=new long double[numprts];
  if(thresh==-1){ //should always coincide with fn3 being ""
    g.threshold = std::numeric_limits<long double>::max(); 
    //above and below equivalent to not imposing a threshold distance
    for(int i=0; i<numprts; i++)
      sd[i] = -1;
  }
  else{
    g.threshold = thresh;
    FILE *f5 = fopen(fn5,"r"); //dijshortest.txt
    int orfid;
    long double shortest;
    while(fscanf(f5,"%d %Lf",&orfid,&shortest)==2){
      sd[orfid] = shortest;
    }//while fscanf
  }//else
  g.sdist = sd;
  g.mult = new long double*[g.numterms];
  if(g.mult==NULL){
    printf("Out of memory - readGraph_ppisim, g.mult\n");
    exit(1);
  }
  for(int i=0; i<g.numterms; i++){
    g.mult[i] = new long double[4];
  }//for i
  FILE *f4 = fopen(fn4,"r"); //protedgpars.txt
  int termid;
  long double comp0;
  long double comp1;
  long double comp2;
  long double comp3;
  while(fscanf(f4,"%d %Lf %Lf %Lf %Lf",&termid,&comp0,&comp1,&comp2,&comp3)==5){
    g.mult[termid][0] = comp0;
    g.mult[termid][1] = comp1;
    g.mult[termid][2] = comp2;
    g.mult[termid][3] = comp3;
  }//while fscanf
  if(fclose(f4)!=0){
    printf("error closing f4\n");
    exit(EXIT_FAILURE);
  }
  char line[65536]; //will have to change size if many annotations in annotations.txt
  f3.getline(line,65536);
  int linecount = 0;
  while(!f3.eof()){
    linecount++;
    char **num = new char*[1000];
    for(int i=0; i<1000; i++)
      num[i] = new char[5]; //biggest possible number is biggest term index, 5 digits

    num[0]=strtok(line," ");
    int numnumsread=0;
    char *temp = "";
    while(temp!=NULL){ 
      numnumsread++;
      temp = strtok(NULL," ");
      if(temp!=NULL)
        if(numnumsread>999){
          printf("more than 999 annotations/NOTs, exiting\n");
          exit(EXIT_FAILURE);
        }//if
        else
          num[numnumsread] = temp;
    }//while(temp!=NULL)
    int prot = atoi(num[0]);
    for(int i=1; i<numnumsread; i++){
      //fprintf(status,"%s\n",&num[i]);
      //fflush(status);
      if((strcmp(num[i-1],"NOT")!=0)&&(strcmp(num[i],"NOT")!=0))  //2nd condition is in case i=1, num[1]="NOT"
        g.evidpos[prot*numgonodes*2+2*atoi(num[i])+1] = true;
      else if(strcmp(num[i-1],"NOT")==0){
        g.evidneg[prot*numgonodes*2+2*atoi(num[i])+1] = true;
      }//else if
    }
    delete(num);
    f3.getline(line,65536);        
  }//while(!f5.eof())
  f3.close();
  fprintf(status,"about to create go edges for all proteins\n");
  fflush(status);
  
  int numcondprobsread = 0;
  int *child = new int[1000];
  int *numparents = new int[1000];
  
  char **num = new char*[1000];
  for(int i=0; i<50; i++)
    num[i] = new char[18];
  int *numsread = new int[1000];
  int **parenttwodim = new int*[1000];
  for(int i=0; i<10; i++)
    parenttwodim[i] = new int[1000];
  long double *condprob = new long double[1000];
  printf("about to open f2p\n");
  ifstream f2p(fn2);
  f2p.getline(line,65536); //reusing variable line
  f2p.getline(line,65536); //first line already read above
  while(!f2p.eof()){
    int numsthisline = 0;
    num[numsthisline] = strtok(line," ");  
    char *temp = "";
    while(temp!=NULL){ 
      numsthisline++;
      temp = strtok(NULL," ");
      if(temp!=NULL){
        num[numsthisline] = temp;
      }//if(temp!=NULL)
    }//while(temp!=NULL)
    child[numcondprobsread] = atoi(num[0]);
    numparents[numcondprobsread] = numsthisline-2;
    numsread[numcondprobsread] = numsthisline;
    for(int j=1; j<=numsthisline-2; j++){
      parenttwodim[numcondprobsread][j-1] = atoi(num[j]);
    }//for j
    condprob[numcondprobsread] = atof(num[numsthisline-1]);
    f2p.getline(line,65536);
    printf("line = %s\n",line);
    numcondprobsread++;
  }//while
  f2p.close();
  for(int prot=0; prot<numprts; prot++){
    printf("prot - %d\n",prot);
    if(sd[prot]<=g.threshold){ //otherwise no need to create go edges  
      for(int j=0; j<numcondprobsread; j++){
        int *par=new int[numparents[j]];
	printf("cp\n");
        int parindex=0;
	printf("cp\n");
        for(int i=numparents[j];i>=1;i--){
	  printf("i = %d\n",i);
          par[parindex]=prot*numgonodes*2+2*parenttwodim[j][i-1]+1; 
                   //rightmost parent corresponds to 0th bin digit -
                   //should have highest index
                   //2*id+1 is mapping to new node #
          parindex++;
        }//for i
        int newchildid=prot*numgonodes*2+2*child[j]+1;
        int funid=prot*numgonodes*2+2*child[j];
        g.parents[funid] = par;
        g.numpars[funid] = numparents[j];
        int nummsgs2 = 2;
        for(int i=0; i<numparents[j]; i++)
          nummsgs2 = nummsgs2*2; //will be 2^(numparents+1)
        long double *msgs = new long double[2+nummsgs2];
        long double *msgs2 = new long double[2+nummsgs2];
        for(int i=0; i<(nummsgs2+2); i++){
          msgs[i]=1.0;
          msgs2[i]=1.0;
        }//for i
        Node *n1ptr = new Node((2+nummsgs2),msgs,condprob[j],newchildid); //2+nummsgs2 so you can have msgs both ways, saves time
        g.add(funid,n1ptr);
        Node *n2ptr = new Node((2+nummsgs2),msgs2,condprob[j],funid); //same as above
        g.add(newchildid,n2ptr);
        for(int i=0; i<numparents[j]; i++){
          msgs = new long double[2+nummsgs2];
          msgs2 = new long double[2+nummsgs2];
          for(int ct=0; ct<(nummsgs2+2); ct++){
            msgs[ct]=1.0;
            msgs2[ct]=1.0;
          }
          n1ptr = new Node((2+nummsgs2),msgs,condprob[j],par[i]);
          n2ptr = new Node((2+nummsgs2),msgs2,condprob[j],funid);
          g.add(funid,n1ptr);
          g.add(par[i],n2ptr);
        }//for i
        //make sure 2nd parent in input file has larger index than first parent, etc. monotone
        g.phi[funid] = new long double[nummsgs2];
        int halfway = nummsgs2/2; //nummsgs is cardinality of funid
        for(int i=0; i<halfway; i++){
          if(i==(halfway-1))
            g.phi[funid][i]=condprob[j]; //calculated from Beta dist, started with B(1,1)
          else
            g.phi[funid][i]=1;
        }//for i
        for(int i=halfway; i<nummsgs2; i++){
          if(i==(nummsgs2-1))
            g.phi[funid][i]=(1-condprob[j]);
          else
            g.phi[funid][i]=0;  
        }//for i
        g.card[newchildid]=2;
        g.card[funid]=nummsgs2;
      }//for j
    }//if sd[prot]
  }//for prot
  delete num;
  int prot1;
  int prot2;
  long double s;
  int numedgesread=0;
  while(fscanf(f1,"%d %d %Lf",&prot1,&prot2,&s)==3){
    if((sd[prot1]<=g.threshold)&&(sd[prot2]<=g.threshold)){
      numedgesread++;
      //now map prot1 and prot2 to all term nodes, link them
      for(int gonode=0; gonode<numgonodes; gonode++){
        long double *msgs12 = new long double[4];
        long double *msgs21 = new long double[4];
        for(int ct=0; ct<4; ct++){
          msgs12[ct] = 1.0;
          msgs21[ct] = 1.0;
        }
        int node1 = prot1*numgonodes*2+2*gonode+1;
        int node2 = prot2*numgonodes*2+2*gonode+1;
        Node *n1ptr = new Node(4,msgs12,0.5,node2,true);
        g.add(node1,n1ptr);
        Node *n2ptr = new Node(4,msgs21,0.5,node1,true);
        g.add(node2,n2ptr);
      }//for(gonode)
    }//if both within threshold distance create edge, otherwise not  
  }//while
  while(fscanf(f0,"%d %d %Lf",&prot1,&prot2,&s)==3){
    if((sd[prot1]<=g.threshold)&&(sd[prot2]<=g.threshold)){
      numedgesread++;
      for(int gonode=0; gonode<numgonodes; gonode++){
        int node1 = prot1*numgonodes*2+2*gonode+1;
        int node2 = prot2*numgonodes*2+2*gonode+1;
        Node *exists = g.findnode(node1,node2);
        if(exists!=NULL){
          exists->s = s; 
          exists = g.findnode(node2,node1);
          if(exists==NULL){
            printf("exists (node1, node2) problem\n");
            exit(1);
          }
          else
            exists->s = s; 
        }//if !NULL
        else{
          long double *msgs12 = new long double[4];
          long double *msgs21 = new long double[4];
          for(int ct=0; ct<4; ct++){
            msgs12[ct] = 1.0;
            msgs21[ct] = 1.0;
          }
          Node *n1ptr = new Node(4,msgs12,s,node2,false); 
	    //false represents no gridedge - not needed for no sim 
            //edge, just make sim's 0.5 by default, if no sim edge, 
            //will stay 0.5
          g.add(node1,n1ptr);
          Node *n2ptr = new Node(4,msgs21,s,node1,false);
          g.add(node2,n2ptr);
        }//else
      }//for(gonode)
    }//if
  }//while
  fprintf(status,"numedgesread = %d\n",numedgesread);
  fflush(status);
  if(fclose(f1)!=0){
    printf("error closing f1\n");
    exit(EXIT_FAILURE);
  }
  if(fclose(f0)!=0){
    printf("error closing f0\n");
    exit(EXIT_FAILURE);
  }
  g.definerestofphi(); //for evidential term nodes
  fprintf(status,"readGraph_ppisim done\n");
  fflush(status);
  return g;
}//readGraph_ppisim


int main(int argc, char **argv){
  FILE *status = fopen(argv[argc-1],"w");
  fprintf(status,"GD cp1\n");
  fprintf(status,"argc = %d\n",argc);
  fflush(status);
  Graph g(0,false); //dummy graph, just to make compiler happy
  fprintf(status,"GD cp1a\n");
  fflush(status);
  char *outfn; //output filename
  if(strcmp(argv[1],"bel_prop")==0){
    if(argc!=9){ 
      printf("in bel_prop mode, ppisim, argc must be 9");
      fprintf(status,"in bel_prop mode, ppisim, argc must be 9");
      exit(EXIT_FAILURE);
    }    
    g=readGraph_ppisim(argv[2],argv[3],argv[4],argv[5],argv[6],"",-1,status);
  }
  else if(strcmp(argv[1],"bel_prop_neigh")==0){
    if(argc!=11){ 
      printf("in bel_prop mode, ppisim, argc must be 11");
      fprintf(status,"in bel_prop mode, ppisim, argc must be 11");
      exit(EXIT_FAILURE);
    }    
    g=readGraph_ppisim(argv[2],argv[3],argv[4],argv[5],argv[6],argv[8],atof(argv[9]),status);
  }//else
  outfn = argv[7];
  printf("outfn = %s\n",outfn);
  fprintf(status,"about to call g.calcsendorder(), g.propagate(status)\n");
  fflush(status);
  //g.calcsendorder(); 
  g.calcsendorder2(); 
  g.propagate(status);
  g.computebeliefs(status);
  fprintf(status,"about to call printresults with outfn\n");
  fflush(status);
  g.printresults(outfn);
  fprintf(status,"exiting GraphDriver.main\n");
  fflush(status);  
  fclose(status);
  return EXIT_SUCCESS;
}//main

