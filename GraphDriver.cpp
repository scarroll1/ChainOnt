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
#include <limits>
#include "Qnode.h"
#include "LinkedQueue.h"
#include "Node.h"
#include "Graph.h"

using namespace std;

Graph readGraph(char *fn1){
    
  FILE *f1=fopen(fn1,"r"); //graph.txt
  int numnodes;
  fscanf(f1,"%d",&numnodes);  
  Graph g(numnodes,false); //first # is # proteins - bool is bayes,
                           //false b/c bayesian GO network not needed
			   //here - see readGraph2    
  //evid not needed, only using this readGraph for dijkstrapath
  //or conn_comp

  //read in edges
  int numedgesread = 0;
  int prot1; int prot2; long double s;
  while(fscanf(f1,"%d %d %Lf",&prot1,&prot2,&s)==3){
    numedgesread++;
    long double *msgs = new long double[2];
    msgs[0]=1.0; msgs[1]=1.0; 
    Node *n1ptr = new Node(2,msgs,s,prot2); //1st # is # messages from prot1 to prot2
    g.add(prot1,n1ptr);
    Node *n2ptr = new Node(2,msgs,s,prot1); //1st # is # messages from prot2 to prot1
    g.add(prot2,n2ptr);  
  }//while
  //printf("numedgesread = %d\n",numedgesread);

  //no need to define phi, only using readgraph for dijkstrapath
  //or conn_comp

  //g.definephi();

  for(int i=0; i<numnodes; i++)
    g.card[i]=2;
    
  if(fclose(f1)!=0){
    printf("error closing f1\n");
    exit(EXIT_FAILURE);
  }
    
  //printf("exiting readGraph\n");

  return g;
    
}//readGraph

Graph readGraph2(char *fn0,char *fn1,char *fn2,char *fn3,long double thresh,FILE *status){

  //fprintf(status,"readGraph2 cp1\n");
  fflush(status);
  FILE *f0 = fopen(fn0,"r"); //pgraph.txt
  int numgonodes;
  FILE *f1 = fopen(fn1,"r"); //gograph.txt
  fscanf(f1,"%d",&numgonodes);
  if(fclose(f1)!=0){ //will reopen later as an ifstream to read
                     //parent-child relationships and cond'l probs
    printf("error closing f1\n");
    exit(EXIT_FAILURE);
  }
  int numprts;
  fscanf(f0,"%d",&numprts);
  fprintf(status,"numprts = %d\n",numprts);
  fprintf(status,"numgonodes = %d\n",numgonodes);
  fflush(status);
  
  Graph g(2*numgonodes*numprts,true); //first # is # nodes in BN, and
                //MN will have one function node above each BN node,
		//so exactly twice as many nodes - bool is bayes,
		//implying replicate bayesian networks
  g.numprots = numprts;
  g.numterms = numgonodes;  

  long double *sd=new long double[numprts]; //will hold shortest
			   //distance to any annotated protein
  if(thresh==-1){ //should always coincide with fn3 being ""
    g.threshold = std::numeric_limits<long double>::max(); 
    //above and below equivalent to not imposing a threshold distance
    for(int i=0; i<numprts; i++)
      sd[i] = -1;
  }
  else{
    g.threshold = thresh;
    FILE *f3 = fopen(fn3,"r"); //dijshortest.txt
    int orfid; //orf = "open reading frame"
    long double shortest;
    while(fscanf(f3,"%d %Lf",&orfid,&shortest)==2)
      sd[orfid] = shortest;
  }
  g.sdist = sd; //otherwise undefined  

  char line[65536]; //will have to change size if too many 
		    //annotations in annotations.txt
  ifstream f2(fn2); //annotations.txt
  f2.getline(line,65536);
  int linecount = 0;
  while(!f2.eof()){
    linecount++;
    char **num = new char*[1000]; //assumes no more than 999
				  //annotations for any protein
    for(int i=0; i<1000; i++)
      num[i] = new char[5]; //biggest possible number is biggest
		            //term index, 5 digits in GO
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
        else //if(temp!="NOT") 
          num[numnumsread] = temp; //will get segmentation fault if > 999 annotations (including repeats)
    }//while(temp!=NULL)
    int prot = atoi(num[0]);
    for(int i=1; i<numnumsread; i++){
      fprintf(status,"%s\n",&num[i]);
      fflush(status);
      if((strcmp(num[i-1],"NOT")!=0)&&(strcmp(num[i],"NOT")!=0))  //2nd condition is in case i=1, num[1]="NOT"
        g.evidpos[prot*numgonodes*2+2*atoi(num[i])+1] = true;
      else if(strcmp(num[i-1],"NOT")==0)
        g.evidneg[prot*numgonodes*2+2*atoi(num[i])+1] = true;
    }
    delete(num);
    f2.getline(line,65536);        
  }//while(!f2.eof())
  f2.close();
  //fprintf(status,"readGraph2 cp3\n");
  //fflush(status);
  fprintf(status,"about to create go edges for all proteins\n");
  fflush(status);
  
  int numcondprobsread = 0;
  int *child = new int[1000];
  int *numparents = new int[1000];
  
  //now create numprots GO BN's
  char **num = new char*[1000];
  for(int i=0; i<1000; i++)
    num[i] = new char[18];
  int *numsread = new int[1000];
  int **parenttwodim = new int*[1000];
  for(int i=0; i<1000; i++)
    parenttwodim[i] = new int[1000];
  long double *condprob = new long double[1000];
  ifstream f1p(fn1);
	   //p stands for prime, cannot use f1 again
  f1p.getline(line,65536); //reusing variable line
  f1p.getline(line,65536); //first line already read, number of gonodes
  while(!f1p.eof()){
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
    f1p.getline(line,65536);
    numcondprobsread++;
  }//while
  f1p.close();
  int numinthresh = 0;
  //note - if not imposing a threshold, can make g.threshold somewhere
  //between max possible distance and long
  //double's maximum allowable value, which represents infinity - then
  //unnecessary edges won't be created here, saving space - in the
  //case of similarity-based distance, any value over 1.0 will suffice
  for(int prot=0; prot<numprts; prot++){
    if(sd[prot]<=g.threshold)
      numinthresh++;  
    for(int j=0; j<numcondprobsread; j++){
      int *par=new int[numparents[j]];
      int parindex=0;
      for(int i=numparents[j];i>=1;i--){
        par[parindex]=prot*numgonodes*2+2*parenttwodim[j][i-1]+1; //rightmost parent corresponds to 0th bin digit
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
    
      //define phi for function node - needs to be of right cardinality
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
  }//for prot

/*  
  printf("thresh = %Lf\n",thresh);
  printf("numinthresh = %d\n",numinthresh);
  fprintf(status,"thresh = %Lf\n",thresh);
  fprintf(status,"numinthresh = %d\n",numinthresh);
*/

  delete num;

  //fprintf(status,"readGraph2 cp5\n");
  //fflush(status);

  int prot1;
  int prot2;
  long double s;
  int numedgesread=0;
  //now add edges b/t proteins GO nodes based on pedges.txt
  while(fscanf(f0,"%d %d %Lf",&prot1,&prot2,&s)==3){
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
        Node *n1ptr = new Node(4,msgs12,s,node2); //1st # is # messages from prot1 to prot2
        g.add(node1,n1ptr);
        Node *n2ptr = new Node(4,msgs21,s,node1); //1st # is # messages from prot2 to prot1
        g.add(node2,n2ptr);
      }//for(gonode)
    }//if both within threshold distance create edge, otherwise not  
  }//while
  
  fprintf(status,"numedgesread = %d\n",numedgesread);
  //fprintf(status,"readGraph2 cp6\n");
  //fflush(status);

  if(fclose(f0)!=0){
    printf("error closing f0\n");
    exit(EXIT_FAILURE);
  }
  
  g.definerestofphi();
  
  fprintf(status,"readGraph2 done\n");
  fflush(status);

  return g;
}//readGraph2
  
int main(int argc, char **argv){
  Graph g(0,false); //dummy graph, just to make compiler happy
  char *outfn; //output filename
  FILE *status = fopen(argv[argc-1],"w");
  fprintf(status,"argc = %d\n",argc);
  fflush(status);
  if((strcmp(argv[1],"bel_prop")==0)||(strcmp(argv[1],"bel_prop_neigh")==0)){
    //each protein is represented by a GO BN
    if(strcmp(argv[1],"bel_prop")==0){
      printf("bel_prop\n");
      if(argc!=7){ 
        printf("in bel_prop mode, argc must be 7");
        fprintf(status,"in bel_prop mode, argc must be 7");
        exit(EXIT_FAILURE);
      }    
      printf("argc is 7\n");
                    //argv[2] = pgraph.txt - has numproteins and
      //specifies edges / similarities
                    //argv[3] = gograph.txt - has number of terms in
      //GO and has parents and cond prob's
                    //argv[4] = annotations.txt - has annotations - first number each line is protein number; following numbers are GO indices, in orignal order - includes all ancestors
                    //argv[5] = gdout.txt - output file
		    //argv[6] = status.txt

      g=readGraph2(argv[2],argv[3],argv[4],"",-1,status);
			// "",-1 indicates not in neighborhood mode
    }
    else if(strcmp(argv[1],"bel_prop_neigh")==0){
      printf("bel_prop_neigh\n");
      if(argc!=9){
        printf("in bel_prop_neigh mode, argc must be 9");
        fprintf(status,"in bel_prop_neigh mode, argc must be 9");
        exit(EXIT_FAILURE);
      }    
      printf("argc is 9\n");
                    //argv[2] = pgraph.txt - has numproteins and
      //specifies edges / similarities
                    //argv[3] = gograph.txt - has number of terms in
      //GO and has parents and cond prob's
                    //argv[4] = annotations.txt - has annotations - first number each line is protein number; following numbers are GO indices, in orignal order - includes all ancestors
                    //argv[5] = gdout.txt - output file
		    //argv[6] = dijsh.txt - contains one line for each
      //protein, first number on line is protein number, 2nd is distance
      //to closest annotated protein
		    //argv[7] = threshold distance - only perform
      //belief propagation within this distance of annotated proteins
		    //argv[8] = status.txt

      g=readGraph2(argv[2],argv[3],argv[4],argv[6],atof(argv[7]),status);
    }

    outfn = argv[5];
    fprintf(status,"about to call g.calcsendorder(), g.propagate(status)\n");
    fflush(status);
    //g.calcsendorder(); //down or down-up
    g.calcsendorder2(); //simple schedule, interspersed 
    fprintf(status,"done with calcsendorder\n");    
    fflush(status);
    g.propagate(status);
    fprintf(status,"done with propagate\n");
    fflush(status);
    g.computebeliefs(status);
    fprintf(status,"about to call printresults with outfn\n");
    fflush(status);
    g.printresults(outfn); //include out later
  }
  else if(strcmp(argv[1],"conn_comp")==0){
    printf("conn_comp\n");
    // ./a.out conn_comp pgraph.txt annotations.txt conn_comp.txt status.txt
    if(argc!=6){
      printf("in conn_comp mode, argc must be 6");
      fprintf(status,"in conn_comp mode, argc must be 6");
      exit(EXIT_FAILURE);
    }    
    printf("argc==6\n");
    g=readGraph(argv[2]); //pgraph.txt
    outfn = argv[4];
    FILE *ann=fopen(argv[3],"r");
    printf("about to call conn_comp\n");
    g.conn_comp(status,outfn,ann);
  }
  else if((strcmp(argv[1],"short_sim")==0)||(strcmp(argv[1],"short_edg")==0)){
    if(argc!=6){
      printf("in short modes, argc must be 6");
      fprintf(status,"in short modes, argc must be 6");
      exit(EXIT_FAILURE);
    }    
    printf("argc==6\n");
    g=readGraph(argv[2]); //pgraph.txt
    printf("done with readGraph\n");
    outfn = argv[3];
    //printf("done with readgraph\n");
    //printf("about to call dijkstrapath\n");
    bool numedges;
    if(strcmp(argv[1],"short_edg")==0){
      printf("short_edg\n");
      numedges = true;
    }
    else{ //using similarity based distance
      printf("short_sim\n");
      numedges = false;
    }
    g.dijkstrapath(status,outfn,atoi(argv[4]),numedges);
    printf("after call to dijkstrapath\n");    
  }
  fprintf(status,"exiting GraphDriver.main\n");
  fflush(status);  
  fclose(status);
  return EXIT_SUCCESS;
}//main

