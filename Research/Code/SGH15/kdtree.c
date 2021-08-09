// See the main() function at the bottom for examples of use for each function

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_rng.h> // only needed for testing
#include <gsl/gsl_randist.h> // only needed for testing


typedef struct _node {
  char type;      // branch or leaf
  int d;          // number of dimensions
  int maxnleaf; // even number
  int dsplit;     // dimension on which (to) split when become branch
  double xsplit;  // value of x on which to split
  struct _node *pChildLo;
  struct _node *pChildHi;
  struct _node *pParent;
  int nleaf;    // number of leaf values
  double *pxleaves; // pointer to array of leaf xvecs
  double *passoc;    // pointer to vector of associated scalars
  int *pn; // vector of number of times this x (or close by) has occured
} Node;

typedef struct nearlist {
  int d;
  int k;              // # nearest neighbours to store
  double *pxstar;     // point of interest
  char current;     // either all A variables or all B variables  
  double *passocA; // # current list of nearest k assocs (nearest first)
  double *passocB; 
  int *pnA;        // # assoc is the mean of n
  int *pnB;     
  double *prsqA;   // # list of squared distances
  double *prsqB;  
  int ind_nearest; // index of nearest point in leaf
  Node *pspecpos_nearest; // pointer to leaf containing nearest point
} NearList;



static void printdvec(int n, double *pd) {
  int i;
  for (i=0;i<n;i++) {
    printf("%f ",pd[i]);
  }
  printf("\n");
}

static void Indent(int depth) {
  int i;
  for (i=0;i<depth;i++) {
    printf("  ");
  }
}
static void PrintNodes(Node *pNode,int depth) {
  if (pNode->type=='L') {
    int i,j;
    Indent(depth);
    printf("Leaf at depth %d with %d entries (split at %d)\n",depth,pNode->nleaf,pNode->maxnleaf);
    for (i=0;i<pNode->nleaf;i++) {
      Indent(depth);
      for (j=0;j<pNode->d;j++) {
	printf("%f  ",pNode->pxleaves[i*pNode->d+j]);
      }
      printf(": %f (%d)\n",pNode->passoc[i],pNode->pn[i]);
    }
  }
  else if (pNode->type=='B') {
    Indent(depth);
    printf("Branch (dim=%d) splitting at x_%d=%f\n",pNode->d,pNode->dsplit+1,pNode->xsplit);
    PrintNodes(pNode->pChildLo,depth+1);
    PrintNodes(pNode->pChildHi,depth+1);
  }
}

void KD_PrintTree(void *pvroot) {
  Node *proot=(Node*)pvroot;
  
  PrintNodes(proot,0);
}

void DepthNodes(Node *pnode, int depth, int *pdpthct) {
  if (pnode->type=='L') {
    pdpthct[depth] ++;
  }
  else if (pnode->type=='B') {
    DepthNodes(pnode->pChildLo,depth+1,pdpthct);
    DepthNodes(pnode->pChildHi,depth+1,pdpthct);
  }
}
void MaxDepth(Node *pnode, int depth, int *pmax_depth) {
  if (pnode->type=='L') {
    if (depth> *pmax_depth) {
      *pmax_depth=depth;
    }
  }
  else if (pnode->type=='B') {
    MaxDepth(pnode->pChildLo,depth+1,pmax_depth);
    MaxDepth(pnode->pChildHi,depth+1,pmax_depth);
  }
}

void KD_TreeDepthDiag(void *pvroot) {
  Node *pnode=(Node*)pvroot;
  int max_depth=0;
  int *pdpthct;
  int i;
  int totnleaf=0,totdepth=0;

  printf("=========================\n");
  printf("= Tree depth diagnostic =\n");
  printf("=========================\n");
  printf("Depth\tFreq.\n");

  MaxDepth(pnode,0,&max_depth);

  pdpthct=(int*)malloc((max_depth+1)*sizeof(int));
  for (i=0;i<(max_depth+1);i++) {
    pdpthct[i]=0;
  }

  DepthNodes(pnode,0,pdpthct);

  for (i=0;i<(max_depth+1);i++) {
    printf("%d\t%d\n",i,pdpthct[i]);
    totnleaf+=pdpthct[i];
    totdepth+= i*pdpthct[i];
  }
  printf("Mean depth %f\n",(double) totdepth / (double) totnleaf);

  free(pdpthct);
}

static void AllocForLeaf(Node *pNode) {
  pNode->pxleaves=(double*)malloc(pNode->maxnleaf*pNode->d*sizeof(double));
  pNode->passoc=(double*)malloc(pNode->maxnleaf*sizeof(double));
  pNode->pn=(int*)malloc(pNode->maxnleaf*sizeof(int));
}
static void FreeForLeaf(Node *pNode) {
  free(pNode->pxleaves);
  free(pNode->passoc);
  free(pNode->pn);
}


static void CommonNewLeafNode(Node *pNode) {
  pNode->type='L';
  pNode->xsplit=0;
  pNode->pChildLo=0;
  pNode->pChildHi=0;
};


void *pCreateRootNode(int d, int maxnleaf, char alloc) {
  Node *pNew= (Node*)malloc(sizeof(Node));
  assert(maxnleaf>1);
  pNew->d=d;
  pNew->maxnleaf=maxnleaf;
  pNew->dsplit=0;
  pNew->nleaf=0;
  pNew->pParent=0;

  CommonNewLeafNode(pNew);
  if (alloc) {
    AllocForLeaf(pNew);
  }

  return(pNew);
}
void *pKD_CreateRootNode(int d, int maxnleaf) {
  return (void*) pCreateRootNode(d,maxnleaf,1);
}

static Node *pCreateChildLeafNode(void *pvParent, char alloc) {
  Node *pParent = (Node*)pvParent;
  Node *pNew = (Node*)malloc(sizeof(Node));

  pNew->d=pParent->d;
  pNew->maxnleaf=pParent->maxnleaf;
  pNew->dsplit= (pParent->dsplit+1) % pParent->d;
  pNew->nleaf=0;
  pNew->pParent=pParent;

  CommonNewLeafNode(pNew);
  if (alloc) {
    AllocForLeaf(pNew);
  }

  return(pNew);
}


static void DestroyNode(Node *pthis) {
  if (pthis->type=='B') {
    assert(pthis->pChildLo!=0);
    DestroyNode(pthis->pChildLo);
    assert(pthis->pChildHi !=0);
    DestroyNode(pthis->pChildHi);
  }
  else if (pthis->type=='L') {
    assert(pthis->pxleaves!=0);
    assert(pthis->passoc!=0);
    assert(pthis->pn!=0);
    free(pthis->pxleaves);
    free(pthis->passoc);
    free(pthis->pn);
  }
  else {
    assert(1==2);
  }
}

void KD_DestroyTree(Node *pvroot) {
  Node *proot=(Node*)pvroot;
  DestroyNode(proot);
}

// Decide on hi/lo status of each vector
static int subtle_split(double *px, int sz, int d, int dsp, 
			char *philo, double *pxsplit) {
  size_t *p= (size_t*)malloc(sz*sizeof(size_t));
  // If sz even, nlo=nhi, else nhi=nlo-1.
  int nlo = (sz+1)/2, nhi = sz-nlo;
  int nlo_curr=0, i;

  gsl_sort_index(p,px+dsp,d,sz);
  // If sz even then average of nearest two else central value
  *pxsplit = 0.5*(px[dsp+p[nlo-1]*d]+px[dsp+p[nhi]*d]);
  for (i=0;i<sz;i++) {
    if (px[dsp+i*d]< *pxsplit) {
      philo[i]=0;
      nlo_curr++;
    }
    else if (px[dsp+i*d]> *pxsplit) {
      philo[i]=1;
    } 
    else {
      philo[i]=2; // exactly equal
    }
  }
  for (i=0;i<sz;i++) {
    if (philo[i]==2) {
      if (nlo_curr<nlo) {
	philo[i]=0;
	nlo_curr++;
      }
      else {
	philo[i]=1;
      }
    }
  }  
  free(p);

  return nlo;
} 


// Do not re-order when allocating to children as we always want the last entry 
// in a particular node to be the one most recently added
static Node *pSortAndSplit(Node *pParent, Node *pChildLo, Node *pChildHi) {
  int i,ilo=0,ihi=0, nlo;
  int d=pParent->d;
  int maxnleaf=pParent->maxnleaf;
  char *philo = (char*)malloc(maxnleaf*sizeof(char));
  Node *paddedto;
  int last;

  assert(pParent->nleaf==maxnleaf);
  nlo=subtle_split(pParent->pxleaves, maxnleaf, d, pParent->dsplit, philo, 
		   &pParent->xsplit);
  assert(nlo==maxnleaf/2);

  for (i=0;i<maxnleaf;i++) {
    if (philo[i]==0) {
      memcpy(pChildLo->pxleaves+ilo*d,pParent->pxleaves+i*d,d*sizeof(double));
      pChildLo->passoc[ilo]=pParent->passoc[i];
      pChildLo->pn[ilo]=pParent->pn[i];
      ilo++;
      last=0;
    }
    else {
      memcpy(pChildHi->pxleaves+ihi*d,pParent->pxleaves+i*d,d*sizeof(double));
      pChildHi->passoc[ihi]=pParent->passoc[i];
      pChildHi->pn[ihi]=pParent->pn[i];
      ihi++;
      last=1;
    }
  }
  pChildLo->nleaf=nlo;
  pChildHi->nleaf=maxnleaf-nlo;

  if (last == 0) {
    paddedto=pChildLo;
  }
  else {
    paddedto=pChildHi;
  }

  free(philo);

  return paddedto;
}

static void LeafToBranch(Node *pParent, char freeleaf) {
  pParent->type='B';
  pParent->nleaf=0;
  pParent->maxnleaf=0;
  if (freeleaf) {
    FreeForLeaf(pParent);
  }
}

static Node *pGiveBirth(Node *pParent) {
  Node *paddedto;
  assert(pParent->type=='L');
  assert(pParent->pChildLo==0);
  assert(pParent->pChildHi==0);
  assert(pParent->pxleaves!=0);
  assert(pParent->passoc!=0);
  assert(pParent->pn!=0);

  pParent->pChildLo=pCreateChildLeafNode(pParent,1);
  pParent->pChildHi=pCreateChildLeafNode(pParent,1);

  paddedto=pSortAndSplit(pParent,pParent->pChildLo,pParent->pChildHi);

  LeafToBranch(pParent,1);

  return paddedto;
}

static Node* pAddToTree(Node *pNode, int depth, double *px, double assoc) {
  Node *paddedto, *pchild;
  int dsp=pNode->dsplit;

  if (pNode->type=='B') {
    if (px[dsp]<pNode->xsplit) {
      pchild=pNode->pChildLo;
    }
    else if (px[dsp]>pNode->xsplit) {
      pchild=pNode->pChildHi;
    }
    else { 
      int k=rand();
      
      if ((k&1) == 0) {
	pchild=pNode->pChildLo;
      }
      else {
	pchild=pNode->pChildHi;
      }
    }
    paddedto=pAddToTree(pchild,depth+1,px,assoc);
  }
  else {
    assert(pNode->type=='L');
    assert(pNode->nleaf<pNode->maxnleaf);

    memcpy(pNode->pxleaves+pNode->nleaf*pNode->d,px,pNode->d*sizeof(double));
    pNode->passoc[pNode->nleaf]=assoc;
    pNode->pn[pNode->nleaf]=1;
    pNode->nleaf++;
    if (pNode->nleaf==pNode->maxnleaf) {
      paddedto=pGiveBirth(pNode);
    }
    else {
      paddedto=pNode;
    }

  }
  return paddedto;
}

void *pKD_AddElement(void *pvroot, double *px, double assoc) {
  Node *pNode = (Node*)pvroot;
  Node *paddedto=pAddToTree(pNode,0,px,assoc);
  return (void*)paddedto;
}


static double reaverage(double av, int n, double new) {
  double dn= (double)n;
  double newav = av+log(dn/(dn+1) + exp(new-av)/(dn+1));
  return newav;
}


void *pKD_CreateNearList(int d, int k) {
  NearList *pnear=(NearList*)malloc(sizeof(NearList));
  pnear->d=d;
  pnear->k=k;
  pnear->passocA=(double*)malloc(k*sizeof(double));
  pnear->pnA=(int*)malloc(k*sizeof(int));
  pnear->prsqA=(double*)malloc(k*sizeof(double));
  pnear->passocB=(double*)malloc(k*sizeof(double));
  pnear->pnB=(int*)malloc(k*sizeof(int));
  pnear->prsqB=(double*)malloc(k*sizeof(double));
  pnear->pxstar=(double*)malloc(d*sizeof(double));
  return (void*)pnear;
}

void KD_DestroyNearList(void *pvnear) {
  NearList *pnear=(NearList*)pvnear;
  free(pnear->passocA);
  free(pnear->pnA);
  free(pnear->prsqA);
  free(pnear->passocB);
  free(pnear->pnB);
  free(pnear->prsqB);
  free(pnear->pxstar);
  free(pnear);
}

static void NearSwap(NearList *pnear) {
  if (pnear->current=='A') {
    pnear->current='B';
  }
  else {
    pnear->current='A';
  }
}

static void CurrNear(NearList *pnear, double **pprsq, double **ppassoc, 
		     int **ppn) {
  if (pnear->current=='A') {
    *pprsq=pnear->prsqA;
    *ppassoc=pnear->passocA;
    *ppn=pnear->pnA;
  }
  else {
    *pprsq=pnear->prsqB;
    *ppassoc=pnear->passocB;
    *ppn=pnear->pnB;
  }
}
static void NextNear(NearList *pnear, double **pprsq, double **ppassoc, 
		     int **ppn) {
  if (pnear->current=='B') {
    *pprsq=pnear->prsqA;
    *ppassoc=pnear->passocA;
    *ppn=pnear->pnA;
  }
  else {
    *pprsq=pnear->prsqB;
    *ppassoc=pnear->passocB;
    *ppn=pnear->pnB;
  }
}

void KD_PrintNearList(void *pvnear) {
  NearList *pnear=(NearList*)pvnear;
  double *passoc,*prsq;
  int *pn,i;

  printf("NearList: k=%d, d=%d (%c)\n",pnear->k,pnear->d,pnear->current);
  printf("xstar=");
  printdvec(pnear->d,pnear->pxstar);

  printf("assoc (n):\n");
  CurrNear(pnear,&prsq,&passoc,&pn);
  for (i=0;i<pnear->k;i++) {
    printf("%f (%d); rsq=%f\n",passoc[i],pn[i],prsq[i]);
  }
  printf("Nearest is index %d on Node %p\n",pnear->ind_nearest,(void*)pnear->pspecpos_nearest);
}

// squared distance to all leaves on a leaf node
static void SqDist(Node *pNode, double *pxstar, double *psqdist) {
  double *px=pNode->pxleaves;
  double r2,r;
  int d=pNode->d;
  int nleaf=pNode->nleaf;
  int i,j;

  for (i=0;i<nleaf;i++) {
    r2=0;
    for (j=0;j<d;j++) {
      r=pxstar[j]-px[i*d+j];
      r2+= r*r;
    }
    psqdist[i]=r2;
  }
}

static Node *pFillFromExactQuadrant(Node *pNode, NearList *pNear) {
  double *pxstar=pNear->pxstar;
  if (pNode->type=='B') { //  Branch
    if (pxstar[pNode->dsplit]<=pNode->xsplit) {
      pNode=pFillFromExactQuadrant(pNode->pChildLo,pNear);
    }
    else {
      pNode=pFillFromExactQuadrant(pNode->pChildHi,pNear);
    }
  }
  else {
    int nleaf=pNode->nleaf;
    size_t *p= (size_t*)malloc(nleaf*sizeof(size_t));
    double *psqdist = (double*)malloc(nleaf*sizeof(double));
    int ind,i;

    assert(pNode->type=='L');
    SqDist(pNode, pxstar, psqdist);
    gsl_sort_index(p,psqdist,1,nleaf);

    pNear->current='A';

    for (i=0;i<pNear->k;i++) {
      ind=p[i]; // i+1 th closest
      pNear->passocA[i]=pNode->passoc[ind];
      pNear->prsqA[i]=psqdist[ind];
      pNear->pnA[i]=pNode->pn[ind];
    }
    pNear->ind_nearest=p[0];
    pNear->pspecpos_nearest=pNode;
    //    PrintNodes(pNode,0);
  }

  return pNode;
}

// If any leaves on new node are nearer than furthest in nearlist
// then merge with the nearlist
static void MergeNearest(NearList *pNear, Node *pNode) {
  int nleaf=pNode->nleaf;
  size_t *p= (size_t*)malloc(nleaf*sizeof(size_t));
  double *psqdist = (double*)malloc(nleaf*sizeof(double));
  double *prsq_curr, *passoc_curr;
  int *pn_curr;

  assert(pNode->type=='L');
  SqDist(pNode, pNear->pxstar, psqdist);
  gsl_sort_index(p,psqdist,1,nleaf);

  CurrNear(pNear, &prsq_curr, &passoc_curr, &pn_curr);

  // Only interleave if nearest in new quadrant is nearer than
  // furthest current.

  if (psqdist[p[0]]<prsq_curr[pNear->k-1]) {
    double *prsq_next, *passoc_next;
    int i_next,i_curr=0,ip_thisQ=0, *pn_next;
    NextNear(pNear, &prsq_next, &passoc_next, &pn_next);

    for (i_next=0;i_next<pNear->k;i_next++) {
      int i_thisQ=p[ip_thisQ];
      if (prsq_curr[i_curr]<=psqdist[i_thisQ]) {
	prsq_next[i_next]=prsq_curr[i_curr];
	passoc_next[i_next]=passoc_curr[i_curr];
	pn_next[i_next]=pn_curr[i_curr];
	i_curr++;
      }
      else {
	prsq_next[i_next]=psqdist[i_thisQ];
	passoc_next[i_next]=pNode->passoc[i_thisQ];
	pn_next[i_next]=pNode->pn[i_thisQ];
	ip_thisQ++;
      }
    }

    // If very nearest has changed then update the
    // additional info on this.
    if (prsq_next[0]<prsq_curr[0]) {
      pNear->ind_nearest=p[0];
      pNear->pspecpos_nearest=pNode;
    }

    NearSwap(pNear);
  }

  free(p);
  free(psqdist);
}

static double maxradsq(NearList *pnear) {
  double maxradsq;

  if (pnear->current=='A') {
    maxradsq=pnear->prsqA[pnear->k-1];
  }
  else {
    maxradsq=pnear->prsqB[pnear->k-1];
  }

  return maxradsq;
}
static double closeradsq(int d, double *pclosest) {
  int i;
  double closestrsq=0;

  for (i=0;i<d;i++) {
    closestrsq += pclosest[i]*pclosest[i];
  }

  return closestrsq;
}

static Node *pChildNear(Node*pNode, NearList *pnear) {
  Node *pcn;
  if (pnear->pxstar[pNode->dsplit]<=pNode->xsplit) {
    pcn=pNode->pChildLo;
  }
  else {
    pcn=pNode->pChildHi;
  }
  return pcn;
}
static Node *pChildFar(Node*pNode, NearList *pnear) {
  Node *pcf;
  if (pnear->pxstar[pNode->dsplit]<=pNode->xsplit) {
    pcf=pNode->pChildHi;
  }
  else {
    pcf=pNode->pChildLo;
  }
  return pcf;
}

// Recursive search down from root-node for any closer than exact quadrant.
static void OtherQuadrants(Node *pNode, NearList *pNear,Node *pexact, double *pclosest) {
  if (pNode->type=='B') {
    Node *pcn=pChildNear(pNode,pNear), *pcf=pChildFar(pNode,pNear);
    double closersq=closeradsq(pNear->d,pclosest);
    double maxrsq=maxradsq(pNear);
    double nearsep=pclosest[pNode->dsplit];
    double farsep=pNode->xsplit-pNear->pxstar[pNode->dsplit];

    //    printf("maxrsq %f; closersq %f, farrsq %f\n",maxrsq,closersq,closersq-nearsep*nearsep+farsep*farsep);
    assert(closersq<maxrsq);

    // Always try the nearer child
    OtherQuadrants(pcn,pNear,pexact,pclosest);
    maxrsq=maxradsq(pNear);

    // Try the farther child if there is intersection
    if (closersq-nearsep*nearsep+farsep*farsep<maxrsq) {
      double *pnewclosest=(double *)malloc(pNear->d*sizeof(double));
      memcpy(pnewclosest,pclosest,pNear->d*sizeof(double));

      //      printf("newclosest b4: ");
      //printdvec(pNear->d,pnewclosest);

      pnewclosest[pNode->dsplit]=farsep;

      //      printf("newclosest after: ");
      //printdvec(pNear->d,pnewclosest);

      OtherQuadrants(pcf,pNear,pexact,pnewclosest);
      free(pnewclosest);
    }
  }
  else {
    if (pNode != pexact)  {
      MergeNearest(pNear,pNode);
    }
  }

}

// Recursive search up and down the node that would contain x*
static void OtherQuadrants2(Node *pNode, NearList *pNear, double *pclosest, char ascending) {
  //CSRM  printf("Node: %p\n",pNode);
  if (ascending) { // Must be a branch node
    Node *pcf=pChildFar(pNode,pNear);
    double maxrsq=maxradsq(pNear);
    double farsep=pNode->xsplit-pNear->pxstar[pNode->dsplit];
    //CSRM    printf("Ascended %f %f\n",maxrsq,farsep);
    // Never try the nearer child - we've just come from there!
    // Try the farther child if there is intersection
    if (farsep*farsep<maxrsq) {
      double *pnewclosest=(double *)malloc(pNear->d*sizeof(double));
      int i;
      // No distance except that induced by this far child
      for (i=0;i<pNear->d;i++) {
	pnewclosest[i]=0;
      }
      pnewclosest[pNode->dsplit]=farsep;

      //CSRM      printf("Descending to %p\n",pcf);
      OtherQuadrants2(pcf,pNear,pnewclosest,0); // descend far child
      free(pnewclosest);
    }
    if (pNode->pParent != 0) { // not the root
      // free(pclosest);
      OtherQuadrants2(pNode->pParent,pNear,0,1); // ascend
    }
  }
  else { // descending
    // CSRM    printf("Descended\n");
    if (pNode->type=='B') {
      Node *pcn=pChildNear(pNode,pNear), *pcf=pChildFar(pNode,pNear);
      double closersq=closeradsq(pNear->d,pclosest);
      double maxrsq=maxradsq(pNear);
      double nearsep=pclosest[pNode->dsplit];
      double farsep=pNode->xsplit-pNear->pxstar[pNode->dsplit];
      //CSRM      printf("Branch %f %f %f\n",maxrsq,nearsep,farsep);
      assert(closersq<maxrsq);

    // Always try the nearer child
      OtherQuadrants2(pcn,pNear,pclosest,0);
      maxrsq=maxradsq(pNear);

      // Try the farther child if there is intersection
      if (closersq-nearsep*nearsep+farsep*farsep<maxrsq) {
	double *pnewclosest=(double *)malloc(pNear->d*sizeof(double));
	memcpy(pnewclosest,pclosest,pNear->d*sizeof(double));

	pnewclosest[pNode->dsplit]=farsep;

	OtherQuadrants2(pcf,pNear,pnewclosest,0);
	free(pnewclosest);
      }
    }
    else { // It's a leaf, and it can never be the exact leaf
      MergeNearest(pNear,pNode);
    }
  } // end descending
}


void KD_FindNearest(void *pvroot, void *pvnear, double *pxstar, char method) {
  NearList *pnear=(NearList*)pvnear;
  Node *proot=(Node*)pvroot;
  Node *pexact;
  double *pclosest=(double *)malloc(pnear->d*sizeof(double));
  int i;

  memcpy(pnear->pxstar,pxstar,pnear->d*sizeof(double));
  pexact=pFillFromExactQuadrant(proot,pnear); //  Node that contains x^*

  if (method==0) {
    for (i=0;i<pnear->d;i++) {
      pclosest[i]=0; // closes possible point in top `quadrant' is at the point
    }
    OtherQuadrants(proot,pnear,pexact,pclosest); // descend from root
  }
  else {
    if (pexact->pParent != 0) { // pexact is not the root
      // CSRM      printf("pExact %p\n",pexact);
      OtherQuadrants2(pexact->pParent,pnear,0,1); // ascend/descend
    }
  }

  free(pclosest);
}

double KD_MeanNearest(void *pvnear,double wmax) {
  NearList *pnear=(NearList*)pvnear;
  double sum=0,sumw=0;
  double wmaxinv=1/wmax;
  double w,r;
  double *prsq,*passoc;
  int *pn;
  int i;

  CurrNear(pnear,&prsq,&passoc,&pn);

  for (i=0;i<pnear->k;i++) {
    r=sqrt(prsq[i]);
    if (r<wmaxinv) {
      w=wmax*(double)(pn[i]);;
    }
    else {
      w=1/r*(double)(pn[i]);
    }
    sumw += w;
    sum += passoc[i]*w;
  }

  return sum/sumw;
}

int KD_MergeWithVeryNearest(void *pvnear,double assoc) {
  NearList *pnear=(NearList*)pvnear;
  Node *pnode=pnear->pspecpos_nearest;
  int ind=pnear->ind_nearest;
  double n=pnode->pn[ind];

  pnode->pn[ind] = n+1;
  pnode->passoc[ind] = reaverage(pnode->passoc[ind],n,assoc);
  if (pnear->current=='A') {
    pnear->passocA[0]=pnode->passoc[ind];
    pnear->pnA[0]=pnode->pn[ind];
  }
  else {
    pnear->passocB[0]=pnode->passoc[ind];
    pnear->pnB[0]=pnode->pn[ind];
  }
  return(pnode->pn[ind]);
}

double KD_DistToVeryNearest(void *pvnear) {
  NearList *pnear=(NearList*)pvnear;
  double dnear2;
  if (pnear->current=='A') {
    dnear2=pnear->prsqA[0];
  }
  else {
    dnear2=pnear->prsqB[0];
  }
  return(sqrt(dnear2));
}

static void RecursiveBulkFill(Node *pNode, int sz, double *px, double *passoc, int *pn) {
  int d=pNode->d;

  if (sz>=pNode->maxnleaf) {
    char *philo= (char*)malloc(sz*sizeof(char));
    int szlo, szhi, ilo=0, ihi=0, i;
    double *pxlo, *pxhi, *passoclo, *passochi;
    int *pnlo, *pnhi;

    szlo=subtle_split(px, sz, d, pNode->dsplit, philo, &pNode->xsplit);
    szhi=sz-szlo; 

    pxlo=(double*)malloc(szlo*d*sizeof(double));
    pxhi=(double*)malloc(szhi*d*sizeof(double));
    passoclo=(double*)malloc(szlo*sizeof(double));
    passochi=(double*)malloc(szhi*sizeof(double));
    pnlo=(int*)malloc(szlo*sizeof(int));
    pnhi=(int*)malloc(szhi*sizeof(int));
      
    for (i=0;i<sz;i++) {
      if (philo[i]==0) {
	memcpy(pxlo+ilo*d,px+i*d,d*sizeof(double));
	passoclo[ilo]=passoc[i];
	pnlo[ilo]=pn[i];
	ilo++;
      }
      else {
	memcpy(pxhi+ihi*d,px+i*d,d*sizeof(double));
	passochi[ihi]=passoc[i];
	pnhi[ihi]=pn[i];
	ihi++;
      }
    }
    assert(ilo==szlo);
    assert(ihi==szhi);
    free(philo);

    pNode->pChildLo=pCreateChildLeafNode(pNode,0);
    pNode->pChildHi=pCreateChildLeafNode(pNode,0);
    LeafToBranch(pNode,0);

    RecursiveBulkFill(pNode->pChildHi, szhi, pxhi, passochi, pnhi);
    free(pxhi);
    free(passochi);
    free(pnhi);
    RecursiveBulkFill(pNode->pChildLo, szlo, pxlo, passoclo, pnlo);
    free(pxlo);
    free(passoclo);
    free(pnlo);
  }
  else {
    AllocForLeaf(pNode);
    pNode->nleaf=sz;
    memcpy(pNode->pxleaves,px,sz*d*sizeof(double));
    memcpy(pNode->passoc,passoc,sz*sizeof(double));
    memcpy(pNode->pn,pn,sz*sizeof(int));
  }
}

void *pKD_CreateAndFillFromInitialVec(int d, int maxnleaf, int n,double *px, double *passoc, int *pn) {
  if (n>0) {
    Node *proot = pCreateRootNode(d,maxnleaf,0);
    RecursiveBulkFill(proot, n, px, passoc, pn);
    return (void*) proot;
  }
  else {
    printf("Failed: vector size <= 0.");
    return (void*)1;
  }
}

// see kdtree.h for a description of all the user-callable functions
// below is some code used to test each of the functions
int main() {
  // d=dimension of the parameter
  // ksplit=2b so  must be even
  // knear= number of nearest neighbours. In this implementation
  // you must have knear <=b
  int d=3; int ksplit=10,knear=5; 
  double *pxstar=(double*)malloc(d*sizeof(double));
  int n=10000, n1=1000000, n2=1000, n3=3;
  double assoc=11, mn, mn1=0, mn2=0, dist;
  void *proot=pKD_CreateRootNode(d,ksplit);
  void *pnear=pKD_CreateNearList(d,knear);
  void *pposn;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  int i,j;

  gsl_rng_set(r,1234);

  // Create some vectors and associated scalars
  // and add them individually to the KD-tree
  for (i=0;i<n;i++) {
    assoc=0;
    for (j=0;j<d;j++) {
      double z=gsl_ran_gaussian(r,1.0);
      pxstar[j]=z;
      assoc += (j+1)*z;
    }
    assoc += gsl_ran_gaussian(r,1.0)/4;
    pposn=pKD_AddElement(proot,pxstar,assoc);
  }
  
  printf("Create done\n");

  // Create a large number of new vectors, and for each one, find the 
  // k nearest neighbours
  mn=0;
  for (i=0;i<n;i++) {
    void *pnearcmp;
    if ((i % 1000) == 0) {
      printf("%d\n",i);
    }
    for (j=0;j<d;j++) {
      pxstar[j]=gsl_ran_gaussian(r,1.0);
    }
    pnearcmp = pKD_CreateNearList(d,knear);
    KD_FindNearest(proot,pnearcmp,pxstar,1);
    // Every 10th iteration merge the scalar associated with
    // the nearest theta with a new N(0,1) scalar
    // think of this scalar as a new log-likelihood
    if ((i % 10) == 0 ) {
      KD_MergeWithVeryNearest(pnearcmp,gsl_ran_gaussian(r,1.0));
    }
    mn+=KD_MeanNearest(pnearcmp,1000);
    KD_DestroyNearList(pnearcmp);
  }
  printf("%f\n",mn);

  // Add values at a single point
  printf("Stage 1\n");  
  pxstar[0]=3;pxstar[1]=3;pxstar[2]=3;
  for (i=0;i<200;i++)  {
    pKD_AddElement(proot,pxstar,0.8);
  }
  KD_PrintTree(proot);
  KD_TreeDepthDiag(proot);

  // Add values all over the place
  printf("Stage 2\n");
  for (i=0;i<n1;i++) {
    assoc=0;
    for (j=0;j<d;j++) {
      double z=gsl_ran_gaussian(r,1.0);
      pxstar[j]=z;
      assoc += (j+1)*z;
    }
    //    sc += gsl_ran_gaussian(r,1.0)/4;
    pposn=pKD_AddElement(proot,pxstar,assoc);
  }
  KD_TreeDepthDiag(proot);

  // find nearest neighbours to each of a bunch of points
  printf("Stage3\n");
  mn1=0;
  mn2=0;
  for (i=0;i<n2;i++) {
    for (j=0;j<d;j++) {
      pxstar[j]=gsl_ran_gaussian(r,1.0);
    }
    // Method 0 was easier to code but is less computationally efficient
    KD_FindNearest(proot,pnear,pxstar,0);
    mn1+=KD_MeanNearest(pnear,100000000);
    KD_FindNearest(proot,pnear,pxstar,1);
    mn2+=KD_MeanNearest(pnear,100000000);
    //printf("Mean=%f\n",mn);
  }
  printf("Method means %f %f\n",mn1,mn2); // same answer
  KD_TreeDepthDiag(proot);

  // Find nearest neighbours and distance to the very nearest point
  // then merge a new scalar with the nearest point
  printf("Stage4\n");
  for (i=0;i<n3;i++) {
    int k;
    for (j=0;j<d;j++) {
      pxstar[j]=gsl_ran_gaussian(r,1.0);
    }
    for (k=1;k<5;k++) {
      KD_FindNearest(proot,pnear,pxstar,0);
      printf("Method 0\n");
      //	KD_PrintNearList(pnear);
      KD_FindNearest(proot,pnear,pxstar,1);
      printf("Method 1\n");
      //	KD_PrintNearList(pnear);
      dist=KD_DistToVeryNearest(pnear);
      printf("Dist to nearest = %f\n",dist);
      KD_MergeWithVeryNearest(pnear,17.0);
    }      
  }
  KD_TreeDepthDiag(proot);

  KD_DestroyNearList(pnear);
  KD_DestroyTree(proot);
  
  // Create a balanced tree from an initial vector
  {
    double pxall[600], passoc[200];
    int pn[200];
    void *proot;
    void *pnear=pKD_CreateNearList(3,10);
    for (i=0; i<200;i++) {
      passoc[i]=0;
      for (j=0;j<3;j++) {
	double z=gsl_ran_gaussian(r,1.0);
	pxall[i*3+j]=z;
	passoc[i] += (j+1)*z;
	pn[i]=1;
      }
    }
    proot=pKD_CreateAndFillFromInitialVec(3,20,200,pxall,passoc,pn);
    printf("Created\n");
    KD_PrintTree(proot);
    printf("Printed\n");
  
    for (i=0;i<10000;i++) {
      assoc=0;
      for (j=0;j<d;j++) {
	double z=gsl_ran_gaussian(r,1.0);
	pxstar[j]=z;
	assoc += (j+1)*z;
      }
      //    sc += gsl_ran_gaussian(r,1.0)/4;
      pposn=pKD_AddElement(proot,pxstar,assoc);
    }
    KD_PrintTree(proot);
    printf("tree printed\n");

    for (i=0;i<10;i++) {
      for (j=0;j<3;j++) {
	pxstar[j]=gsl_ran_gaussian(r,1.0);
      }
      
      KD_FindNearest(proot,pnear,pxstar,0);
      printf("Method 0\n");
      KD_PrintNearList(pnear);
      KD_FindNearest(proot,pnear,pxstar,1);
      printf("Method 1\n");
      KD_PrintNearList(pnear);
    }
    KD_TreeDepthDiag(proot);
    KD_DestroyTree(proot);
    KD_DestroyNearList(pnear);
  }

  return 0;
}
