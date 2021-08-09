// KDtree functions for use with Delayed Acceptance algorithm
//
// pvroot = pointer to root node defining a tree on R^d with maximum of
//          maxnleaf leaves on any leaf node.
//
// px     = x-coord in R^d.
// assoc  = associated scalar value (the log-likelihood).
// 
// pKD_AddElement returns a pointer to information that allows the user to
//                insert another estimate of the assoc value at exactly the
//                same x-coord, with merging by
//                new=log(n*exp(old)/(n+1)+exp(assoc)/(n+1))
// KD_MergeWithVeryNearest inserts the assoc at the nearest point using
//                the above formula
// 
//
// k      = number of nearest neighbours to use must have 2k <= maxnleaf
// wmax   = maximum weight (weight is propto 1/r up to this maximum)


void *pKD_CreateRootNode(int d, int maxnleaf); // returns pvroot
void KD_DestroyTree(void *pvroot);
void KD_PrintTree(void *pvroot); // prints out the whole tree
void KD_TreeDepthDiag(void *pvroot);  // diagnostic: looks at the depths of the different leaf nodes
void *pKD_AddElement(void *pvroot, double *px, double assoc); // since element: vector position and scalar (log-like) 
void *pKD_CreateAndFillFromInitialVec(int d, int maxnleaf, int n,double *px, double *passoc); // array of horizontal vectors and vector of associated scalars

void *pKD_CreateNearList(int d, int k); // returns pvnear
void KD_DestroyNearList(void *pvnear);
void KD_PrintNearList(void *pvnear);
// method=1 more efficient than method=0 as 0 starts from top of tree
// whereas method=0 starts from leaf node that contains the nearest point 
void KD_FindNearest(void *pvroot, void *pvnear, double *pxstar, char method); // find k nearest to xstar
// average with inverse distance weighting up to a maximum weight of wmax
double KD_MeanNearest(void *pvnear,double wmax);
int KD_MergeWithVeryNearest(void *pvnear,double assoc);
double KD_DistToVeryNearest(void *pvnear); // so you can see if it's less than epsilon and decide whether or not to merge 

