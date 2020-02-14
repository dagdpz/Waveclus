//Copyright (c) 2011 ashelly.myopenid.com under <http://www.opensource.org/licenses/mit-license>

//Code from https://gist.github.com/ashelly/5665911
 
#include <stdlib.h>
#include <mex.h>
 
//Customize for your data Item type
typedef double Item;
#define ItemLess(a,b)  ((a)<(b))
#define ItemMean(a,b)  (((a)+(b))/2)
 
typedef struct Mediator_t
{
   Item* data;  //circular queue of values
   int*  pos;   //index into `heap` for each value
   int*  heap;  //max/median/min heap holding indexes into `data`.
   int   N;     //allocated size.
   int   idx;   //position in circular queue
   int   ct;    //count of items in queue
} Mediator;
 
/*--- Helper Functions ---*/
 
#define minCt(m) (((m)->ct-1)/2) //count of items in minheap
#define maxCt(m) (((m)->ct)/2)   //count of items in maxheap 
 
//returns 1 if heap[i] < heap[j]
int mmless(Mediator* m, int i, int j)
{
   return ItemLess(m->data[m->heap[i]],m->data[m->heap[j]]);
}
 
//swaps items i&j in heap, maintains indexes
int mmexchange(Mediator* m, int i, int j)
{
   int t = m->heap[i];
   m->heap[i]=m->heap[j];
   m->heap[j]=t;
   m->pos[m->heap[i]]=i;
   m->pos[m->heap[j]]=j;
   return 1;
}
 
//swaps items i&j if i<j;  returns true if swapped
int mmCmpExch(Mediator* m, int i, int j)
{
   return (mmless(m,i,j) && mmexchange(m,i,j));
}
 
//maintains minheap property for all items below i/2.
void minSortDown(Mediator* m, int i)
{
   for (; i <= minCt(m); i*=2)
   {  if (i>1 && i < minCt(m) && mmless(m, i+1, i)) { ++i; }
      if (!mmCmpExch(m,i,i/2)) { break; }
   }
}
 
//maintains maxheap property for all items below i/2. (negative indexes)
void maxSortDown(Mediator* m, int i)
{
   for (; i >= -maxCt(m); i*=2)
   {  if (i<-1 && i > -maxCt(m) && mmless(m, i, i-1)) { --i; }
      if (!mmCmpExch(m,i/2,i)) { break; }
   }
}
 
//maintains minheap property for all items above i, including median
//returns true if median changed
int minSortUp(Mediator* m, int i)
{
   while (i>0 && mmCmpExch(m,i,i/2)) i/=2;
   return (i==0);
}
 
//maintains maxheap property for all items above i, including median
//returns true if median changed
int maxSortUp(Mediator* m, int i)
{
   while (i<0 && mmCmpExch(m,i/2,i))  i/=2;
   return (i==0);
}
 
/*--- Public Interface ---*/
 
 
//creates new Mediator: to calculate `nItems` running median. 
//mallocs single block of memory, caller must free.
Mediator* MediatorNew(int nItems)
{
   int size = sizeof(Mediator)+nItems*(sizeof(Item)+sizeof(int)*2);
   Mediator* m=  malloc(size);
   m->data= (Item*)(m+1);
   m->pos = (int*) (m->data+nItems);
   m->heap = m->pos+nItems + (nItems/2); //points to middle of storage.
   m->N=nItems;
   m->ct = m->idx = 0;
   while (nItems--)  //set up initial heap fill pattern: median,max,min,max,...
   {  m->pos[nItems]= ((nItems+1)/2) * ((nItems&1)?-1:1);
      m->heap[m->pos[nItems]]=nItems;
   }
   return m;
}
 
 
//Inserts item, maintains median in O(lg nItems)
void MediatorInsert(Mediator* m, Item v)
{
   int isNew=(m->ct<m->N);
   int p = m->pos[m->idx];
   Item old = m->data[m->idx];
   m->data[m->idx]=v;
   m->idx = (m->idx+1) % m->N;
   m->ct+=isNew;
   if (p>0)         //new item is in minHeap
   {  if (!isNew && ItemLess(old,v)) { minSortDown(m,p*2);  }
      else if (minSortUp(m,p)) { maxSortDown(m,-1); }
   }
   else if (p<0)   //new item is in maxheap
   {  if (!isNew && ItemLess(v,old)) { maxSortDown(m,p*2); }
      else if (maxSortUp(m,p)) { minSortDown(m, 1); }
   }
   else            //new item is at median
   {  if (maxCt(m)) { maxSortDown(m,-1); }
      if (minCt(m)) { minSortDown(m, 1); }
   }
}
 
//returns median item (or average of 2 when item count is even)
Item MediatorMedian(Mediator* m)
{
   Item v= m->data[m->heap[0]];
   if ((m->ct&1)==0) { v= ItemMean(v,m->data[m->heap[-1]]); }
   return v;
}
 
 
/*--- Test Code ---*/
#include <stdio.h>
void PrintMaxHeap(Mediator* m)
{
   int i;
   if(maxCt(m))
      printf("Max: %3d",m->data[m->heap[-1]]);
   for (i=2;i<=maxCt(m);++i)
   {
      printf("|%3d ",m->data[m->heap[-i]]);
      if(++i<=maxCt(m)) printf("%3d",m->data[m->heap[-i]]);
   }
   printf("\n");
}
void PrintMinHeap(Mediator* m)
{
   int i;
   if(minCt(m))
      printf("Min: %3d",m->data[m->heap[1]]);
   for (i=2;i<=minCt(m);++i)
   {
      printf("|%3d ",m->data[m->heap[i]]);
      if(++i<=minCt(m)) printf("%3d",m->data[m->heap[i]]);
   }
   printf("\n");
}
 
void ShowTree(Mediator* m)
{
   PrintMaxHeap(m);
   printf("Mid: %3d\n",m->data[m->heap[0]]);
   PrintMinHeap(m);
   printf("\n");
}
 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for proper number of input/output arguments
    if (nrhs != 2) {
        mexErrMsgTxt("Two input arguments required.");
    } else if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    //Extract mxArrays pointers from input argument
    const mxArray *a;
    a=prhs[0];
    
    const mxArray *b;
    b=prhs[1];
    
    //Check validity of inputs
    if(mxGetM(a)!=1) {
        mexErrMsgTxt("Input a must be a row vector.");
    }

    if (!mxIsDouble(a) || mxIsComplex(a)) {
        mexErrMsgTxt("Input a must be array of doubles");
    }
    
    if ((mxGetM(b)!=1) || (mxGetN(b)!=1)) {
        mexErrMsgTxt("Input b must be a scalar.");
    }
    
    
    //Extract sizes
    size_t a_ncols = mxGetN(a);
    

    
    //Extract C arrays from input
    double *a_array = mxGetPr(a);
    double *filter_order=mxGetPr(b);
    
    // create the output matrix
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)a_ncols,mxREAL );
    
    //Get C array for output matrix
    double *out_array= mxGetPr(plhs[0]);
    
    
    int ind_start=0;
    
    Mediator* m = MediatorNew((int)*filter_order);
    //printf("%d\n",(int)*filter_order);
    
    for (int k=0;k<a_ncols;k++) 
    {
        MediatorInsert(m,a_array[k]);
        out_array[k]=MediatorMedian(m);
    }
    
}


// int main(int argc, char* argv[])
// {
//    int i,v;
//    Mediator* m = MediatorNew(15);
//  
//    for (i=0;i<30;i++)
//    {
//       v = rand()&127;
//       printf("Inserting %3d \n",v);
//       MediatorInsert(m,v);
//       v=MediatorMedian(m);
//       printf("Median = %3d.\n\n",v);
//       ShowTree(m);
//    }
// }