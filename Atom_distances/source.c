#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#define DIMENSIONS 3

struct atom{
	float x,y,z;
	char type[4];
	char name[4];
	int id;
	int naa;
};

struct Sphere{
	struct atom center;
	float radious;
};
typedef struct Node{
	struct atom split_atom;
	int level;
	int id;
	struct Node *left;
	struct Node *right;
} KDNode;

struct Pset{
	struct atom *atoms;
	int size;
};
struct KeyValuePair {
    char key[4];
    int value;
};
struct KeyValuePair *include;

KDNode *BuildKDTree(struct Pset ,int level);
int compare_x(const void* a, const void* b);
int compare_y(const void* a, const void* b);
int compare_z(const void* a, const void* b);
float distanceSquare(struct atom,struct atom);
void search_InSphere(const struct Sphere sphere, const KDNode *root, struct Pset pts, int *i);
float distance(struct atom a1,struct atom a2);
struct Sphere search(const KDNode *root, float *theresshold, int *naa);
KDNode *bruteForceSearch(KDNode *root , int target) ;

int main()
{
	char filename[30];
	printf("Enter the filename: ");
    scanf("%s", filename); //1hfr.pdb
    FILE *fp=fopen(filename, "r"); 
    int data_size,aa1,aa2;
    char ch;
    int i=0;
    while(!feof(fp)){
    	fscanf(fp,"ATOM %d",&data_size);
    	ch=getc(fp);
	}
	rewind(fp);
	struct Pset pts;
	pts.size=data_size;
	pts.atoms = (struct atom*)malloc(pts.size * sizeof(struct atom));
	include = (struct KeyValuePair *)malloc(data_size * sizeof(struct KeyValuePair));
	while(!feof(fp)){
		if(ch=='\n' && fscanf(fp,"ATOM %d %s %s %*s %d %f %f %f ",&pts.atoms[i].id,&pts.atoms[i].type,&pts.atoms[i].name,&pts.atoms[i].naa,&pts.atoms[i].x,&pts.atoms[i].y,&pts.atoms[i].z)){
            include[i].value=0;
            strcpy(include[i].key, pts.atoms[i].type);
			i++;
		}
		ch=fgetc(fp);
	}
	fclose(fp);
    int max_naa=pts.atoms[i-2].naa;
	int level=0;
    KDNode *root=BuildKDTree(pts,level);
    free(pts.atoms);
    int a_1,ind;;
	float theresshold;
    
    double cpu_time_used;
    while(1){
		printf("Enter the number of  amino acid  (or -1 to quit): ");
		scanf("%d",&a_1);
		while(a_1>186) {
		    printf("Enter an aminoacid id between 1 and %d, or else -1 to quit\n",max_naa);
			scanf("%d",&a_1);
		}
		if(a_1==-1) break;
        KDNode *v;
        v=bruteForceSearch(root, a_1) ;
		printf("Enter the thresshold: ");
		scanf("%f",&theresshold);
		struct Pset pset;
	    pset.size=data_size;
		pset.atoms = (struct atom*)malloc(pset.size * sizeof(struct atom));
    	struct Sphere sphere;
    	sphere.center.x=v->split_atom.x;
		sphere.center.y=v->split_atom.y;
		sphere.center.z=v->split_atom.z;
	    sphere.radious=theresshold;
   		int iindex=0;
   		clock_t start, end;
        //double cpu_time_used;
        start = clock();
    	search_InSphere(sphere, root, pset, &iindex);
    	end = clock();
        cpu_time_used += ((double) (end - start)) / CLOCKS_PER_SEC;

    	int l;
    	for(l=0;l<iindex;l++){
    		printf("DISTANCE:%f between %s of %d aminoacid and %s of %d aminoacid \n",
			sqrt(distanceSquare(v->split_atom,pset.atoms[l])),pset.atoms[l].type,
			pset.atoms[l].naa,v->split_atom.type,v->split_atom.naa);
		}
		
		free(pset.atoms);
		for (ind=0;ind<data_size;ind++)
		    include[ind].value=0;
		    
		printf("Execution time: %f seconds\n", cpu_time_used);
	}
    free(include);
    return 0;
}

KDNode *BuildKDTree(struct Pset pts ,int level){
    KDNode *v=(KDNode *)malloc(sizeof(KDNode));
    const int axis = level % DIMENSIONS;
    if(axis== 0 && pts.size!=1){
    	qsort(pts.atoms, pts.size, sizeof(pts.atoms[0]), compare_x);
	}
	else if(axis== 1 && pts.size!=1) {
		qsort(pts.atoms, pts.size, sizeof(pts.atoms[0]), compare_y);
	}
	else if(axis== 2 && pts.size!=1){
		qsort(pts.atoms, pts.size, sizeof(pts.atoms[0]), compare_z);
	}
	int i_median = pts.size/ 2 + pts.size % 2-1;
	v->level = level;
    v->split_atom = pts.atoms[i_median];
    v->id=pts.atoms[i_median].id;
	if(pts.size<=1){
		v->left=NULL;
	    v->right=NULL;
		return v;
	}    
	else {
		int i;
		struct Pset P1,P2;
		P1.size=pts.size/ 2+pts.size%2;
		P2.size=pts.size-P1.size;
		P1.atoms = (struct atom*)malloc((P1.size) * sizeof(struct atom));
		P2.atoms = (struct atom*)malloc((P2.size) * sizeof(struct atom));
		for (i=0;i<P1.size;i++){
			P1.atoms[i]=pts.atoms[i];
		}
		for (i=0;i<P2.size;i++){
			P2.atoms[i]=pts.atoms[i+P1.size];
		}
		v->left=BuildKDTree(P1,level+1);
       	v->right=BuildKDTree(P2,level+1);
       	free(P1.atoms);
       	free(P2.atoms);
	    return v;
	}
}

int compare_x(const void* a, const void* b) {
	float fa = (*(struct atom*)a).x;
    float fb = (*(struct atom*)b).x;
    return (fa > fb) - (fa < fb);
}
int compare_y(const void* a, const void* b) {
	float fa = (*(struct atom*)a).y;
    float fb = (*(struct atom*)b).y;
	return ((*(struct atom*)a).y-(*(struct atom*)b).y);
}
int compare_z(const void* a, const void* b) {
	float fa = (*(struct atom*)a).z;
    float fb = (*(struct atom*)b).z;
	return ((*(struct atom*)a).z-(*(struct atom*)b).z);
}

void search_InSphere(const struct Sphere sphere, const KDNode *root, struct Pset pts, int *i)
{
    if (!root) return;
    const float d = distanceSquare(sphere.center,root->split_atom);
    const int axis = root->level % DIMENSIONS;
    const float d_split = ((axis == 0) ? root->split_atom.x - sphere.center.x : (axis == 1) ? root->split_atom.y - sphere.center.y : root->split_atom.z - sphere.center.z);
	const bool right_of_split = d_split <= 0;
    if (d < sphere.radious * sphere.radious)
    {
        if (strcmp(include[root->id-1].key, root->split_atom.type) == 0 && !include[root->id-1].value) {
        include[root->id-1].value = 1;  
        pts.atoms[*i]=(root->split_atom);
        (*i)++;
    	}
    }
    search_InSphere(sphere, right_of_split ? root->right : root->left, pts , i);
    if (d_split * d_split >= sphere.radious * sphere.radious) return;
    search_InSphere(sphere, right_of_split ? root->left : root->right, pts , i);
}

float distanceSquare(struct atom Sc,struct atom p){
	float distance;
	distance=(Sc.x-p.x)*(Sc.x-p.x)+(Sc.y-p.y)*(Sc.y-p.y)+(Sc.z-p.z)*(Sc.z-p.z);
	return distance;
}
float distance(struct atom a1,struct atom a2){
	float distance;
	distance=(a1.x-a2.x)*(a1.x-a2.x)+(a1.y-a2.y)*(a1.y-a2.y)+(a1.z-a2.z)*(a1.z-a2.z);
	distance=sqrt(distance);
	return distance;
}

KDNode *bruteForceSearch(KDNode *root, int target) {
    if (root == NULL || (root->split_atom.naa == target && strcmp(root->split_atom.type,"CA")==0)) {
        return root;
    }
    KDNode * leftResult = bruteForceSearch(root->left, target);
    KDNode * rightResult = bruteForceSearch(root->right, target);

    return (leftResult != NULL) ? leftResult : rightResult;
}
