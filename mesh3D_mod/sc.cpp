/*
Class to project face normals. 
Written by Benjamin Irving 20111102
*/

#include "sc.h"

//Contructor takes vertices and meshes as input
//initialisation lists required for passing previously uninitialised references
MeshProject::MeshProject(std::vector<double> & Fa1, std::vector<double> & Fb1, std::vector<double> & Fc1, std::vector<double> & Vx1, std::vector<double> & Vy1, std::vector<double> & Vz1): Fa(Fa1), Fb(Fb1), Fc(Fc1), Vx(Vx1), Vy(Vy1), Vz(Vz1) {
    
    FacesN=Fa.size();
    VerticesN=Vx.size();
    
    //Allocate memory for Face normals)
    FNx = new double [FacesN]; 
    FNy = new double [FacesN];
    FNz = new double [FacesN];
    
    //Allocating memory for Face centres
    // (Might be easier int future to use vectors).
    FCx = new double [FacesN];
    FCy = new double [FacesN];
    FCz = new double [FacesN];
    
   //Populating the vectors with -1 so can tell if just a place holder
   // These are the faces and vertices of each edge
    E_Face1.resize(2*FacesN, -1);
    E_Face2.resize(2*FacesN, -1);
    E_Vertex1.resize(2*FacesN, -1);
    E_Vertex2.resize(2*FacesN, -1);
}

//Deconstructor (called when class goes out of scope)
MeshProject::~MeshProject(){

// free allocated memory of Face Normals
   delete [] FNx;
   delete [] FNy;
   delete [] FNz;
   
   delete [] FCx;
   delete [] FCy;
   delete [] FCz;
   
}

int MeshProject::update_vertices(std::vector<double> & Vx1, std::vector<double> & Vy1, std::vector<double> & Vz1){
// coping new vector input to original vector
// we are actually copying the vector using 2 references but references behave
// exactly the same as variables.
//Luckily vectors can be copied directly using =

Vx=Vx1;
Vy=Vy1;
Vz=Vz1;

return 1;
}

/* Outputting the viewing fn */
double MeshProject::test1() {
    return ViewPos[2];
}

/* Another test function */
int MeshProject::another_test(){
    int xx=1;
    return xx;
}


/* FACE NORMALS
Calculate the normals for each face
*/
int MeshProject::normals_face(){
    
    // Counter
    int ii;
    int i0, i1, i2;
    int Er=1;
    
    //edges
    double e0x, e0y, e0z, e0l;
    double e1x, e1y, e1z, e1l;
    double e2x, e2y, e2z, e2l;
    
    for (ii=0; ii<FacesN; ii++) {

        // indices of the face vertices
        i0=(int)Fa.at(ii);
        i1=(int)Fb.at(ii);
        i2=(int)Fc.at(ii);
        
        
        // edge 1
        e0x=Vx.at(i0)-Vx.at(i1);
        e0y=Vy.at(i0)-Vy.at(i1);
        e0z=Vz.at(i0)-Vz.at(i1);
        
        // edge 2
        e1x=Vx.at(i1)-Vx.at(i2);
        e1y=Vy.at(i1)-Vy.at(i2);
        e1z=Vz.at(i1)-Vz.at(i2);
        
        // edge 3
        e2x=Vx.at(i2)-Vx.at(i0);
        e2y=Vy.at(i2)-Vy.at(i0);
        e2z=Vz.at(i2)-Vz.at(i0);
        
        // Normalise the vectors
        e0l = sqrt(e0x*e0x+e0y*e0y+e0z*e0z)+1e-14;
        e0x/=e0l; e0y/=e0l; e0z/=e0l;
        e1l = sqrt(e1x*e1x+e1y*e1y+e1z*e1z)+1e-14; 
        e1x/=e1l; e1y/=e1l; e1z/=e1l;
        e2l = sqrt(e2x*e2x+e2y*e2y+e2z*e2z)+1e-14; 
        e2x/=e2l; e2y/=e2l; e2z/=e2l;
        

        /* 
        Face normals
        Cross production and assuming vertice order is consistant with one direction
        */
        
        FNx[ii]=e0y*(e2z) - e0z * (e2y);
        FNy[ii]=e0z*(e2x) - e0x * (e2z);
        FNz[ii]=e0x*(e2y) - e0y * (e2x);
        
        /*
        Face centres
        */
        
        FCx[ii]=(Vx.at(i0) + Vx.at(i1) + Vx.at(i2))/3;
        FCy[ii]=(Vy.at(i0) + Vy.at(i1) + Vy.at(i2))/3;
        FCz[ii]=(Vz.at(i0) + Vz.at(i1) + Vz.at(i2))/3;
        
       // mexPrintf("FNx: %f, FNy: %f, FNz: %f \n", FNx[ii], FNy[ii], FNz[ii]);
        //mexPrintf("Square: %f \n", FNx[ii]*FNx[ii]+FNy[ii]*FNy[ii]+FNz[ii]*FNz[ii]);
        
    }
    return Er;
}

/*
UNIQUE EDGES (slow)
Create a list of unique edges and the related vertices, faces and face normals
*/
int MeshProject::unique_edges(){

    // Counter
    int ii, jj, kk;
    int i0, i1, i2;
    int Ecur [3][2];
    int found1;
    int Er=1;
    
    EList=0;
    
    //for each face
    for (ii=0; ii<FacesN; ii++) {
        
        // indices of the face vertices
        i0=(int)Fa.at(ii);
        i1=(int)Fb.at(ii);
        i2=(int)Fc.at(ii);
        
        Ecur[0][0]=i0; Ecur[0][1]=i1; 
        Ecur[1][0]=i1; Ecur[1][1]=i2;
        Ecur[2][0]=i2; Ecur[2][1]=i0;
        
        //for each edge of the face
        for (jj=0; jj<3; jj++) {
        
            found1=0;
        
            // for each edge that is already found
            // consider using stl::list to monitor values that are already found
            for (kk=0; kk<EList; kk++) {
                
                // check if the vertices of the edge match the current vertex
                if ((Ecur[jj][0]==E_Vertex1[kk] && Ecur[jj][1]==E_Vertex2[kk]) 
                    || (Ecur[jj][0]==E_Vertex2[kk] && Ecur[jj][1]==E_Vertex1[kk])) {
                
                    // add face to the matching edge
                    E_Face2[kk]=ii;
                    found1=1;
                
                }
                
            }
            if (found1==0) {
                
                // add face and vertexes to the end of the list
                E_Vertex1.at(EList)=Ecur[jj][0];
                E_Vertex2.at(EList)=Ecur[jj][1];
                E_Face1.at(EList)=ii;
                
                //mexPrintf("Face: %i, Edge: %i, Add to list %i \n", ii,jj, EList);
                
                EList++;
                if (EList >= (2*FacesN)) {
                // length of Values is greater than previously assigned
                   Er=2;
                }
            }
            
        }
        
    }
    
    // erasing unused space in the vectors
    E_Face1.erase(E_Face1.begin()+EList, E_Face1.end());
    E_Face2.erase(E_Face2.begin()+EList, E_Face2.end());
    E_Vertex1.erase(E_Vertex1.begin()+EList, E_Vertex1.end());
    E_Vertex2.erase(E_Vertex1.begin()+EList, E_Vertex2.end());
    return Er;
}

/*
SILHOUETTE

Finds the edges that are part of the silhouette given a viewing point [x, y, z]

*/
void MeshProject::silhouette_edges(const std::vector<double> & ViewP){
    
    //View position
    ViewPos[0]=ViewP.at(0); 
    ViewPos[1]=ViewP.at(1); 
    ViewPos[2]=ViewP.at(2);
    
    //Projection vector
    double ProjDir[3];
    double F1Dot, F2Dot;
    
    //Faces of edge of interest
    int ff1, ff2;
    //counters
    int ii;
    SList=0;
    
    // empties the vector in case new inputs are used
    Sil.clear(); 
    //assigning space to the vector to hopefully make it run faster
    Sil.reserve(EList/5); 
    
    // for each edge
    for (ii=0; ii<EList; ii++) {
                
        //faces of interest
        ff1=E_Face1.at(ii);
        ff2=E_Face2.at(ii);
        
        //Including boundary points
        // if there is only one face attached to the edge (i.e. a boundary point)
        if ( (ff1==-1) || (ff2==-1) ){
            //mexPrintf("-1 found \n");
            
            // making sure these variables are included
            F1Dot=1;
            F2Dot=-1;
        }
        
        else {
            // First face connected to edge
            ProjDir[0]=FCx[ff1]-ViewPos[0]; ProjDir[1]=FCy[ff1]-ViewPos[1]; ProjDir[2]=FCz[ff1]-ViewPos[2];
            F1Dot=FNx[ff1]*ProjDir[0] + FNy[ff1]*ProjDir[1] + FNz[ff1]*ProjDir[2];
            
            // Second face connected to edge
            ProjDir[0]=FCx[ff2]-ViewPos[0]; ProjDir[1]=FCy[ff2]-ViewPos[1]; ProjDir[2]=FCz[ff2]-ViewPos[2];
            F2Dot=FNx[ff2]*ProjDir[0] + FNy[ff2]*ProjDir[1] + FNz[ff2]*ProjDir[2];
        }
        
        //Opposite signs implies that this edge is a silhouette
        if (F1Dot*F2Dot < 0) {
            
            //assigning new value to v
            Sil.push_back(ii);
            SList++;
            
        }
        
    }
}

/*
SILHOUETTE

Finds the edges that are part of the silhouette given a viewing point [x, y, z]

*/
void MeshProject::silhouette_edges_lnoe(const std::vector<double> & ViewP){
    
    //View position
    ViewPos[0]=ViewP.at(0); 
    ViewPos[1]=ViewP.at(1); 
    ViewPos[2]=ViewP.at(2);
    
    //Projection vector
    double ProjDir[3];
    double ViewTo[3];
    double ProjDirN;
    double F1Dot, F2Dot;
    
    //Faces of edge of interest
    int ff1, ff2;
    //counters
    int ii;
    SList=0;
    
    // empties the vector in case new inputs are used
    Sil.clear(); 
    //assigning space to the vector to hopefully make it run faster
    Sil.reserve(EList/5); 
    
    // for each edge
    for (ii=0; ii<EList; ii++) {
                
        //faces of interest
        ff1=E_Face1.at(ii);
        ff2=E_Face2.at(ii);
        
        // Not including boundary points
        // if there is only one face attached to the edge (i.e. a boundary point)
        if ( (ff1==-1) || (ff2==-1) ){
            //mexPrintf("-1 found \n");
            
            // making sure these variables are not included
            F1Dot=1;
            F2Dot=1;
        }
        
        else {
            // View to point (position between to face centroids)
            ViewTo[0]=(FCx[ff1] + FCx[ff2])/2; 
            ViewTo[1]=(FCy[ff1] + FCy[ff2])/2;
            ViewTo[2]=(FCz[ff1] + FCz[ff2])/2;

            // View direction
            ProjDir[0]=ViewTo[0]-ViewPos[0]; ProjDir[1]=ViewTo[1]-ViewPos[1]; ProjDir[2]=0;
            // Normalising the vector (easier for future analysis)
            ProjDirN=std::sqrt(ProjDir[0]*ProjDir[0] + ProjDir[1]*ProjDir[1] + ProjDir[2]*ProjDir[2]);
            ProjDir[0]=ProjDir[0]/ProjDirN; ProjDir[1]=ProjDir[1]/ProjDirN; ProjDir[2]=ProjDir[2]/ProjDirN;             

            // First face connected to edge
            F1Dot=FNx[ff1]*ProjDir[0] + FNy[ff1]*ProjDir[1] + FNz[ff1]*ProjDir[2];
            
            // Second face connected to edge
            F2Dot=FNx[ff2]*ProjDir[0] + FNy[ff2]*ProjDir[1] + FNz[ff2]*ProjDir[2];
        }
        
        //Opposite signs implies that this edge is a silhouette
       
        if (F1Dot*F2Dot <= 0.005) {
        // Some vertices are missed with cases that are almost parrallel. 
        // Allowing  up to a 0.005 seems to reduce this problem without introducing issues
            
            //assigning new value to v
            Sil.push_back(ii);
            SList++;
            
        }
        
    }
}


//Output the number of edges that make up the silhouette
int MeshProject::edge_silhouette_num(){
    return SList;
}

/*
std::vector<int> MeshProject::OutF1() {
    return E_Face1;
}

std::vector<int> MeshProject::OutF2() {
    return E_Face2;
}

std::vector<int> MeshProject::OutE1() {
    return E_Vertex1;
}

std::vector<int> MeshProject::OutE2() {
    return E_Vertex2;
}
*/


void MeshProject::silhouette_out( std::vector<int> & F1, std::vector<int> & F2, std::vector<int> & V1, std::vector<int> & V2){
    
    int ii;
    
    for (ii=0; ii<SList; ii++) {
        // Each silhouette edge (Sil)
        // The vectors for that edge
        
        F1[ii]=E_Face1.at(Sil.at(ii));
        F2[ii]=E_Face2.at(Sil.at(ii));
        V1[ii]=E_Vertex1.at(Sil.at(ii));
        V2[ii]=E_Vertex2.at(Sil.at(ii));
    }
      
}

