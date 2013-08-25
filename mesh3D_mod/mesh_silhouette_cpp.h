#ifndef MESHP_H
#define MESHP_H 

#include <cmath>
#include <vector>
#include <iostream>

using std::vector;

class MeshProject {

    private:
        //inputs
        std::vector<double> Fa, Fb, Fc, Vx, Vy, Vz;
        int FacesN, VerticesN;
        
        //face normals
        double * FNx, * FNy, * FNz;
        
        //face centres
        double * FCx, * FCy, * FCz;
        
        //Edges (related vertices and faces)
        int EList; //Number of vertices found.
        vector<int> E_Vertex1, E_Vertex2, E_Face1, E_Face2;
        
        //Silhouette
        double ViewPos [3];
        //Silhouette vertex number, 
        //(Note that using vector class handles dynamic memory management and is safer if vector methods are used)
        vector<int> Sil;
        
        //Number of silhouette edges
        int SList; 
        

    public:
        //contructor
        MeshProject(std::vector<double> & Fa1, std::vector<double> & Fb1, std::vector<double> & Fc1, std::vector<double> & Vx1, std::vector<double> & Vy1, std::vector<double> & Vz1);
        //deconstructor
        ~MeshProject();
        
        int update_vertices( std::vector<double> & Vx1, std::vector<double> & Vy1, std::vector<double> & Vz1);

        //Find face normals
        int normals_face();
        
        //Find unique edges and related faces
        int unique_edges();
        
        //Find silhouette
        void silhouette_edges(const std::vector<double> & i);
        
        //Find silhouette (no end points and lodox projection style)
        void silhouette_edges_lnoe(const std::vector<double> & i);
        
        int edge_silhouette_num();
        
        int another_test();
        
        double test1 ();
        
        //Output face normals
        
        void silhouette_out( std::vector<int> & F1, std::vector<int> & F2, std::vector<int> & V1, std::vector<int> & V2);
        /*
        std::vector<int> OutF1();
        std::vector<int> OutF2();
        std::vector<int> OutE1();
        std::vector<int> OutE2();
        */
        //void silhouette_out( double * a, double * b, double * c, double * d);
        

};

#endif
