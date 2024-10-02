#ifndef OBJ_CLASS_H
#define OBJ_CLASS_H
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream> // Use GLM for vector math
#include<glad/glad.h>
#include <glm/glm.hpp> 
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <set>
#include <map>
#include <array>
struct Particle{
glm::vec3  position ;
glm::vec3 velocity ; 
float invMass ; 
}; 

struct Edge {
    unsigned int v1, v2;
    float restLen; 

    Edge(unsigned int v1, unsigned int v2  ,float restLen) : v1(v1), v2(v2) ,restLen(restLen) {
        if (v1 > v2) std::swap(v1, v2);  
    }

    bool operator==(const Edge& other) const {
        return v1 == other.v1 && v2 == other.v2;
    }

    bool operator<(const Edge& other) const {
        return std::tie(v1, v2) < std::tie(other.v1, other.v2);
    }
};
struct Face {
    std::array<unsigned int, 3> vertices;  // Vertices defining the triangular face
    std::vector<Edge> edges;               // Edges of the face
    std::vector<int> adjacentFaces;        // Indices of adjacent faces

    Face(const std::array<unsigned int, 3>& verts, const std::vector<Edge>& edgeList)
        : vertices(verts), edges(edgeList) {}

    // Add an adjacent face
    void addAdjacentFace(int faceIndex) {
        adjacentFaces.push_back(faceIndex);
    }

    // Check if this face shares an edge with another face
    bool isAdjacentTo(const Face& otherFace) const {
        for (const Edge& edge : edges) {
            for (const Edge& otherEdge : otherFace.edges) {
                if (edge == otherEdge) {
                    return true;  // They share an edge
                }
            }
        }
        return false;
    }
};
struct tetrahedron{
    unsigned int indices[4] ; 
    float restVolume ; 

        tetrahedron(const std::array<unsigned int, 4>& vertices, float volume) {
        std::copy(vertices.begin(), vertices.end(), indices);
        restVolume = volume;
    }
} ;  
class ObjLoader {
 public:
    glm::vec3 position;
    // std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec2> texCoords;
    std::vector<unsigned int> indices;
    std::vector<Particle> particles;
    std::set<Edge> edges;
    std::vector<tetrahedron> tetrahedrons;
    std::vector<Face> faces;
    float gravity = -9.81f; // Gravity constant
    float groundLevel = 0.0f; // Height of the flat surface (y = 0)


    // static const float gravity = -9.81f; // Gravity constant

    GLuint VBO, VAO, EBO;
ObjLoader(const std::string& filepath) ;  
void loadOBJ(const std::string& filepath)  ; 

void setupBuffers() ; 

void updatePhysics(float deltaTime) ;

 void draw(GLuint shaderProgram)  ;
 
 void deleteBuffers() ;

 float calculateTetrahedronVolume(const std::array<unsigned int, 4>& tetraVertices) ;
 void  createEdges(unsigned int vertexIndex[4])  ;
void addEdge( unsigned int v1, unsigned int v2, int face) ;
void findTetrahedrons() ;
void addFace(const std::array<unsigned int, 3>& vertices) ; 
void updateVBO() ;
} ;

#endif