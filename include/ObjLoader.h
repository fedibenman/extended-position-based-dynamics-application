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

struct Particle{
float  position[3] ;
float velocity[3] ; 
float invMass ; 
}; 

struct Edges{
    unsigned int indices[2] ;
    float restLen ; 
} ; 

struct tetrahedron{
    unsigned int indices[4] ; 
    float restVolume ; 
} ;  
class ObjLoader {
 public:
    glm::vec3 position;
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec2> texCoords;
    std::vector<unsigned int> indices;
    std::vector<Particle> particles;
    std::vector<Edges> edges;
    std::vector<tetrahedron> tetrahedrons;
    float gravity = -9.81f; // Gravity constant
    float groundLevel = 0.0f; // Height of the flat surface (y = 0)


    // static const float gravity = -9.81f; // Gravity constant

    GLuint VBO, VAO, EBO;
ObjLoader(const std::string& filepath) ;  
void loadOBJ(const std::string& filepath)  ; 

void setupBuffers() ; 

void updatePhysics(float deltaTime);

 void draw(GLuint shaderProgram)  ;
 
 void deleteBuffers() ;

 float calculateRestLength(unsigned int index1, unsigned int index2) ;
 float calculateVolume(unsigned int indices[4]) ;
 void  createEdges(unsigned int vertexIndex[4])  ;
} ;

#endif