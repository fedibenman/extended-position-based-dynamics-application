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
float  postion[3] ;
float velocity[3] ; 
float invMass ; 
}; 

struct Edges{
    int indices[2] ;
    float restLen ; 
} ; 

struct tetrahedron{
    int indices[4] ; 
    float restVolume ; 
} ;  
class ObjLoader {
 public:
    glm::vec3 position;
    glm::vec3 velocity;
    glm::vec3 force;
    float mass;
    float invMass;
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec2> texCoords;
    std::vector<unsigned int> indices;
    float stiffness;  
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
void solveDistanceConstraint( ObjLoader& obj2, float restLength, float compliance, float dt) ; 
} ;

#endif