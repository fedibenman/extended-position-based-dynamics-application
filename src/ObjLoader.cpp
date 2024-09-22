
#include "ObjLoader.h"



    void ObjLoader::loadOBJ(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            std::cerr << "Could not open OBJ file: " << filepath << std::endl;
            return;
        }

        std::string line;
    std::cout << line << std::endl ;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string prefix;
            iss >> prefix;

            if (prefix == "v") {
            glm::vec3 vertex;
                iss >> vertex.x >> vertex.y >> vertex.z;
// print the vertex variable

                vertices.push_back(vertex);
            }
            else if (prefix == "vt") {
                glm::vec2 texCoord;
                iss >> texCoord.x >> texCoord.y;
                texCoords.push_back(texCoord);
            }
            // else if (prefix == "vn") {
            //     glm::vec3 normal;
            //     iss >> normal.x >> normal.y >> normal.z;
            //     normals.push_back(normal);
            // }
            else if (prefix == "f") {
                unsigned int vertexIndex[3], texCoordIndex[3], normalIndex[3];
                char slash;
                for (int i = 0; i < 3; i++) {
                    iss >> vertexIndex[i] >> slash >> texCoordIndex[i] >> slash >> normalIndex[i];
                    indices.push_back(vertexIndex[i] - 1);  // OBJ indices start at 1
                }
            }
        }
    }


ObjLoader::ObjLoader(const std::string& filepath) {
    // Load the OBJ and set up buffers
    loadOBJ(filepath);
    setupBuffers();
    
    // Initialize the object's position above the ground
    position = glm::vec3(0.0f, 10.0f, 0.0f); // Start 10 units above the ground
    velocity = glm::vec3(0.0f, 0.0f, 0.0f);  // Initial velocity
}

void ObjLoader::setupBuffers() {
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    // Combine vertex data (positions only) into a single VBO
    std::vector<float> vertexData;
    for (size_t i = 0; i < vertices.size(); ++i) {
    glm::vec3 vertex = vertices[i];
    vertexData.push_back(vertex.x);
    vertexData.push_back(vertex.y);
    vertexData.push_back(vertex.z);

    // Add texture coordinates (uv)
    glm::vec2 texCoord = texCoords[i];  // Make sure indices match tex coords
    vertexData.push_back(texCoord.x);
    vertexData.push_back(texCoord.y);

    }
    // Vertex buffer
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertexData.size() * sizeof(float), vertexData.data(), GL_STATIC_DRAW);

    // Element buffer (for indices)
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

    // Vertex Positions (location = 0), now with correct stride
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // Since you're not using normals or texture coordinates, you don't need to set them up
    // Unbind
    glBindVertexArray(0);
}



void ObjLoader::draw(GLuint shaderProgram) {
    // Create a model matrix for the object's position
    glm::mat4 model = glm::translate(glm::mat4(1.0f), position);

    // Send the model matrix to the shader
    GLuint modelLoc = glGetUniformLocation(shaderProgram, "model");
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

    // Draw the object
    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}

void ObjLoader::updatePhysics(float deltaTime) {
    // Apply gravity to the y-velocity
    velocity.y += gravity * deltaTime;

    // Update the position based on velocity
    position += velocity * deltaTime;

    // Check for collision with the ground
    if (position.y <= groundLevel) {
        position.y = groundLevel;  // Prevent the object from going below the ground
        velocity.y = 0.0f;         // Reset the velocity to 0 when hitting the ground
    }
}



// void ObjLoader::solveDistanceConstraint(ObjLoader& obj2, float restLength, float compliance, float dt) {
//     glm::vec3 deltaPos = position - obj2.position;
//     float currentLength = glm::length(deltaPos);
//     float constraint = currentLength - restLength;

//     // Compute the lagrange multiplier correction
//     float correction = constraint / (invMass + obj2.invMass + compliance / dt / dt);

//     // Apply position corrections
//     glm::vec3 direction = glm::normalize(deltaPos);
//     position -= invMass * correction * direction;
//     obj2.position += obj2.invMass * correction * direction;
// }



    void ObjLoader::deleteBuffers() {
        glDeleteVertexArrays(1, &VAO);
        glDeleteBuffers(1, &VBO);
        glDeleteBuffers(1, &EBO);
    }
