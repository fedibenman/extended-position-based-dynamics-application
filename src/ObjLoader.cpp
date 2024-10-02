
#include "ObjLoader.h"



    void ObjLoader::loadOBJ(const std::string& filepath) {

 
      int faceCounter = 0; 
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
            // vertices.push_back(vertex);
            particles.push_back({ {vertex.x, vertex.y, vertex.z}, {0.0f, 0.0f, 0.0f}, 1.0f });  // Create particles with invMass
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
            unsigned int vertexIndex[4], texCoordIndex[4], normalIndex[4]; // Arrays for indices
            char slash;
            for (int i = 0; i < 4; i++) {
                iss >> vertexIndex[i] >> slash >> texCoordIndex[i] >> slash >> normalIndex[i];
                // Store the vertex index in the indices array, converting to 0-based index
                // indices.push_back(vertexIndex[i] - 1);
            }
                indices.push_back(vertexIndex[0] - 1);
                indices.push_back(vertexIndex[1] - 1);
                indices.push_back(vertexIndex[3] - 1);

    // Triangle 2: v0, v2, v3
                indices.push_back(vertexIndex[1] -1);
                indices.push_back(vertexIndex[3] -1);
                indices.push_back(vertexIndex[2] -1);


            // Create edges for this face (3 unique edges for a triangular face)
            addEdge(vertexIndex[0], vertexIndex[1], faceCounter);
            addEdge(vertexIndex[1], vertexIndex[3], faceCounter);
            addEdge(vertexIndex[0], vertexIndex[3], faceCounter);
            std::array<unsigned int, 3> firstTriangle ={vertexIndex[0] , vertexIndex[1] , vertexIndex[3] } ; 
            addFace(firstTriangle);
            ++faceCounter;

            addEdge(vertexIndex[1], vertexIndex[3], faceCounter);
            addEdge(vertexIndex[1], vertexIndex[2], faceCounter);
            addEdge(vertexIndex[2], vertexIndex[3], faceCounter);
            std::array<unsigned int, 3> secondTriangle ={vertexIndex[1] , vertexIndex[3] , vertexIndex[2] } ;
            addFace(secondTriangle);

            // Check for a possible tetrahedron formation when an edge is shared
            // checkForTetrahedron();
            ++faceCounter;

    std::cout <<  "v1 "<<vertexIndex[0] << "v2 "<<vertexIndex[1] <<"v3 "<<vertexIndex[2] <<"v4 "<<vertexIndex[3] << "nb faces " << faceCounter << std::endl ;

        }
        }
        ObjLoader::findTetrahedrons();

    }


ObjLoader::ObjLoader(const std::string& filepath) {
    // Load the OBJ and set up buffers
    loadOBJ(filepath);
    setupBuffers();
    
    // Initialize the object's position above the ground
    position = glm::vec3(0.0f, 5.0f, 0.0f); // Start 10 units above the ground
}

void ObjLoader::setupBuffers() {
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    // Combine particle position data and texture coordinates into a single VBO
    std::vector<float> particleData;
    for (size_t i = 0; i < particles.size(); ++i) {
        // Access particle position
        glm::vec3 position = particles[i].position; // Assuming particles[i] has a `position` field
        particleData.push_back(position.x);
        particleData.push_back(position.y);
        particleData.push_back(position.z);

        // Add texture coordinates (uv), assuming the texture coordinates match particle indices
        glm::vec2 texCoord = texCoords[i]; // Ensure that the texCoords and particles are aligned
        particleData.push_back(texCoord.x);
        particleData.push_back(texCoord.y);
    }

    // Vertex buffer (particle positions + texture coordinates)
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, particleData.size() * sizeof(float), particleData.data(), GL_STATIC_DRAW);

    // Element buffer (for indices)
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

    // Particle Positions (location = 0), with correct stride
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Texture Coordinates (location = 1), stored after the 3 position floats
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // Unbind VAO
    glBindVertexArray(0);
}

void ObjLoader::updateVBO() {
        std::vector<float> particleData;
    for (size_t i = 0; i < particles.size(); ++i) {
        // Access particle position
        glm::vec3 position = particles[i].position; // Assuming particles[i] has a `position` field
        particleData.push_back(position.x);
        particleData.push_back(position.y);
        particleData.push_back(position.z);

        // Add texture coordinates (uv), assuming the texture coordinates match particle indices
        glm::vec2 texCoord = texCoords[i]; // Ensure that the texCoords and particles are aligned
        particleData.push_back(texCoord.x);
        particleData.push_back(texCoord.y);

        
    }
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, particleData.size() * sizeof(float), particleData.data(), GL_STATIC_DRAW);
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
    glm::vec3 gravity = glm::vec3(0.0f, -9.81f, 0.0f);  // Gravity force pulling particles downwards

for (const Edge& edge : edges) {
    glm::vec3& pos1 = particles[edge.v1].position;
    glm::vec3& pos2 = particles[edge.v2].position;

    float currentLength = glm::distance(pos1, pos2);
    float stretch = currentLength - edge.restLen;

    std::cout << "First particle pos: (" << pos1.x << ", " << pos1.y << ", " << pos1.z << ")" << std::endl;
    std::cout << "Second particle pos: (" << pos2.x << ", " << pos2.y << ", " << pos2.z << ")" << std::endl;
    std::cout << "Current Length: " << currentLength << ", Rest Length: " << edge.restLen << std::endl;
    std::cout << "Stretch: " << stretch << std::endl;

    if (fabs(stretch) < 1e-4f) {
        std::cout << "Stretch is too small, skipping correction." << std::endl;
        continue;
    }

    glm::vec3 correctionDir = pos2 - pos1;
    if (glm::length(correctionDir) < 1e-4f) {
        std::cout << "Correction direction is too small, skipping correction." << std::endl;
        continue;
    }

    correctionDir = glm::normalize(correctionDir);
    glm::vec3 correction = (stretch * correctionDir) / 2.0f;

    // Check for NaN
    if (std::isnan(correction.x) || std::isnan(correction.y) || std::isnan(correction.z) ||
        std::isnan(pos1.x) || std::isnan(pos1.y) || std::isnan(pos1.z) ||
        std::isnan(pos2.x) || std::isnan(pos2.y) || std::isnan(pos2.z)) {
        std::cerr << "NaN detected in correction or particle positions!\n";
        continue;  // Skip this iteration
    }

    // Cap the correction
    const float MAX_CORRECTION = 0.1f;  // Example maximum correction
    correction = glm::clamp(correction, -MAX_CORRECTION, MAX_CORRECTION);

    // Apply the correction to particle positions
    particles[edge.v1].position += correction;
    particles[edge.v2].position -= correction;

    // Print updated positions
    std::cout << "Updated pos1: (" << particles[edge.v1].position.x << ", " 
              << particles[edge.v1].position.y << ", " 
              << particles[edge.v1].position.z << ")" << std::endl;
    std::cout << "Updated pos2: (" << particles[edge.v2].position.x << ", " 
              << particles[edge.v2].position.y << ", " 
              << particles[edge.v2].position.z << ")" << std::endl;
}

//     // Step 2: Pre-solve Volume Constraints (for each tetrahedron)
// for (tetrahedron& tetra : tetrahedrons) {
//     glm::vec3& pos1 = particles[tetra.indices[0]].position;
//     glm::vec3& pos2 = particles[tetra.indices[1]].position;
//     glm::vec3& pos3 = particles[tetra.indices[2]].position;
//     glm::vec3& pos4 = particles[tetra.indices[3]].position;

//     float currentVolume = glm::dot(pos1 - pos4, glm::cross(pos2 - pos4, pos3 - pos4)) / 6.0f;

//     // Check for degenerate tetrahedron
//     if (currentVolume != 0.0f) {  // Avoid correcting tetrahedrons with zero volume
//         float volumeError = currentVolume - tetra.restVolume;
//         glm::vec3 centroid = (pos1 + pos2 + pos3 + pos4) / 4.0f;

//         glm::vec3 correctionFactor = (volumeError / 4.0f) * glm::normalize(centroid);

//         particles[tetra.indices[0]].position += correctionFactor;
//         particles[tetra.indices[1]].position += correctionFactor;
//         particles[tetra.indices[2]].position += correctionFactor;
//         particles[tetra.indices[3]].position += correctionFactor;
//     }
// }

for (Particle& particle : particles) {
        // Apply gravity to the particle's velocity
        particle.velocity += 2.0f * gravity * deltaTime;

        // Update particle's position based on velocity
        particle.position += particle.velocity * deltaTime;

        // Check for collision with the ground (y = groundLevel)
        if (particle.position.y < this->groundLevel) {
            // Handle ground collision: Reset position to ground level
            particle.position.y = this->groundLevel;

            // Reverse or dampen the velocity to simulate a bounce or friction
            // You can tweak the restitution (bounce factor) to adjust the bounce behavior
            float restitution = 0.5f;  // Adjust this value (between 0 and 1) for bounce effect
            particle.velocity.y *= -restitution;  // Reverse y-velocity with damping
        }
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


float ObjLoader::calculateTetrahedronVolume(const std::array<unsigned int, 4>& tetraVertices) {
    glm::vec3 v0 = particles[tetraVertices[0]].position;
    glm::vec3 v1 = particles[tetraVertices[1]].position;
    glm::vec3 v2 = particles[tetraVertices[2]].position;
    glm::vec3 v3 = particles[tetraVertices[3]].position;

    // Compute the vectors for three edges
    glm::vec3 vec1 = v1 - v0;
    glm::vec3 vec2 = v2 - v0;
    glm::vec3 vec3 = v3 - v0;

    // Calculate the volume using the scalar triple product
    float volume = glm::dot(vec1, glm::cross(vec2, vec3)) / 6.0f;

    return std::abs(volume);  // Return the absolute value of the volume
}





void ObjLoader::addEdge( unsigned int v1, unsigned int v2, int face) {
    glm::vec3 vertex1 = particles[v1].position;
    glm::vec3 vertex2 = particles[v2].position;

    float restLen = glm::distance(vertex1, vertex2);
    Edge edge(v1, v2, restLen);
    auto it = edges.find(edge);
    
    if (it == edges.end()) {
        // Edge doesn't exist yet, insert it with the current face
        edges.insert(edge);
    } 
}


void ObjLoader::findTetrahedrons() {
    for (int i = 0; i < faces.size(); ++i) {
        const Face& face = faces[i];

        // Check if the face has at least two adjacent faces
        if (face.adjacentFaces.size() >= 2) {
            for (int j = 0; j < face.adjacentFaces.size(); ++j) {
                for (int k = j + 1; k < face.adjacentFaces.size(); ++k) {
                    const Face& adjFace1 = faces[face.adjacentFaces[j]];
                    const Face& adjFace2 = faces[face.adjacentFaces[k]];

                    // Check if the two adjacent faces are also adjacent to each other
                    if (adjFace1.isAdjacentTo(adjFace2)) {
                        // We've found a set of faces that form a tetrahedron
                        std::set<unsigned int> uniqueVertices;
                        uniqueVertices.insert(face.vertices.begin(), face.vertices.end());
                        uniqueVertices.insert(adjFace1.vertices.begin(), adjFace1.vertices.end());
                        uniqueVertices.insert(adjFace2.vertices.begin(), adjFace2.vertices.end());

                        // If there are exactly 4 unique vertices, we have a tetrahedron
                        if (uniqueVertices.size() == 4) {
                            // Convert the set to an array for easy handling
                            std::array<unsigned int, 4> tetraVertices;
                            std::copy(uniqueVertices.begin(), uniqueVertices.end(), tetraVertices.begin());

                            // Calculate the rest volume of the tetrahedron
                            float restVolume = calculateTetrahedronVolume(tetraVertices);

                            // Create and store the tetrahedron
                            tetrahedrons.push_back(tetrahedron(tetraVertices, restVolume));

                            // Output for debugging
                            std::cout << "Tetrahedron found with vertices: ";
                            for (unsigned int v : tetraVertices) {
                                std::cout << v << " ";
                            }
                            std::cout << " and volume: " << restVolume << std::endl;
                        }
                    }
                }
            }
        }
    }
}









    // void ObjLoader::checkForTetrahedron() {
    // std::map<std::pair <int,int>, std::vector<Edge>> sharedFaceEdges; // Map face pairs to shared edges

    // for (const Edge& edge : edges) {
    //     // if (edge.isShared()) {
    //     //             for (int face : edge.faces) {
    //     // std::cout << face << " ";
    //     //     }
    //     //         std::cout << "," << edge.v1 << "," << edge.v2 << std::endl;
    //         // Loop through all pairs of faces sharing this edge
    //         for (size_t i = 0; i < edge.faces.size(); ++i) {
    //             for (size_t j = i + 1; j < edge.faces.size(); ++j) {
    //                 int face1 = edge.faces[i];
    //                 int face2 = edge.faces[j];

    //                 // Store the edge under the face pair (face1, face2)
    //                 sharedFaceEdges[{face1, face2}].push_back(edge);
    //             }
    //         }
    //     }
    // }

    // // Now find valid tetrahedrons from the shared face edges
    // for (const auto& entry : sharedFaceEdges) {
    //     const std::vector<Edge>& edges = entry.second;

    //     // Check if there are 6 shared edges (which should form a tetrahedron)
    //     if (edges.size() == 6) {
    //         // Collect unique vertices from these 6 edges
    //         std::set<unsigned int> uniqueVertices;
    //         for (const Edge& edge : edges) {
    //             uniqueVertices.insert(edge.v1);
    //             uniqueVertices.insert(edge.v2);
    //         }

    //         // A tetrahedron should have exactly 4 unique vertices
    //         if (uniqueVertices.size() == 4) {
    //             // Convert the set of unique vertices to an array
    //             unsigned int tetraVertices[4];
    //             int i = 0;
    //             for (unsigned int v : uniqueVertices) {
    //                 tetraVertices[i++] = v;
    //             }

    //             // Store the tetrahedron and create its edges
    //             tetrahedrons.push_back({{tetraVertices[0], tetraVertices[1], tetraVertices[2], tetraVertices[3]} , calculateVolume(tetraVertices)});
    //         }
    //     }
    // }




void ObjLoader::addFace(const std::array<unsigned int, 3>& vertices) {
    glm::vec3 vertex1 = this->particles[vertices[0]].position;
    glm::vec3 vertex2 = this->particles[vertices[1]].position;
    glm::vec3 vertex3 = this->particles[vertices[2]].position;

    float restLen1 = glm::distance(vertex1, vertex2);
    float restLen2 = glm::distance(vertex2, vertex3);
    float restLen3 = glm::distance(vertex3, vertex1);

    // Create the edges with their respective rest lengths
    std::vector<Edge> faceEdges = {
        Edge(vertices[0], vertices[1], restLen1),
        Edge(vertices[1], vertices[2], restLen2),
        Edge(vertices[2], vertices[0], restLen3)
    };
    Face newFace(vertices, faceEdges);

    // Check for adjacency with previously added faces
    for (int i = 0; i < faces.size(); ++i) {
        if (newFace.isAdjacentTo(faces[i])) {
            newFace.addAdjacentFace(i);
            faces[i].addAdjacentFace(faces.size());  // Mutual adjacency
        }
    }

    // Add the new face to the list of faces
    faces.push_back(newFace);
}
