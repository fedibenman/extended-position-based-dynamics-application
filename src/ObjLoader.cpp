
#include "ObjLoader.h"



void ObjLoader::loadOBJ(const std::string& baseFilepath) {
    // File paths
    std::string nodeFile = baseFilepath + ".node";
    std::string eleFile = baseFilepath + ".ele";
    std::string edgeFile = baseFilepath + ".edge";  // New edge file
    std::string faceFile = baseFilepath + ".face";  // Face file for surface mesh
    
    std::ifstream nodeStream(nodeFile);
    std::ifstream eleStream(eleFile);
    std::ifstream edgeStream(edgeFile);
    std::ifstream faceStream(faceFile);
    
    if (!nodeStream.is_open()) {
        std::cerr << "Could not open TetGen node file: " << nodeFile << std::endl;
        return;
    }
    if (!eleStream.is_open()) {
        std::cerr << "Could not open TetGen element file: " << eleFile << std::endl;
        return;
    }
    if (!edgeStream.is_open()) {
        std::cerr << "Could not open TetGen edge file: " << edgeFile << std::endl;
        return;
    }
    if (!faceStream.is_open()) {
        std::cerr << "Could not open TetGen face file: " << faceFile << std::endl;
        return;
    }

    // --- Read the .node file (particles: vertex positions) ---
    std::string line;
    int numVertices, dimensions, numAttributes, numBoundaryMarkers;
    
    std::getline(nodeStream, line);  // First line: number of vertices and file info
    std::istringstream nodeHeader(line);
    nodeHeader >> numVertices >> dimensions >> numAttributes >> numBoundaryMarkers;

    for (int i = 0; i < numVertices; ++i) {
        unsigned int index;
        glm::vec3 vertex;
        std::getline(nodeStream, line);
        std::istringstream vertexStream(line);
        vertexStream >> index >> vertex.x >> vertex.y >> vertex.z;

        // Create particles with position, velocity, and inverse mass
        particles.push_back({ vertex, glm::vec3(0.0f), 1.0f});
        
        // Generate spherical texture coordinates for the vertex
        glm::vec2 texCoord = generateSphericalTexCoord(vertex);
        texCoords.push_back(texCoord);
        
        std::cout << "Loaded Particle " << index << ": Position (" << vertex.x << ", " << vertex.y << ", " << vertex.z << "), TexCoord (" << texCoord.x << ", " << texCoord.y << ")" << std::endl;
    }

    // --- Read the .ele file (tetrahedrons) ---
    int numTetrahedrons, pointsPerTet, numAttributesPerTet;
    
    std::getline(eleStream, line);  // First line: number of tetrahedrons and file info
    std::istringstream eleHeader(line);
    eleHeader >> numTetrahedrons >> pointsPerTet >> numAttributesPerTet;

    if (pointsPerTet != 4) {
        std::cerr << "Expected 4 vertices per tetrahedron, found: " << pointsPerTet << std::endl;
        return;
    }

    for (int i = 0; i < numTetrahedrons; ++i) {
        unsigned int index, v1, v2, v3, v4;
        std::getline(eleStream, line);
        std::istringstream tetStream(line);
        tetStream >> index >> v1 >> v2 >> v3 >> v4;

        // Add tetrahedron
        std::array<unsigned int, 4> vertices = {v1, v2, v3, v4};  // Adjust to 0-based indexing
        float volume = calculateTetrahedronVolume(particles[v1].position, particles[v2].position, particles[v3].position, particles[v4].position);
        tetrahedrons.push_back(tetrahedron(vertices, volume));

        std::cout << "Loaded Tetrahedron " << index << ": vertices (" << v1 << ", " << v2 << ", " << v3 << ", " << v4 << "), volume: " << volume << std::endl;
    }

    // --- Read the .edge file (edges) ---
    int numEdges, boundaryMarkerFlag;
    std::getline(edgeStream, line);  // First line: number of edges and file info
    std::istringstream edgeHeader(line);
    edgeHeader >> numEdges >> boundaryMarkerFlag;

    for (int i = 0; i < numEdges; ++i) {
        unsigned int index, v1, v2;
        std::getline(edgeStream, line);
        std::istringstream edgeStream(line);
        edgeStream >> index >> v1 >> v2;

        // Add edge
        float restLength = glm::distance(particles[v1].position, particles[v2].position);  // Compute edge length
        edges.insert(Edge(v1 , v2 , restLength));

        std::cout << "Loaded Edge " << index << ": vertices (" << v1 << ", " << v2 << "), rest length: " << restLength << std::endl;
    }

    // --- Read the .face file (surface mesh) ---
    int numFaces;
    std::getline(faceStream, line);  // First line: number of faces and file info
    std::istringstream faceHeader(line);
    faceHeader >> numFaces;

    for (int i = 0; i < numFaces; ++i) {
        unsigned int index, v1, v2, v3;
        std::getline(faceStream, line);
        std::istringstream faceStream(line);
        faceStream >> index >> v1 >> v2 >> v3;

        indices.push_back(v1);
        indices.push_back(v2);
        indices.push_back(v3);

        std::cout << "Loaded Face " << index << ": vertices (" << v1 << ", " << v2 << ", " << v3 << ")" << std::endl;
    }

    std::cout << "Mesh loading complete." << std::endl;
}


// Function to compute the volume of a tetrahedron given its 4 vertices
float ObjLoader::calculateTetrahedronVolume(const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& v3, const glm::vec3& v4) {
    glm::vec3 vec1 = v1 - v4;
    glm::vec3 vec2 = v2 - v4;
    glm::vec3 vec3 = v3 - v4;
    return std::abs(glm::dot(vec1, glm::cross(vec2, vec3)) / 6.0f);

}


glm::vec2 ObjLoader::generateSphericalTexCoord(const glm::vec3& vertex) {
    // Normalize the vertex position
    glm::vec3 normalizedPos = glm::normalize(vertex);

    // Calculate spherical coordinates (longitude and latitude)
    float u = 0.5f + (atan2(normalizedPos.z, normalizedPos.x) / (2.0f * M_PI));  // Longitude
    float v = 0.5f - (asin(normalizedPos.y) / M_PI);  // Latitude

    return glm::vec2(u, v);  // Return as texture coordinates
}


ObjLoader::ObjLoader(const std::string& baseFilepath) {
    // Load the OBJ and set up buffers
    loadOBJ(baseFilepath);
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
    glm::vec3 gravity = glm::vec3(0.0f, -9.81f, 0.0f);  
for (const Edge& edge : edges) {

    glm::vec3& pos1 = particles[edge.v1].position;
    glm::vec3& pos2 = particles[edge.v2].position;
 
    float currentLength = glm::distance(pos1, pos2);
    float stretch = currentLength - edge.restLen;

    if (fabs(stretch) < 1e-4f) {
        std::cout << "Stretch is too small, skipping correction." << std::endl;
        continue;
    }

    glm::vec3 correctionDir = pos2 - pos1;
    // if (glm::length(correctionDir) < 1e-4f) {
    //     std::cout << "Correction direction is too small, skipping correction." << std::endl;
    //     continue;
    // }

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


}




for (tetrahedron& tetra : tetrahedrons) {
        glm::vec3& pos1 = particles[tetra.indices[0]].position;
        glm::vec3& pos2 = particles[tetra.indices[1]].position; 
        glm::vec3& pos3 = particles[tetra.indices[2]].position;
        glm::vec3& pos4 = particles[tetra.indices[3]].position;

        // Calculate current volume using the determinant formula
        float currentVolume = std::abs(glm::dot(pos1 - pos4, glm::cross(pos2 - pos4, pos3 - pos4))) / 6.0f;

        // Calculate the volume error relative to the rest volume
        float volumeError = (currentVolume - tetra.restVolume);

        // Compute gradients of the volume constraint with respect to each position
        glm::vec3 grad_p1 = glm::cross(pos4 - pos2, pos3 - pos2);
        glm::vec3 grad_p2 = glm::cross(pos3 - pos1, pos4 - pos1);
        glm::vec3 grad_p3 = glm::cross(pos4 - pos1, pos2 - pos1);
        glm::vec3 grad_p4 = glm::cross(pos2 - pos1, pos3 - pos1);

        // Compute individual particle corrections
        float stiffness = 0.5e14f;  // Adjust this value for desired stiffness
        float minDistance = 0.01f;  // Minimum distance between particles

        glm::vec3 correction_p1 = stiffness * grad_p1 / (glm::abs(glm::dot(grad_p1, grad_p1)) + minDistance);
        glm::vec3 correction_p2 = stiffness * grad_p2 / (glm::abs(glm::dot(grad_p2, grad_p2)) + minDistance);
        glm::vec3 correction_p3 = stiffness * grad_p3 / (glm::abs(glm::dot(grad_p3, grad_p3)) + minDistance);
        glm::vec3 correction_p4 = stiffness * grad_p4 / (glm::abs(glm::dot(grad_p4, grad_p4)) + minDistance);

        // Apply corrections with spatial relationship consideration
        glm::vec3 avgCorrection = (correction_p1 + correction_p2 + correction_p3 + correction_p4) / 4.0f;
        
        // Normalize the average correction vector
        glm::vec3 normalizedAvgCorrection = glm::normalize(avgCorrection);

        // Calculate final corrections considering spatial relationships
        float spatialFactor = 0.001f;  // Adjust based on tetrahedron orientation
        glm::vec3 finalCorrection_p1 = normalizedAvgCorrection * glm::length(pos1 - pos4) * spatialFactor;
        glm::vec3 finalCorrection_p2 = normalizedAvgCorrection * glm::length(pos2 - pos4) * spatialFactor;
        glm::vec3 finalCorrection_p3 = normalizedAvgCorrection * glm::length(pos3 - pos4) * spatialFactor;
        glm::vec3 finalCorrection_p4 = normalizedAvgCorrection * glm::length(pos4 - pos1) * spatialFactor;

        // Clamp corrections to prevent excessive movement
        const float MAX_CORRECTION = 1.0f;
        finalCorrection_p1 = glm::clamp(finalCorrection_p1, -MAX_CORRECTION, MAX_CORRECTION);
        finalCorrection_p2 = glm::clamp(finalCorrection_p2, -MAX_CORRECTION, MAX_CORRECTION);
        finalCorrection_p3 = glm::clamp(finalCorrection_p3, -MAX_CORRECTION, MAX_CORRECTION);
        finalCorrection_p4 = glm::clamp(finalCorrection_p4, -MAX_CORRECTION, MAX_CORRECTION);

        // Apply corrections to particle positions
        particles[tetra.indices[0]].position += finalCorrection_p1;
        particles[tetra.indices[1]].position += finalCorrection_p2;
        particles[tetra.indices[2]].position -= finalCorrection_p3;
        particles[tetra.indices[3]].position -= finalCorrection_p4;

        std::cout << "After update" << std::endl;
        std::cout << pos1.x << ", " << pos1.y << ", " << pos1.z << std::endl;
        std::cout << pos2.x << ", " << pos2.y << ", " << pos2.z << std::endl;
        std::cout << pos3.x << ", " << pos3.y << ", " << pos3.z << std::endl;
        std::cout << pos4.x << ", " << pos4.y << ", " << pos4.z << std::endl;
}


for (Particle& particle : particles) {  
    // Apply gravity to the particle's velocity
    particle.velocity +=  gravity * deltaTime;

    // Update particle's position based on velocity
    particle.position += particle.velocity * deltaTime;

    // Check for collision with the ground (y = groundLevel)
    if (particle.position.y < this->groundLevel) {
        // Handle ground collision: Reset position to ground level
        particle.position.y = this->groundLevel;

        float restitution = 0.5f;  // Adjust this value (between 0 and 1) for bounce effect

        // If the vertical velocity is below a small threshold, stop the particle
        if (fabs(particle.velocity.y) < 0.1f) {
            particle.velocity.y = 0.0f;  // Stop bouncing in the y-axis
        } else {
            particle.velocity.y *= -restitution;  // Reverse y-velocity with damping
        }

        // Apply friction to the horizontal movement (x, z)
        // float friction = 0.95f;  // Friction factor to slow down horizontal movement
        // particle.velocity.x *= friction;
        // particle.velocity.z *= friction;

        // If horizontal velocity is very small, stop the particle
        if (glm::length(glm::vec2(particle.velocity.x, particle.velocity.z)) < 0.5f) {
            particle.velocity.x = 0.0f;
            particle.velocity.z = 0.0f;
        }
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





