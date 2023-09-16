#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>


float cameraSensitivity = 0.1f;
double lastX = 400.0f;
double lastY = 300.0f;
float cameraSpeed = 0.05f;
bool isMousePressed = false;
double lastMouseX = 0.0; // Initial X position of the mouse (center of the window)
double lastMouseY = 0.0; // Initial Y position of the mouse (center of the window)
bool mouseButtonPressed = false;   // Define sphere vertices
const int sectors = 6;
const int stacks = sectors / 2;
const int numberObjects = 25;


const char* vertexShaderSource = R"(
	
	#version 450 core
        
	layout (location = 0) in vec3 aPos;
	layout (location = 1) in vec4 aColor;
        
	uniform mat4 model;
    uniform mat4 view;
    uniform mat4 projection;

	out vec4 FragColor;

	void main() {
 	
		gl_Position = projection * view * model * vec4(aPos, 1.0);
		FragColor = aColor;
    
    }
    
)";

const char* fragmentShaderSource = R"(

	#version 450 core

	in vec4 FragColor;
    
    out vec4 FragColorOutput;

    void main() {
    
        FragColorOutput = vec4(FragColor);

	}
    
)";


std::vector<float> order(float n1, float n2, float n3) {

	std::vector<float> tempVector;

	tempVector.push_back(n1);
	tempVector.push_back(n2);
	tempVector.push_back(n3);

	for (int i = 0; i < tempVector.size() - 1; i++) {

		for (int j = 0; j < tempVector.size() - i - 1; j++) {
					
			if (tempVector[j] > tempVector[j + 1]) {

				float temp = tempVector[j];
				tempVector[j] = tempVector[j + 1];
				tempVector[j + 1] = temp;

			}
				
		}

	}

	return tempVector;

}

class Model {

	public:

		std::string filename;	
		std::vector<float> vertices;
    	std::vector<unsigned int> indices;
		std::vector<float> color;
		float mass = 10000000000;
		float e = 1.0;
		float minX = 0.0;
		float maxX = 0.0;
		float minY = 0.0;
		float maxY = 0.0;
		float minZ = 0.0;
		float maxZ = 0.0;
		GLuint vertexVAO, vertexVBO, vertexEBO, colorVAO, colorVBO, colorEBO;


		void initModel() {

			// Load the Leviathan.obj model using Assimp
			Assimp::Importer importer;
    		const aiScene* scene = importer.ReadFile(filename, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_GenNormals);

    		if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
   
   				std::cerr << "Assimp error: " << importer.GetErrorString() << std::endl;
   
   			}	

			aiMesh* mesh = scene->mMeshes[0];


			for (unsigned int i = 0; i < mesh->mNumVertices; ++i) {
   
   				vertices.push_back(mesh->mVertices[i].x);
        		vertices.push_back(mesh->mVertices[i].y);
        		vertices.push_back(mesh->mVertices[i].z);

				if (mesh->mVertices[i].x > maxX) {

					maxX = mesh->mVertices[i].x;

				}
				 
				if (mesh->mVertices[i].x < minX) {

					minX = mesh->mVertices[i].x;

				}

				if (mesh->mVertices[i].y > maxY) {

					maxY = mesh->mVertices[i].y;

				}
				 
				if (mesh->mVertices[i].y < minY) {

					minY = mesh->mVertices[i].y;

				}

				if (mesh->mVertices[i].z > maxZ) {

					maxZ = mesh->mVertices[i].z;

				}
				 
				if (mesh->mVertices[i].z < minZ) {

					minZ = mesh->mVertices[i].z;

				}

				color.push_back(1.0);
				color.push_back(1.0); 
				color.push_back(1.0);
				color.push_back(1.0);
   
   			}

    		for (unsigned int i = 0; i < mesh->mNumFaces; ++i) {
   
   				aiFace face = mesh->mFaces[i];
   	
   				for (unsigned int j = 0; j < face.mNumIndices; ++j) {
   
   					indices.push_back(face.mIndices[j]);
   
   				}
   
   			}	

    		glGenVertexArrays(1, &vertexVAO);
    		glBindVertexArray(vertexVAO);
			glGenBuffers(1, &vertexVBO);
			glGenBuffers(1, &vertexEBO);
    		glBindBuffer(GL_ARRAY_BUFFER, vertexVBO);
    		glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

	 		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vertexEBO);
    		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

    		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    		glEnableVertexAttribArray(0);


			glGenBuffers(1, &colorVBO);
    		glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
    		glBufferData(GL_ARRAY_BUFFER, color.size() * sizeof(float), color.data(), GL_STATIC_DRAW);

       		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), (GLvoid*)0);
     	   	glEnableVertexAttribArray(1);

			glBindVertexArray(0);

		}


		void displayModel() {

			glBindVertexArray(vertexVAO);
        	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);	
			glBindVertexArray(0);

		}


		void cleanUpModel() {

			glDeleteVertexArrays(1, &vertexVAO);
    		glDeleteBuffers(1, &vertexVBO);
    		glDeleteBuffers(1, &vertexEBO);

		}

};


class Axis {

	public: 

		float mass = 1000000000;

		GLfloat verticesAxis[18] = {

			  0.0,  0.0,  0.0, 
			 10.0,  0.0,  0.0,

			  0.0,  0.0,  0.0,
			  0.0, 10.0,  0.0,

			  0.0,  0.0,  0.0,
			  0.0,  0.0, 10.0

		};

		GLfloat colorAxis[24] {

			1.0f, 0.0f, 0.0f, 1.0f,
			1.0f, 0.0f, 0.0f, 1.0f,

			0.0f, 1.0f, 0.0f, 1.0f,
			0.0f, 1.0f, 0.0f, 1.0f,

			0.0f, 0.0f, 1.0f, 1.0f,
			0.0f, 0.0f, 1.0f, 1.0f

		};

		GLuint VAO, vertexVBO, colorVBO;


		void initAxis() {
	
			glGenVertexArrays(1, &VAO);
    	    glBindVertexArray(VAO);

        	glGenBuffers(1, &vertexVBO);
        	glBindBuffer(GL_ARRAY_BUFFER, vertexVBO);
        	glBufferData(GL_ARRAY_BUFFER, sizeof(verticesAxis), verticesAxis, GL_STATIC_DRAW);
        	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
	        glEnableVertexAttribArray(0);

        	glGenBuffers(1, &colorVBO);
        	glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
        	glBufferData(GL_ARRAY_BUFFER, sizeof(colorAxis), colorAxis, GL_STATIC_DRAW);
       		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), (GLvoid*)0);
     	   	glEnableVertexAttribArray(1);

        	glBindVertexArray(0);
		
		}


		void displayAxis() {

			glBindVertexArray(VAO);
			glDrawArrays(GL_LINES, 0, 2);
			glDrawArrays(GL_LINES, 2, 2);
			glDrawArrays(GL_LINES, 4, 2);
			glBindVertexArray(0);

		}


		void cleanUpAxis() {

			glDeleteVertexArrays(1, &VAO);

		}

};


class ModelBox {

	public: 

		std::vector<float> vertex;
		std::vector<float> color;
		GLuint VAO, vertexVBO, colorVAO, colorVBO;
		float minX, maxX, minY, maxY, minZ, maxZ;
		float mass = 100000000;
		float e = 1.0;

		void initModelBox(Model& modelObject) {

			minX = modelObject.minX;
			maxX = modelObject.maxX;
			minY = modelObject.minY;
			maxY = modelObject.maxY;
			minZ = modelObject.minZ;
			maxZ = modelObject.maxZ;

			vertex.push_back(minX);
			vertex.push_back(minY);
			vertex.push_back(minZ);

			vertex.push_back(minX);
			vertex.push_back(maxY);
			vertex.push_back(minZ);
			
			vertex.push_back(minX);
			vertex.push_back(maxY);
			vertex.push_back(minZ);
			
			vertex.push_back(maxX);
			vertex.push_back(maxY);
			vertex.push_back(minZ);
			
			vertex.push_back(maxX);
			vertex.push_back(maxY);
			vertex.push_back(minZ);

			vertex.push_back(maxX);
			vertex.push_back(minY);
			vertex.push_back(minZ);

			vertex.push_back(maxX);
			vertex.push_back(minY);
			vertex.push_back(minZ);

			vertex.push_back(minX);
			vertex.push_back(minY);
			vertex.push_back(minZ);


			vertex.push_back(minX);
			vertex.push_back(minY);
			vertex.push_back(maxZ);

			vertex.push_back(minX);
			vertex.push_back(maxY);
			vertex.push_back(maxZ);
			
			vertex.push_back(minX);
			vertex.push_back(maxY);
			vertex.push_back(maxZ);
			
			vertex.push_back(maxX);
			vertex.push_back(maxY);
			vertex.push_back(maxZ);
			
			vertex.push_back(maxX);
			vertex.push_back(maxY);
			vertex.push_back(maxZ);

			vertex.push_back(maxX);
			vertex.push_back(minY);
			vertex.push_back(maxZ);

			vertex.push_back(maxX);
			vertex.push_back(minY);
			vertex.push_back(maxZ);

			vertex.push_back(minX);
			vertex.push_back(minY);
			vertex.push_back(maxZ);


			vertex.push_back(minX);
			vertex.push_back(minY);
			vertex.push_back(minZ);

			vertex.push_back(minX);
			vertex.push_back(minY);
			vertex.push_back(maxZ);
			
			vertex.push_back(minX);
			vertex.push_back(maxY);
			vertex.push_back(minZ);
			
			vertex.push_back(minX);
			vertex.push_back(maxY);
			vertex.push_back(maxZ);
			
			vertex.push_back(maxX);
			vertex.push_back(maxY);
			vertex.push_back(minZ);

			vertex.push_back(maxX);
			vertex.push_back(maxY);
			vertex.push_back(maxZ);

			vertex.push_back(maxX);
			vertex.push_back(minY);
			vertex.push_back(minZ);

			vertex.push_back(maxX);
			vertex.push_back(minY);
			vertex.push_back(maxZ);

			for (int i = 0; i < vertex.size(); i++) {

				color.push_back(1.0);
				color.push_back(0.6);
				color.push_back(0.1);
				color.push_back(1.0);

			}

			glGenVertexArrays(1, &VAO);
    	    glBindVertexArray(VAO);

        	glGenBuffers(1, &vertexVBO);
        	glBindBuffer(GL_ARRAY_BUFFER, vertexVBO);
        	glBufferData(GL_ARRAY_BUFFER, vertex.size() * sizeof(float), vertex.data(), GL_STATIC_DRAW);
        	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (GLvoid*)0);
	        glEnableVertexAttribArray(0);

        	glGenBuffers(1, &colorVBO);
        	glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
        	glBufferData(GL_ARRAY_BUFFER, color.size() * sizeof(float), color.data(), GL_STATIC_DRAW);
       		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (GLvoid*)0);
     	   	glEnableVertexAttribArray(1);

        	glBindVertexArray(0);

		}


		void displayModelBox() {

			glBindVertexArray(VAO);
			glDrawArrays(GL_LINES, 0, 2);
			glDrawArrays(GL_LINES, 2, 2);
			glDrawArrays(GL_LINES, 4, 2);
			glDrawArrays(GL_LINES, 6, 2);
			glDrawArrays(GL_LINES, 8, 2);
			glDrawArrays(GL_LINES, 10, 2);
			glDrawArrays(GL_LINES, 12, 2);
			glDrawArrays(GL_LINES, 14, 2);
			glDrawArrays(GL_LINES, 16, 2);
			glDrawArrays(GL_LINES, 18, 2);
			glDrawArrays(GL_LINES, 20, 2);
			glDrawArrays(GL_LINES, 22, 2);
			glBindVertexArray(0);			

		}

};


class SimBox {

	public: 

		float mass = 1000000000;
		float e = 1.0;
		std::vector<float> vertex;
		std::vector<float> color;
		int increment = 3;
		GLuint vertexVAO, vertexVBO, colorVAO, colorVBO;
		float minX, maxX, minY, maxY, minZ, maxZ;


		void initSimBox(Model& modelObject, ModelBox& modelBox) {

			for (int i = 0; i < modelBox.vertex.size(); i++) {

				vertex.push_back(increment * modelBox.vertex[i]);
				color.push_back(1.0);
				color.push_back(1.0);
				color.push_back(0.0);
				color.push_back(0.0);

			}

			minX = increment * modelObject.minX;
			maxX = increment * modelObject.maxX;
			minY = increment * modelObject.minY;
			maxY = increment * modelObject.maxY;
			minZ = increment * modelObject.minZ;
			maxZ = increment * modelObject.maxZ;

			glGenVertexArrays(1, &vertexVAO);
    	    glBindVertexArray(vertexVAO);

        	glGenBuffers(1, &vertexVBO);
        	glBindBuffer(GL_ARRAY_BUFFER, vertexVBO);
        	glBufferData(GL_ARRAY_BUFFER, vertex.size() * sizeof(float), vertex.data(), GL_STATIC_DRAW);
        	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (GLvoid*)0);
	        glEnableVertexAttribArray(0);

        	glGenBuffers(1, &colorVBO);
        	glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
        	glBufferData(GL_ARRAY_BUFFER, color.size() * sizeof(float), color.data(), GL_STATIC_DRAW);
       		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (GLvoid*)0);
     	   	glEnableVertexAttribArray(1);

        	glBindVertexArray(0);

		}


		void displaySimBox() {

			glBindVertexArray(vertexVAO);
			glDrawArrays(GL_LINES, 0, 2);
			glDrawArrays(GL_LINES, 2, 2);
			glDrawArrays(GL_LINES, 4, 2);
			glDrawArrays(GL_LINES, 6, 2);
			glDrawArrays(GL_LINES, 8, 2);
			glDrawArrays(GL_LINES, 10, 2);
			glDrawArrays(GL_LINES, 12, 2);
			glDrawArrays(GL_LINES, 14, 2);
			glDrawArrays(GL_LINES, 16, 2);
			glDrawArrays(GL_LINES, 18, 2);
			glDrawArrays(GL_LINES, 20, 2);
			glDrawArrays(GL_LINES, 22, 2);
			glBindVertexArray(0);						

		}

};


class Sphere {

	public: 


		GLfloat vertexData[(stacks + 1) * (sectors + 1) * 3];
		std::vector<GLfloat> color;
		GLuint indexData[stacks][sectors * 2 + 2];
    	GLuint vertexVAO, vertexVBO, vertexEBO, colorVAO, colorVBO;  	  // Create VAO, VBO, and EBO
		float radius;
		float coordinatesX, coordinatesY, coordinatesZ;
		float velocityX, velocityY, velocityZ;		
		float accelerationX, accelerationY, accelerationZ;
		float deltaTime = 0.0, frameTime = 0.00;
		float mass;
		float e;
		float V_c;


		void init(GLuint& shaderProgram) {

			int count = 0;
		
			for (int i = 0; i <= stacks; i++) {

			    double stackAngle = M_PI / 2 - i * M_PI / stacks;
    			double xy = radius * cos(stackAngle);
    			double z = (radius * sin(stackAngle));

    			for (int j = 0; j <= sectors; ++j) {

			        double sectorAngle = j * 2 * M_PI / sectors;
        			GLfloat x = (xy * cos(sectorAngle));
        			GLfloat y = (xy * sin(sectorAngle));

        			vertexData[count++] = x + coordinatesX;
        			vertexData[count++] = y + coordinatesY;
        			vertexData[count++] = z + coordinatesZ;

					color.push_back(1.0);
					color.push_back(0.0);
					color.push_back(0.0);

    			}

			}


			for (int i = 0; i < stacks; ++i) {

			    for (int j = 0; j <= sectors; ++j) {

					int first = i * (sectors + 1) + j;
        			int second = first + sectors + 1;

        			indexData[i][j * 2] = second;
        			indexData[i][j * 2 + 1] = first;
    			
				}

    			// Close the triangle strip loop for each stack
   	 			indexData[i][sectors * 2] = indexData[i][sectors * 2 + 1] = indexData[i][0];
			
			}

    		glGenVertexArrays(1, &vertexVAO);
    		glGenBuffers(1, &vertexEBO);

    		glBindVertexArray(vertexVAO);

    		glGenBuffers(1, &vertexVBO);
    		glBindBuffer(GL_ARRAY_BUFFER, vertexVBO);
    		glBufferData(GL_ARRAY_BUFFER, sizeof(vertexData), vertexData, GL_STATIC_DRAW);
    		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vertexEBO);
    		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indexData), indexData, GL_STATIC_DRAW);
    		glEnableVertexAttribArray(0);
    		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);

			glGenBuffers(1, &colorVBO);
        	glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
        	glBufferData(GL_ARRAY_BUFFER, color.size() * sizeof(float), color.data(), GL_STATIC_DRAW);
       		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (GLvoid*)0);
     	   	glEnableVertexAttribArray(1);

    		glBindBuffer(GL_ARRAY_BUFFER, 0);
    		glBindVertexArray(0);		

		}


		void borderCollision(SimBox& simBox) {

			if (coordinatesX + radius > simBox.maxX || coordinatesX - radius < simBox.minX) {

				velocityX =	((velocityX) * ((mass) - (simBox.mass * simBox.e))) / (mass + simBox.mass);

			}
					
			if (coordinatesY + radius > simBox.maxY || coordinatesY - radius < simBox.minY) {

				velocityY =	((velocityY) * ((mass) - (simBox.mass * simBox.e))) / (mass + simBox.mass);

			}

			else if (coordinatesZ + radius > simBox.maxZ || coordinatesZ - radius < simBox.minZ) {

				velocityZ =	((velocityZ) * ((mass) - (simBox.mass * simBox.e))) / (mass + simBox.mass);

			}

		}

		void objectCollision(std::vector<Sphere>& sphereObjects) {
			
			for (int i = 0; i < sphereObjects.size(); i++) {

				for (int j = i + 1; j < sphereObjects.size(); j++) {

					float dx = sphereObjects[j].coordinatesX - sphereObjects[i].coordinatesX - (sphereObjects[j].radius - sphereObjects[i].radius);
					float dy = sphereObjects[j].coordinatesY - sphereObjects[i].coordinatesY - (sphereObjects[j].radius - sphereObjects[i].radius);
					float dz = sphereObjects[j].coordinatesZ - sphereObjects[i].coordinatesZ - (sphereObjects[j].radius - sphereObjects[i].radius);
					float distance = sqrt((dx * dx) + (dy * dy) + (dz * dz));

					float limit = sphereObjects[i].radius + sphereObjects[j].radius;

					if (distance <= limit) {

						/* float velX = objects[i].velocityX; */

						/* velocityX = ((objects[i].mass * objects[i].velocityX) + (objects[j].mass * objects[j].velocityX) - (objects[j].mass * objects[i].e * (objects[i].velocityX - objects[j].velocityX))) / (objects[i].mass + objects[j].mass); */
						/* objects[j].velocityX = ((objects[i].mass * velX) + (objects[j].mass * objects[j].velocityX) + (objects[i].mass * objects[j].e * (velX - objects[j].velocityX))) / (objects[i].mass + objects[j].mass); */

						float normalX = dx / distance;
			            float normalY = dy / distance;
			            float normalZ = dz / distance;

                		float relativeVelocityX = sphereObjects[j].velocityX - sphereObjects[i].velocityX;
               		 	float relativeVelocityY = sphereObjects[j].velocityY - sphereObjects[i].velocityY;
               		 	float relativeVelocityZ = sphereObjects[j].velocityZ - sphereObjects[i].velocityZ;
                		float dotProduct = (relativeVelocityX * normalX) + (relativeVelocityY * normalY) + (relativeVelocityZ * normalZ);

	           		    sphereObjects[i].velocityX += normalX * dotProduct;
    	           		sphereObjects[i].velocityY += normalY * dotProduct;
    	           		sphereObjects[i].velocityZ += normalZ * dotProduct;

        	       		sphereObjects[j].velocityX -= normalX * dotProduct;
            	   		sphereObjects[j].velocityY -= normalY * dotProduct;
	           	   		sphereObjects[j].velocityZ -= normalZ * dotProduct;

						/* float impulse = ((1 + sphereObjects[i].e) * dotProduct) / (sphereObjects[i].mass + sphereObjects[j].mass); */

	          		    /* sphereObjects[i].velocityX -= (impulse * sphereObjects[i].mass * normalX) * V_c; */
    	           		/* sphereObjects[i].velocityY -= (impulse * sphereObjects[i].mass * normalY) * V_c; */
    	           		/* sphereObjects[i].velocityZ -= (impulse * sphereObjects[i].mass * normalZ) * V_c; */

        	       		/* sphereObjects[j].velocityX += (impulse * sphereObjects[i].mass * normalX) * V_c; */
            	   		/* sphereObjects[j].velocityY += (impulse * sphereObjects[i].mass * normalY) * V_c; */ 
	           	   		/* sphereObjects[j].velocityZ += (impulse * sphereObjects[i].mass * normalZ) * V_c; */

					}	

				}

			}

		}


		void modelCollision(std::vector<Sphere>& sphereObject, Model& modelObject) {

			for (int i = 0; i < modelObject.vertices.size(); i += 9) {

				std::vector<float> Xvector = order(modelObject.vertices[i], modelObject.vertices[i + 3], modelObject.vertices[i + 6]);	
				std::vector<float> Yvector = order(modelObject.vertices[i + 1], modelObject.vertices[i + 4], modelObject.vertices[i + 7]);
				std::vector<float> Zvector = order(modelObject.vertices[i + 3], modelObject.vertices[i + 5], modelObject.vertices[i + 8]);

				if (coordinatesX > Xvector[0] && coordinatesX < Xvector[2] && coordinatesY > Yvector[0] && coordinatesY < Yvector[2] && coordinatesZ > Zvector[0] && coordinatesZ < Zvector[2]) {

					float dx = -radius;
					float dy = -radius;
					float dz = -radius;
					float distance = sqrt((dx * dx) + (dy * dy) + (dz * dz));
					
					float normalX = dx / distance;
		    	    float normalY = dy / distance;
		        	float normalZ = dz / distance;

       			  	float dotProduct = (-velocityX * normalX) + (-velocityY * normalY) + (-velocityZ * normalZ);
					/* float impulse = ((1 + e) * dotProduct) / (mass + modelObject.mass); */
	
					/* velocityX -= (impulse * mass * normalX) * V_c; */
					/* velocityY -= (impulse * mass * normalY) * V_c; */
					/* velocityZ -= (impulse * mass * normalZ) * V_c; */

          			velocityX += normalX * dotProduct;
	           		velocityY += normalY * dotProduct;
    	       		velocityZ += normalZ * dotProduct;

				}

			}

		}


		void update(std::vector<Sphere>& sphereObjects, SimBox& simBox, Model& modelObject) {

			GLfloat totalTime = glfwGetTime();
			deltaTime = totalTime - frameTime;
			frameTime = totalTime;

			velocityX += accelerationX * deltaTime;
			velocityY += accelerationY * deltaTime;
			velocityZ += accelerationZ * deltaTime;

			coordinatesX = (coordinatesX) + (velocityX * deltaTime) + ((0.5) * (accelerationX) * (deltaTime * deltaTime));	
			coordinatesY = (coordinatesY) + (velocityY * deltaTime) + ((0.5) * (accelerationY) * (deltaTime * deltaTime));
			coordinatesZ = (coordinatesZ) + (velocityZ * deltaTime) + ((0.5) * (accelerationZ) * (deltaTime * deltaTime));

			int count = 0;
		
			for (int i = 0; i <= stacks; i++) {

			    double stackAngle = M_PI / 2 - i * M_PI / stacks;
    			double xy = radius * cos(stackAngle);
    			double z = (radius * sin(stackAngle));

    			for (int j = 0; j <= sectors; ++j) {

			        double sectorAngle = j * 2 * M_PI / sectors;
        			GLfloat x = (xy * cos(sectorAngle));
        			GLfloat y = (xy * sin(sectorAngle));

        			vertexData[count++] = x + coordinatesX;
        			vertexData[count++] = y + coordinatesY;
        			vertexData[count++] = z + coordinatesZ;

    			}

			}

			borderCollision(simBox);

			objectCollision(sphereObjects);

			/* modelCollision(sphereObjects, modelObject); */

			// Update vertex buffer data
    		glBindBuffer(GL_ARRAY_BUFFER, vertexVBO);
    		glBufferData(GL_ARRAY_BUFFER, sizeof(vertexData), vertexData, GL_STATIC_DRAW);
    		glBindBuffer(GL_ARRAY_BUFFER, 0);
		
		}

		void display(GLuint& shaderProgram) {
	
			// Draw the sphere
    		glBindVertexArray(vertexVAO);
    		
			for (int i = 0; i < stacks; ++i) {

        		glDrawElements(GL_TRIANGLE_STRIP, sectors * 2, GL_UNSIGNED_INT, (void*)(i * sectors * 2 * sizeof(GLuint)));
    		
			}
    		
			glBindVertexArray(0);

		}	

		void cleanUp() {

			// Cleanup
    		glDeleteVertexArrays(1, &vertexVAO);
    		glDeleteBuffers(1, &vertexVBO);
    		glDeleteBuffers(1, &vertexEBO);

		}

};

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {

	glViewport(0, 0, width, height);

}


void processInput(GLFWwindow* window, glm::vec3& cameraPos, glm::vec3& cameraFront, glm::vec3& cameraRight, glm::vec3& cameraUp, float& cameraYaw, float& cameraPitch) {

	// Calculate delta time to smooth out camera movement
    static float lastFrame = 0.0f;
    float currentFrame = glfwGetTime();
    float deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;

    float cameraSpeed = 2.5f * deltaTime;

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        cameraPos += cameraSpeed * cameraFront;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        cameraPos -= cameraSpeed * cameraFront;
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;

    // Handle mouse movement
    double mouseX, mouseY;
    glfwGetCursorPos(window, &mouseX, &mouseY);

    // Calculate the change in mouse position
    double deltaX = mouseX - lastMouseX;
    double deltaY = lastMouseY - mouseY; // Note the inversion since Y-coordinates increase upwards

    // Update last known mouse position
    lastMouseX = mouseX;
    lastMouseY = mouseY;

    // Sensitivity factor for mouse movement
    float sensitivity = 0.1f;

	    // Check if the mouse button is currently held down
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
        
		if (!mouseButtonPressed) {

            // Initialize the last known position to the current position to avoid sudden jumps
            lastMouseX = mouseX;
			lastMouseY = mouseY;
    
		}

		else {
		
			mouseButtonPressed = true;
	
			cameraYaw += deltaY * sensitivity;
    		cameraPitch += deltaX * sensitivity;

	    	// Limit the camera pitch to avoid flipping
	    	if (cameraPitch > 89.0f)
    	    	cameraPitch = 89.0f;
    			if (cameraPitch < -89.0f)
        		cameraPitch = -89.0f;

    		// Calculate the new front vector
    		glm::vec3 front;
    		front.x = cos(glm::radians(cameraYaw)) * cos(glm::radians(cameraPitch));
    		front.y = sin(glm::radians(cameraPitch));
    		front.z = sin(glm::radians(cameraYaw)) * cos(glm::radians(cameraPitch));
    		cameraFront = glm::normalize(front);
	
		}
	}
	
	else {
    
		mouseButtonPressed = false;
    
	}
	
}


void init(Model& modelObject, std::vector<Sphere>& sphereObjects, Axis& axis, ModelBox& modelBox, SimBox& simBox, GLuint& shaderProgram) {

	axis.initAxis();
	modelObject.initModel();
	modelBox.initModelBox(modelObject);
	simBox.initSimBox(modelObject, modelBox);

	float radius = 0.1;

	int circlesZ = floor((simBox.maxZ - (simBox.minZ)) / (2 * radius));
	int circlesX = floor((simBox.maxX - (simBox.minX)) / (2 * radius));

	float incrementX = (simBox.maxX - simBox.minX) / (circlesX + 2);
	float incrementZ = (simBox.maxZ - simBox.minZ) / (circlesZ + 2);

	float positionX = 0.0;
	float positionY = simBox.maxY;
	float positionZ = 0.0;

	for (int i = 0; i < numberObjects; i++) {

		sphereObjects.push_back(Sphere());

		sphereObjects[i].radius = radius;

		sphereObjects[i].coordinatesX = positionX;
		sphereObjects[i].coordinatesY = positionY;
		sphereObjects[i].coordinatesZ = positionZ;

		sphereObjects[i].velocityX = 1.0;
		sphereObjects[i].velocityY = -5.0;
		sphereObjects[i].velocityZ = -1.0;

		sphereObjects[i].accelerationX = 0.0;
		sphereObjects[i].accelerationY = 0.0;
		sphereObjects[i].accelerationZ = 0.0;

		sphereObjects[i].mass = 0.05;
		sphereObjects[i].e = 1.0;
		sphereObjects[i].V_c = 1.0;

		positionX += incrementX;
		positionZ += incrementZ;

		sphereObjects[i].init(shaderProgram);

		/* std::cout << "hi" << i << '\n'; */

	}

}


void update(std::vector<Sphere>& sphereObjects, SimBox& simBox, Model& modelObject) {

	for (int i = 0; i < sphereObjects.size(); i++) {

		sphereObjects[i].update(sphereObjects, simBox, modelObject);

	}

}


void display(Model& modelObject, std::vector<Sphere>& sphereObjects, Axis& axis, ModelBox& modelBox, SimBox& simBox, glm::vec3& cameraPos, glm::vec3& cameraFront, glm::vec3& cameraUp, GLFWwindow* window, GLuint& shaderProgram) {

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(shaderProgram);

	glm::mat4 model = glm::mat4(1.0f);
	glm::mat4 view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
    glm::mat4 projection = glm::perspective(glm::radians(45.0f), 800.0f / 600.0f, 0.1f, 100.0f);

    GLuint modelLoc = glGetUniformLocation(shaderProgram, "model");
    GLuint viewLoc = glGetUniformLocation(shaderProgram, "view");
    GLuint projectionLoc = glGetUniformLocation(shaderProgram, "projection");

    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));

	axis.displayAxis();

	/* modelObject.displayModel(); */

	for (int i = 0; i < sphereObjects.size(); i++) {

		sphereObjects[i].display(shaderProgram);

		/* std::cout << "hi" << i << '\n'; */

	}

	/* modelBox.displayModelBox(); */

	simBox.displaySimBox();
	
	glfwSwapBuffers(window);

}


int main() {

	Model modelObject;

	Axis axis;

	ModelBox modelBox;

	SimBox simBox;

	std::vector<Sphere> sphereObjects;

	modelObject.filename = "Leviathan.obj";


    if (!glfwInit()) {
    
		std::cerr << "GLFW initialization failed" << std::endl;
        return -1;
    
	}

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(800, 600, "OpenGL Leviathan", nullptr, nullptr);

	if (!window) {

		std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
		return -1;

    }

    glfwMakeContextCurrent(window);

    if (glewInit() != GLEW_OK) {
        
		std::cerr << "GLEW initialization failed" << std::endl;
        return -1;
    
	}

    glViewport(0, 0, 800, 600);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	
    GLuint vertexShader, fragmentShader;
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
    glCompileShader(vertexShader);

    // Check for compilation errors and handle them

    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
    glCompileShader(fragmentShader);

    // Check for compilation errors and handle them

    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    // Check for linking errors and handle them

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

	init(modelObject, sphereObjects, axis, modelBox, simBox, shaderProgram);

	// Define camera properties
	glm::vec3 cameraPos = glm::vec3(simBox.maxX, simBox.maxY, simBox.maxZ);
	glm::vec3 cameraFront = glm::vec3(-simBox.maxX, -simBox.maxY, -simBox.maxZ);
	glm::vec3 cameraUp = glm::vec3(0.0f, 0.0f, 1.0f);
	glm::vec3 cameraRight = glm::vec3(0.0f);
	float cameraYaw = -0.0f;   // Initialized to face along negative z-axis
	float cameraPitch = 0.0f;   // Initialized to zero


    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	 while (!glfwWindowShouldClose(window)) {

        processInput(window, cameraPos, cameraFront, cameraRight, cameraUp, cameraYaw, cameraPitch); // Handle camera movement
		update(sphereObjects, simBox, modelObject);
		display(modelObject, sphereObjects, axis, modelBox, simBox, cameraPos, cameraFront, cameraUp, window, shaderProgram);
        glfwPollEvents();

    }
	 
	modelObject.cleanUpModel();

    glDeleteProgram(shaderProgram);

    glfwTerminate();
    return 0;

}
