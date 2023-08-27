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
bool mouseButtonPressed = false;


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


class Model {

	public:

		std::string filename;	
		std::vector<float> vertices;
    	std::vector<unsigned int> indices;
		std::vector<float> color;
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

		void initModelBox(Model& modelObject) {

			vertex.push_back(modelObject.minX);
			vertex.push_back(modelObject.minY);
			vertex.push_back(modelObject.minZ);

			vertex.push_back(modelObject.minX);
			vertex.push_back(modelObject.maxY);
			vertex.push_back(modelObject.minZ);
			
			vertex.push_back(modelObject.minX);
			vertex.push_back(modelObject.maxY);
			vertex.push_back(modelObject.minZ);
			
			vertex.push_back(modelObject.maxX);
			vertex.push_back(modelObject.maxY);
			vertex.push_back(modelObject.minZ);
			
			vertex.push_back(modelObject.maxX);
			vertex.push_back(modelObject.maxY);
			vertex.push_back(modelObject.minZ);

			vertex.push_back(modelObject.maxX);
			vertex.push_back(modelObject.minY);
			vertex.push_back(modelObject.minZ);

			vertex.push_back(modelObject.maxX);
			vertex.push_back(modelObject.minY);
			vertex.push_back(modelObject.minZ);

			vertex.push_back(modelObject.minX);
			vertex.push_back(modelObject.minY);
			vertex.push_back(modelObject.minZ);


			vertex.push_back(modelObject.minX);
			vertex.push_back(modelObject.minY);
			vertex.push_back(modelObject.maxZ);

			vertex.push_back(modelObject.minX);
			vertex.push_back(modelObject.maxY);
			vertex.push_back(modelObject.maxZ);
			
			vertex.push_back(modelObject.minX);
			vertex.push_back(modelObject.maxY);
			vertex.push_back(modelObject.maxZ);
			
			vertex.push_back(modelObject.maxX);
			vertex.push_back(modelObject.maxY);
			vertex.push_back(modelObject.maxZ);
			
			vertex.push_back(modelObject.maxX);
			vertex.push_back(modelObject.maxY);
			vertex.push_back(modelObject.maxZ);

			vertex.push_back(modelObject.maxX);
			vertex.push_back(modelObject.minY);
			vertex.push_back(modelObject.maxZ);

			vertex.push_back(modelObject.maxX);
			vertex.push_back(modelObject.minY);
			vertex.push_back(modelObject.maxZ);

			vertex.push_back(modelObject.minX);
			vertex.push_back(modelObject.minY);
			vertex.push_back(modelObject.maxZ);


			vertex.push_back(modelObject.minX);
			vertex.push_back(modelObject.minY);
			vertex.push_back(modelObject.minZ);

			vertex.push_back(modelObject.minX);
			vertex.push_back(modelObject.minY);
			vertex.push_back(modelObject.maxZ);
			
			vertex.push_back(modelObject.minX);
			vertex.push_back(modelObject.maxY);
			vertex.push_back(modelObject.minZ);
			
			vertex.push_back(modelObject.minX);
			vertex.push_back(modelObject.maxY);
			vertex.push_back(modelObject.maxZ);
			
			vertex.push_back(modelObject.maxX);
			vertex.push_back(modelObject.maxY);
			vertex.push_back(modelObject.minZ);

			vertex.push_back(modelObject.maxX);
			vertex.push_back(modelObject.maxY);
			vertex.push_back(modelObject.maxZ);

			vertex.push_back(modelObject.maxX);
			vertex.push_back(modelObject.minY);
			vertex.push_back(modelObject.minZ);

			vertex.push_back(modelObject.maxX);
			vertex.push_back(modelObject.minY);
			vertex.push_back(modelObject.maxZ);

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


void init(Model& modelObject, Axis& axis, ModelBox& modelBox, SimBox& simBox) {

	axis.initAxis();
	modelObject.initModel();
	modelBox.initModelBox(modelObject);
	simBox.initSimBox(modelObject, modelBox);

}


void display(Model& modelObject, Axis& axis, ModelBox& modelBox, SimBox& simBox, glm::vec3& cameraPos, glm::vec3& cameraFront, glm::vec3& cameraUp, GLFWwindow* window, GLuint& shaderProgram) {

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

	modelObject.displayModel();

	modelBox.displayModelBox();

	simBox.displaySimBox();
	
	glfwSwapBuffers(window);

}


int main() {

	Model modelObject;

	Axis axis;

	ModelBox modelBox;

	SimBox simBox;

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

	init(modelObject, axis, modelBox, simBox);

	// Define camera properties
	glm::vec3 cameraPos = glm::vec3(simBox.maxX, simBox.maxY, simBox.maxZ);
	glm::vec3 cameraFront = glm::vec3(-simBox.maxX, -simBox.maxY, -simBox.maxZ);
	glm::vec3 cameraUp = glm::vec3(0.0f, 0.0f, 1.0f);
	glm::vec3 cameraRight = glm::vec3(0.0f);
	float cameraYaw = -0.0f;   // Initialized to face along negative z-axis
	float cameraPitch = 0.0f;   // Initialized to zero

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

    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	 while (!glfwWindowShouldClose(window)) {

        processInput(window, cameraPos, cameraFront, cameraRight, cameraUp, cameraYaw, cameraPitch); // Handle camera movement
		display(modelObject, axis, modelBox, simBox, cameraPos, cameraFront, cameraUp, window, shaderProgram);
        glfwPollEvents();

    }
	 
	modelObject.cleanUpModel();

    glDeleteProgram(shaderProgram);

    glfwTerminate();
    return 0;

}
