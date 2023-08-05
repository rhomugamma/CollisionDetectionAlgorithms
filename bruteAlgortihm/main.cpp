#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <math.h>
#include <vector>

const float PI = 3.14159265358979323;

const float g = 9.80665;


class box {

	public:

		float mass = 10000000000;

		GLfloat verticesbox[16] = {

			-1.0, -1.0,
			 1.0, -1.0,

			 1.0, -1.0,
			 1.0,  1.0,

			 1.0, 1.0,
			-1.0, 1.0,

			-1.0, 1.0,
			-1.0, -1.0

		};


		GLfloat colorbox[24] {

			0.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.0f,

			0.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.0f,

			0.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.0f,

			0.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.0f,

		};

		GLuint VAO;

		void renderBoxObject() {

	        glGenVertexArrays(1, &VAO);
    	    glBindVertexArray(VAO);

	        GLuint vertexVBO;
        	glGenBuffers(1, &vertexVBO);
        	glBindBuffer(GL_ARRAY_BUFFER, vertexVBO);
        	glBufferData(GL_ARRAY_BUFFER, sizeof(verticesbox), verticesbox, GL_STATIC_DRAW);
        	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), (GLvoid*)0);
	        glEnableVertexAttribArray(0);

        	GLuint colorVBO;
        	glGenBuffers(1, &colorVBO);
        	glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
        	glBufferData(GL_ARRAY_BUFFER, sizeof(colorbox), colorbox, GL_STATIC_DRAW);
       		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
     	   	glEnableVertexAttribArray(1);

        	glBindVertexArray(0);
    	
		}

	    void cleanup() {

    	    glDeleteVertexArrays(1, &VAO);
    
		}

};

class object {

	public:

		//Physical properties
		float radius;
		float mass;
		float weight = (mass) * (g);
		float color1;
		float color2;
		float color3;
		float e;

		//Rendering properties
		int iterations = 15;
		float alpha = (360 / 15) * (2 * PI / 360);
		float beta = alpha;
		/* const static int verticesArraySize = 6 * iterations; */
		/* const static int colorsArraySize = 9 * iterations; */
		GLfloat vertices[90];													
		GLfloat colors[135];
		GLuint VAO;

		//Goberning characteristics
		float coordinatesX;
		float coordinatesY;

		float velocityX;
		float velocityY;
		float velocityModule = sqrt((velocityX * velocityX) + (velocityY * velocityY));

		float accelerationX;
		float accelerationY;
		float accelerationModule = sqrt((accelerationX * accelerationX) + (accelerationY * accelerationY));

		float angularVelocityZ;
		float angularVelocityModule = sqrt(angularVelocityZ * angularVelocityZ);

		float angularAccelerationZ;
		float angularAcceleration = sqrt(angularAccelerationZ * angularAccelerationZ);

		float linearMomentumX = (mass) * (velocityX);
	   	float linearMomentumY = (mass) * (velocityY);
		float linearMomentumModule = sqrt((linearMomentumX * linearMomentumX) + (linearMomentumY * linearMomentumY));
 		float momentInertia = (mass * (radius * radius)) / (2);

		float kinetickEnergy = ((0.5) + (mass) + (velocityModule * velocityModule)) + ((0.5) + (momentInertia) + (angularVelocityZ * angularVelocityZ)); 
		float potentialEnergy = (mass) * (g) * (coordinatesY + 0.5);

		//Time characteristics
		float deltaTime;
		float frameTime;

		void renderSphereObject() {

			GLfloat totalTime = glfwGetTime();
			deltaTime = totalTime - frameTime;
			frameTime = totalTime;


			for (int i = 0; i < 6 * iterations; i += 6) {

				vertices[i] = coordinatesX;
        		vertices[i + 1] = coordinatesY;

		        vertices[i + 2] = radius * cos(alpha) + coordinatesX;
 		       	vertices[i + 3] = radius * sin(alpha) + coordinatesY;

 		       	vertices[i + 4] = radius * cos(alpha + beta) + coordinatesX;
        		vertices[i + 5] = radius * sin(alpha + beta) + coordinatesY;

				alpha += beta;

			};

    		for (int i = 0; i < 9 * iterations; i += 3) {

		        colors[i] = color1;
       			colors[i + 1] = color2;
     		    colors[i + 2] = color3;
    
			};

		}


		void updateObjectPosition(std::vector<object>& objects, box box1) {

		GLfloat totalTime = glfwGetTime();
		deltaTime = totalTime - frameTime;
		frameTime = totalTime;

		borderCollision(box1);

		objectCollision(objects);

		coordinatesX = (coordinatesX) + (velocityX * deltaTime) + ((1/2) * (accelerationX) * (deltaTime * deltaTime));	
		coordinatesY = (coordinatesY) + (velocityY * deltaTime) + ((1/2) * (accelerationY) * (deltaTime * deltaTime));

		velocityX = (velocityX) + ((accelerationX) * (deltaTime));
		velocityY = (velocityY) + ((accelerationY) * (deltaTime));

		
		for (int i = 0; i < 6 * iterations; i += 6) {

				vertices[i] = coordinatesX;
	    	    vertices[i + 1] = coordinatesY;

				vertices[i + 2] = radius * cos(alpha) + coordinatesX;
    		    vertices[i + 3] = radius * sin(alpha) + coordinatesY;

 				vertices[i + 4] = radius * cos(alpha + beta) + coordinatesX;
	        	vertices[i + 5] = radius * sin(alpha + beta) + coordinatesY;

				alpha += beta;

			};
 
		}


		void renderObject() {

	        glGenVertexArrays(1, &VAO);
    	    glBindVertexArray(VAO);

	        GLuint vertexVBO;
        	glGenBuffers(1, &vertexVBO);
        	glBindBuffer(GL_ARRAY_BUFFER, vertexVBO);
        	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
        	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), (GLvoid*)0);
	        glEnableVertexAttribArray(0);

        	GLuint colorVBO;
        	glGenBuffers(1, &colorVBO);
        	glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
        	glBufferData(GL_ARRAY_BUFFER, sizeof(colors), colors, GL_STATIC_DRAW);
       		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
     	   	glEnableVertexAttribArray(1);

        	glBindVertexArray(0);
    	
		}


		void borderCollision(box box1) {


			if (coordinatesX + radius > box1.verticesbox[2] || coordinatesX - radius < box1.verticesbox[0]) {

				/* velocityX = -velocityX; */
				velocityX =	((velocityX) * ((mass) - (box1.mass * e))) / (mass + box1.mass);

			}
				
				
			else if (coordinatesY + radius > box1.verticesbox[7] || coordinatesY - radius < box1.verticesbox[1]) {

				/* velocityY = -velocityY; */	
				velocityY =	((velocityY) * ((mass) - (box1.mass * e))) / (mass + box1.mass);

			}

		}	


		void objectCollision(std::vector<object>& objects) {

			for (int i = 0; i < objects.size(); i++) {

				for (int j = i + 1; j < objects.size(); j++) {

						float dx = objects[j].coordinatesX - objects[i].coordinatesX - (objects[j].radius - objects[i].radius);
						float dy = objects[j].coordinatesY - objects[i].coordinatesY - (objects[j].radius - objects[i].radius);
						float distance = sqrt((dx * dx) + (dy * dy));

						float limit = objects[i].radius + objects[j].radius;

						if (distance <= limit) {

							/* float velX = objects[i].velocityX; */

							/* velocityX = ((objects[i].mass * objects[i].velocityX) + (objects[j].mass * objects[j].velocityX) - (objects[j].mass * objects[i].e * (objects[i].velocityX - objects[j].velocityX))) / (objects[i].mass + objects[j].mass); */
							/* objects[j].velocityX = ((objects[i].mass * velX) + (objects[j].mass * objects[j].velocityX) + (objects[i].mass * objects[j].e * (velX - objects[j].velocityX))) / (objects[i].mass + objects[j].mass); */

							float normalX = dx / distance;
			                float normalY = dy / distance;

                			float relativeVelocityX = objects[j].velocityX - objects[i].velocityX;
               		 		float relativeVelocityY = objects[j].velocityY - objects[i].velocityY;
                			float dotProduct = (relativeVelocityX * normalX) + (relativeVelocityY * normalY);

	            		    objects[i].velocityX += normalX * dotProduct;
    	            		objects[i].velocityY += normalY * dotProduct;
        	        		objects[j].velocityX -= normalX * dotProduct;
            	    		objects[j].velocityY -= normalY * dotProduct;

						}	

				}


			}

		}

	    void cleanup() {

    	    glDeleteVertexArrays(1, &VAO);
    
		}

};

class sample {
	
	public:

		float radius = 0.15;
		float mass = 10000000000;
		float weight = (mass) * (g);

		//Rendering properties
		int iterations = 60;
		float alpha = (360 / 60) * (2 * PI / 360);
		float beta = alpha;
		/* const static int verticesArraySize = 6 * iterations; */
		/* const static int colorsArraySize = 9 * iterations; */
		GLfloat vertices[360];													
		GLfloat colors[540];
		GLuint VAO;

		//Goberning characteristics
		float coordinatesX = -0.75;
		float coordinatesY = -0.00;

		//Time characteristics
		float deltaTime;
		float frameTime;

		void renderSampleObject() {

			GLfloat totalTime = glfwGetTime();
			deltaTime = totalTime - frameTime;
			frameTime = totalTime;


			for (int i = 0; i < 6 * iterations; i += 6) {

				vertices[i] = coordinatesX;
        		vertices[i + 1] = coordinatesY;

		        vertices[i + 2] = radius * cos(alpha) + coordinatesX;
 		       	vertices[i + 3] = radius * sin(alpha) + coordinatesY;

 		       	vertices[i + 4] = radius * cos(alpha + beta) + coordinatesX;
        		vertices[i + 5] = radius * sin(alpha + beta) + coordinatesY;

				alpha += beta;

			};

    		for (int i = 0; i < 9 * iterations; i += 3) {

		        colors[i] = 1.0;
       			colors[i + 1] = 1.0;
     		    colors[i + 2] = 1.0;
    
			};

		}
	
		void renderObject() {

	        glGenVertexArrays(1, &VAO);
    	    glBindVertexArray(VAO);

	        GLuint vertexVBO;
        	glGenBuffers(1, &vertexVBO);
        	glBindBuffer(GL_ARRAY_BUFFER, vertexVBO);
        	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
        	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), (GLvoid*)0);
	        glEnableVertexAttribArray(0);

        	GLuint colorVBO;
        	glGenBuffers(1, &colorVBO);
        	glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
        	glBufferData(GL_ARRAY_BUFFER, sizeof(colors), colors, GL_STATIC_DRAW);
       		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
     	   	glEnableVertexAttribArray(1);

        	glBindVertexArray(0);
    	
		}

	    void cleanup() {

    	    glDeleteVertexArrays(1, &VAO);
    
		}

};

const GLchar* vertexShaderSource = R"(

    #version 450 core

    layout (location = 0) in vec3 position;
	layout (location = 1) in vec3 color;

	out vec3 vertexColor;

    void main() {

		vertexColor = color;
        gl_Position = vec4(position.x, position.y, 0.0, 1.0);

    }

)";

const GLchar* fragmentShaderSource = R"(

    #version 450 core

	in vec3 vertexColor;
	out vec4 fragColor;

    void main() {

        fragColor = vec4(vertexColor, 1.0f);

    }

)";


void init(std::vector<object>& objects, sample& sample1);

void display(std::vector<object>& objects, box box1, sample& sample1, GLFWwindow* window, GLuint& shaderProgram);

void framebufferSizeCallback(GLFWwindow* window, int width, int height);

void rendersphere(object& obt);

void renderbox(box box1);

void rendersample(sample& sample1);

void updateObjectPosition(object& obt);


int main() {

	std::vector<object> objects;

	box box1;

	sample sample1;

	GLuint shaderProgram;


    if (!glfwInit()) {

        std::cout << "GLFW initialization failed!" << std::endl;
        return -1;
    
	}

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    GLFWwindow* window = glfwCreateWindow(3840, 2160, "Dynamic Window", nullptr, nullptr);

    if (!window) {

        std::cout << "GLFW window creation failed!" << std::endl;
        glfwTerminate();
        return -1;

    }

    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;

    if (glewInit() != GLEW_OK) {
		 
        std::cout << "GLEW initialization failed!" << std::endl;
        glfwDestroyWindow(window);
        glfwTerminate();
        return -1;
    
	}

    glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);

    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
    glCompileShader(vertexShader);

    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
    glCompileShader(fragmentShader);

    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    init(objects, sample1);

    while (!glfwWindowShouldClose(window)) {

        glfwPollEvents();
		display(objects, box1, sample1, window, shaderProgram);
    
	}

    glDeleteProgram(shaderProgram);

    glfwTerminate();

    return 0;
}

void init(std::vector<object>& objects, sample& sample1) {

	float xposition = -0.0;
	float yposition = -0.0;
	float yprime = yposition;
	float numberObjects = 1;
	float height = 3.94;
	float increase = height / numberObjects;
	

	for (int i = 0; i < numberObjects; i++) {

		objects.push_back(object());

		objects[i].radius = 0.01;						//0.001538
		objects[i].mass = 4.809 * pow(10, -23);
		objects[i].color1 = 1.0;
		objects[i].color2 = 0.0;
		objects[i].color3 = 0.0;
		objects[i].e = 0.85;

		objects[i].coordinatesX = xposition;
		objects[i].coordinatesY = yposition;
	
		objects[i].velocityX = 0.0;
		objects[i].velocityY = 0.0;

		objects[i].accelerationX = 0.0;
		objects[i].accelerationY = -1.62;

		objects[i].angularVelocityZ = 0.0;

		objects[i].angularAccelerationZ = 0.0;

		objects[i].deltaTime = 0.0;
		objects[i].frameTime = 0.0;

		objects[i].renderSphereObject();
					  
		if ((yposition + increase) < 1.0) {
			
			yposition += increase;

		}

		else {

			yposition = yprime;

		}

	}

	sample1.renderSampleObject();
	
}

void display(std::vector<object>& objects, box box1, sample& sample1, GLFWwindow* window, GLuint& shaderProgram) {


	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT);


    glUseProgram(shaderProgram);


	for (int i = 0; i < objects.size(); i++) {

		objects[i].updateObjectPosition(objects, box1);
		rendersphere(objects[i]);

	}

	renderbox(box1);

	rendersample(sample1);

	glfwSwapBuffers(window);

}

void framebufferSizeCallback(GLFWwindow* window, int width, int height) {

    glViewport(0, 0, width, height);

}

void rendersphere(object& obt) {

	obt.renderObject();
	glBindVertexArray(obt.VAO);
	glDrawArrays(GL_TRIANGLES, 0, 3 * obt.iterations);
	glBindVertexArray(0);
 
}

void renderbox(box box1) {

	box1.renderBoxObject();
	glBindVertexArray(box1.VAO);
	glDrawArrays(GL_LINES, 0, 4);
	glDrawArrays(GL_LINES, 4, 4);
	glDrawArrays(GL_LINES, 8, 4);
	glDrawArrays(GL_LINES, 12, 4);
	glBindVertexArray(0);

}

void rendersample(sample& sample1) {

	sample1.renderObject();
	glBindVertexArray(sample1.VAO);
	glDrawArrays(GL_TRIANGLES, 0, 3 * sample1.iterations);
	glBindVertexArray(0);

}
