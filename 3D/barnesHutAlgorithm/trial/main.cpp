#include <iostream>
#include <vector>
#include <cmath>

// Define a structure to represent a sphere
struct Sphere {
    double x, y, z;   // Center coordinates
    double radius;    // Radius
};

// Define a quadtree node for Barnes-Hut
struct QuadTreeNode {
    double mass;
    double cx, cy, cz;  // Center of mass coordinates
    bool isLeaf;
    Sphere* sphere;     // Only used for leaf nodes
    QuadTreeNode* children[8]; // Subdivision into octants
};

// Function to check if two spheres collide
bool spheresCollide(const Sphere& s1, const Sphere& s2) {
    double dx = s1.x - s2.x;
    double dy = s1.y - s2.y;
    double dz = s1.z - s2.z;
    double distanceSquared = dx * dx + dy * dy + dz * dz;
    double minDistance = s1.radius + s2.radius;
    return distanceSquared <= (minDistance * minDistance);
}

// Function to insert a sphere into the Barnes-Hut quadtree
void insert(QuadTreeNode* node, Sphere* sphere, double x0, double y0, double z0, double size) {
    if (node->isLeaf) {
        // Node is a leaf, insert the sphere here
        if (node->sphere != nullptr) {
            // There's already a sphere here, create a new non-leaf node and insert both spheres
            node->isLeaf = false;
            // Recursively insert the old sphere
            insert(node, node->sphere, x0, y0, z0, size);
            // Now insert the new sphere
            insert(node, sphere, x0, y0, z0, size);
            node->sphere = nullptr; // Remove the old sphere from this node
        } else {
            // Node is empty, simply store the sphere here
            node->sphere = sphere;
        }
    } else {
        // Node is non-leaf, update center of mass and recursively insert the sphere
        double xmid = x0 + size / 2.0;
        double ymid = y0 + size / 2.0;
        double zmid = z0 + size / 2.0;
        int octant = 0;
        if (sphere->x >= xmid) octant |= 1;
        if (sphere->y >= ymid) octant |= 2;
        if (sphere->z >= zmid) octant |= 4;
        
        if (node->children[octant] == nullptr) {
            // Create a new child node if it doesn't exist
            node->children[octant] = new QuadTreeNode;
            node->children[octant]->mass = 0.0;
            node->children[octant]->cx = xmid;
            node->children[octant]->cy = ymid;
            node->children[octant]->cz = zmid;
            node->children[octant]->isLeaf = true;
            node->children[octant]->sphere = nullptr;
        }

        // Update the center of mass and recursively insert the sphere into the child node
        node->mass += sphere->radius;
        insert(node->children[octant], sphere, x0, y0, z0, size / 2.0);
    }
}

// Function to recursively check for collisions with a sphere
void checkCollisions(QuadTreeNode* node, Sphere& sphere) {
    if (node->isLeaf) {
        if (node->sphere != nullptr && node->sphere != &sphere) {
            // Check for collision between the current sphere and the one in this leaf node
            if (spheresCollide(sphere, *(node->sphere))) {
                // Handle the collision here
                std::cout << "Collision detected!" << std::endl;
            }
        }
    } else {
        // Check if the distance between the center of mass and the sphere is greater than
        // some threshold or if the node is empty. If so, treat the node as a single body.
        double dx = node->cx - sphere.x;
        double dy = node->cy - sphere.y;
        double dz = node->cz - sphere.z;
        double distanceSquared = dx * dx + dy * dy + dz * dz;
        double size = 1.0; // Size of the current node (assuming it covers the whole space)
        double theta = 0.5; // Threshold for treating node as a single body (adjust as needed)

        if (size * size / distanceSquared < theta * theta || node->mass == 0.0) {
            // Treat the node as a single body and check for collision
            if (node->sphere != nullptr && node->sphere != &sphere) {
                // Check for collision between the current sphere and the one in this node
                if (spheresCollide(sphere, *(node->sphere))) {
                    // Handle the collision here
                    std::cout << "Collision detected!" << std::endl;
                }
            }
        } else {
            // Recursively check for collisions with child nodes
            for (int i = 0; i < 8; i++) {
                if (node->children[i] != nullptr) {
                    checkCollisions(node->children[i], sphere);
                }
            }
        }
    }
}


// Function to perform collision detection using Barnes-Hut
void collisionDetection(QuadTreeNode* root, std::vector<Sphere>& spheres) {
    // Insert all spheres into the quadtree
    for (Sphere& sphere : spheres) {
        insert(root, &sphere, 0.0, 0.0, 0.0, 1.0); // Assuming root node covers the whole space (0,0,0) to (1,1,1)
    }

    // Perform collision detection
    for (Sphere& sphere : spheres) {
        checkCollisions(root, sphere);
    }
}


int main() {
    // Create a root node for the quadtree
    QuadTreeNode* root = new QuadTreeNode;
    root->mass = 0.0;
    root->cx = 0.0;
    root->cy = 0.0;
    root->cz = 0.0;
    root->isLeaf = true;
    root->sphere = nullptr;

    // Create some spheres for collision detection
    std::vector<Sphere> spheres = {
        {0.2, 0.2, 0.2, 0.1},
        {0.5, 0.5, 0.5, 0.15},
        {0.8, 0.8, 0.8, 0.1}
    };

    // Perform collision detection
    collisionDetection(root, spheres);

    // Clean up
    delete root;

    return 0;
}
