#define MAX_TRAIL 1000 // May depend on length and sizes of objects.
#define TRAIL_SIZE 10

#define SPACE_COLOR BLACK
#define RESOLUTION_WIDTH 1920 // Modify for different resolution
#define RESOLUTION_HEIGHT 1080

#define SCALE 10 // Change the scale of simulation

#define SOFTENING 10 // Modify if the particles have SMALL distances to avoid numerical error
#define DELAY 1 // Microsecond

#define UPDATE_COM 1000 // Time interval for updating centre of mass display in simulation

// Limit amount for objects to display mechanical energy of each object in simulation screen
#define NUM_OBJECT_DISPLAY_ENERGY_LIMIT 10 

#include "n_body_simulation.h"


/**
 * Function: getAcc
 * -----------------------
 * Obtain the current acceleration of an object and store it in a passed-on
 * array. Also update the object's current potential energy.
 *
 * Parameters:
 *   objects    - Object[] - The array of objects whose acceleration are to be obtained.
 *   numObjects - int      - The number of objects in the array of objects.
 *   accX       - double[] - The x component of the acceleration of objects.
 *   accY       - double[] - The y component of the acceleration of objects.
 *   accZ       - double[] - The z component of the acceleration of objects.
 *   G          - double   - Gravitational Constant
 *
 * Example usage:
 *   getAcc(objects, numObjects, accX, accY);
 *
 * Notes:
 *   - Assumes accX and accY arrays are correctly sized to numObjects.
 */
void getAcc(Object objects[], int numObjects, double accX[], double accY[], double accZ[], double G){
    for (int i = 0; i < numObjects; i++) {
        double forceX = 0, 
               forceY = 0, 
               forceZ = 0,
               potEnergy = 0;

        for (int j = 0; j < numObjects; j++) {

            // Refer to the gravity of other objects except the object itself.
            if (i != j){
                double distX = objects[j].posX - objects[i].posX,
                       distY = objects[j].posY - objects[i].posY,
                       distZ = objects[j].posZ - objects[i].posZ,
                       dist = sqrt(distX*distX + distY*distY + distZ*distZ);

                if (dist < SOFTENING){
                    dist = SOFTENING;
                }

                // Avoid division by 0 if the object has the same coordinates as this current one.
                if (dist != 0){
                    // Gravitational force magnitude: 1 / r^2
                    double magForce = G * objects[i].mass * objects[j].mass / (dist*dist*dist);

                    // Add force components
                    forceX += magForce * distX;
                    forceY += magForce * distY;
                    forceZ += magForce * distZ;

                    // Sum potential energy
                    potEnergy += -G * objects[i].mass * objects[j].mass / dist;
                }
            }
        } 

        // Calculate the acceleration components of the object.
        accX[i] = forceX / objects[i].mass;
        accY[i] = forceY / objects[i].mass;
        accZ[i] = forceZ / objects[i].mass;

        // Update the object's potential energy.
        objects[i].potEnergy = potEnergy;
    }
}


/**
 * Function: getCentreOfMomentum
 * -----------------------
 * Obtain the current centre of momentum of a system of objects.
 *
 * Parameters:
 *   objects    - Object[] - The array of objects whose centre of momentum to be determined.
 *   numObjects - int      - The number of objects in the array of objects.
 *
 * Example usage:
 *   getCentreOfMomentum(objects, numObjects);
 */
Vector3D getCentreOfMomentum(Object objects[], int numObjects){
    Vector3D total = {0.0, 0.0, 0.0};
    double totalMass = 0;

    for (int i = 0; i < numObjects; i++) {
        total.posX += objects[i].mass * objects[i].velX;
        total.posY += objects[i].mass * objects[i].velY;
        total.posZ += objects[i].mass * objects[i].velZ;
        totalMass += objects[i].mass;
    }

    total.posX /= totalMass;
    total.posY /= totalMass;
    total.posZ /= totalMass;

    return total;
}


/**
 * Function: kick
 * -----------------------
 * Updates the velocities of objects in an N-body simulation using
 * the given accelerations and a timestep. This implements the
 * velocity "kick" step in the Leapfrog integration method.
 * Also calculates objects' kinetic energy.
 *
 * Parameters:
 *   objects    - Object[] - The array of objects whose velocities are to be updated.
 *   numObjects - int      - The number of objects in the array of objects.
 *   timeStep   - double   - The simulation timestep in seconds.
 *   accX       - double[] - The x component of the acceleration of objects.
 *   accY       - double[] - The y component of the acceleration of objects.
 *   accZ       - double[] - The z component of the acceleration of objects.
 *
 * Notes:
 *   - Part of the Leapfrog integration scheme (velocity kick step).
 *   - Assumes accX and accY arrays are correctly sized to numObjects.
 */
void kick(Object objects[], int numObjects, double timeStep, double accX[], double accY[], double accZ[]){
    Vector3D momentum = getCentreOfMomentum(objects, numObjects);
    for (int i = 0; i < numObjects; i++) {
        objects[i].velX += accX[i] * timeStep;
        objects[i].velY += accY[i] * timeStep;
        objects[i].velZ += accZ[i] * timeStep;
        objects[i].kinEnergy = 0.5 * objects[i].mass *
                    ((objects[i].velX - momentum.posX)*(objects[i].velX - momentum.posX) +
                     (objects[i].velY - momentum.posY)*(objects[i].velY - momentum.posY) +
                     (objects[i].velZ - momentum.posZ)*(objects[i].velZ - momentum.posZ));
    }
}


/**
 * Function: drift
 * -----------------------
 * Updates the position of objects in an N-body simulation using
 * the their kicked velocity and a timestep. This implements the
 * position "drift" step in the Leapfrog integration method.
 *
 * Parameters:
 *   objects    - Object[] - The array of objects whose velocities are to be updated.
 *   numObjects - int      - The number of objects in the array of objects.
 *   timeStep   - double   - The simulation timestep in seconds.
 *
 * Notes:
 *   - Part of the Leapfrog integration scheme (position drift step).
 *   - Assumes accX and accY arrays are correctly sized to numObjects.
 */
void drift(Object objects[], int numObjects, double timeStep){
    for (int i = 0; i < numObjects; i++){
        objects[i].posX += objects[i].velX * timeStep;
        objects[i].posY += objects[i].velY * timeStep;
        objects[i].posZ += objects[i].velZ * timeStep;
    }
}


/**
 * Function: getTotalEnergy
 * -----------------------
 * Get the total energy of objects in a system.
 *
 * Parameters:
 *   objects    - Object[]              - The array of objects whose total energy will be obtained.
 *   numObjects - int                   - The number of objects in the array of objects.
 * 
 * Return:
 *   double     - totalKE + 0.5*totalPE - The total energy of objects in the system.
 * 
 * Example usage:
 *   getTotalEnergy(objects, numObjects);
 */
double getTotalEnergy(Object objects[], int numObjects){
    double totalKE = 0, 
           totalPE = 0;

    for (int i = 0; i < numObjects; i++){
        totalKE += objects[i].kinEnergy;
        totalPE += objects[i].potEnergy;
    }

    return totalKE + 0.5*totalPE;
}


/**
 * Function: get2KEtoPERatio
 * -----------------------
 * Obtain the negative of the ratio between the two times the 
 * kinetic energy to the potential energy of a system of objects.
 *
 * Parameters:
 *   objects    - Object[] - The array of objects.
 *   numObjects - int      - The number of objects in the array of objects.
 *
 * Example usage:
 *   get2KEtoPERatio(objects, numObjects);
 * 
 * Notes:
 *   - Mainly used to determine whether a system lies in virial equilibrium.
 */
double get2KEtoPERatio(Object objects[], int numObjects){
    double totalKE = 0, 
           totalPE = 0;

    for (int i = 0; i < numObjects; i++){
        totalKE += objects[i].kinEnergy;
        totalPE += objects[i].potEnergy;
    }

    return -2*totalKE/(0.5*totalPE);
}


/**
 * Function: getTotalMomentum
 * -----------------------
 * Get the total momentum of objects in a system.
 *
 * Parameters:
 *   objects    - Object[] - The array of objects whose total energy will be obtained.
 *   numObjects - int      - The number of objects in the array of objects.
 * 
 * Return:
 *   Array of the dimension of momentum in 2D.
 * 
 * Example usage:
 *   getCentreOfMass(objects, numObjects);
 */
Vector3D getTotalMomentum(Object objects[], int numObjects){
    Vector3D total = {0.0, 0.0, 0.0};

    for (int i = 0; i < numObjects; i++) {
        total.posX += objects[i].mass * objects[i].velX;
        total.posY += objects[i].mass * objects[i].velY;
        total.posZ += objects[i].mass * objects[i].velZ;
    }

    return total;
}


/**
 * Function: getCenterOfMass
 * -----------------------
 * Obtain the current centre of mass of a system of objects.
 *
 * Parameters:
 *   objects    - Object[] - The array of objects whose center of mass to be determined.
 *   numObjects - int      - The number of objects in the array of objects.
 *
 * Example usage:
 *   getCentreOfMass(objects, numObjects);
 */
Vector3D getCenterOfMass(Object objects[], int numObjects){
    Vector3D com = {0};
    double totalMass = 0.0;

    for (int i = 0; i < numObjects; i++) {
        com.posX += objects[i].posX * objects[i].mass;
        com.posY += objects[i].posY * objects[i].mass;
        com.posZ += objects[i].posZ * objects[i].mass;
        totalMass += objects[i].mass;
    }

    com.posX /= totalMass;
    com.posY /= totalMass;
    com.posZ /= totalMass;

    return com;
}


/**
 * Function: updateObjectsEuler
 * -----------------------
 * Update the objects after a simulation time step using Euler's method.
 *
 * Parameters:
 *   objects     - Object[] - The array of objects whose total energy will be obtained.
 *   numObjects  - int      - The number of objects in the array of objects.
 *   totalEnergy - double   - The total of energy in the array of objects.
 *   G           - double   - Gravitational Constant
 */
void updateObjectsEuler(Object objects[], int numObjects, double timeStep, double G){
    double accX[numObjects];
    double accY[numObjects];
    double accZ[numObjects];

    // Compute acceleration at current sample time
    getAcc(objects, numObjects, accX, accY, accZ, G);

    // Update velocity
    for (int i = 0; i < numObjects; i++) {
        objects[i].velX += accX[i] * timeStep;
        objects[i].velY += accY[i] * timeStep;
        objects[i].velZ += accZ[i] * timeStep;
    }

    // Update position
    for (int i = 0; i < numObjects; i++) {
        objects[i].posX += objects[i].velX * timeStep;
        objects[i].posY += objects[i].velY * timeStep;
        objects[i].posZ += objects[i].velZ * timeStep;
    }
}


/**
 * Function: updateObjects
 * -----------------------
 * Update the objects after a simulation time step using Leapfrog method.
 *
 * Parameters:
 *   objects     - Object[] - The array of objects whose total energy will be obtained.
 *   numObjects  - int      - The number of objects in the array of objects.
 *   totalEnergy - double   - The total of energy in the array of objects.
 *   G           - double   - Gravitational Constant
 */
void updateObjects(Object objects[], int numObjects, double timeStep, double G){
    double accX[numObjects],
           accY[numObjects],
           accZ[numObjects];
    
    // 1. a(t)
    getAcc(objects, numObjects, accX, accY, accZ, G);

    // 2. Kick half step
    kick(objects, numObjects, 0.5*timeStep, accX, accY, accZ);

    // 3. Drift full step
    drift(objects, numObjects, timeStep);

    // 4. a(t+dt)
    getAcc(objects, numObjects, accX, accY, accZ, G);
   
    // 5. Kick second half step
    kick(objects, numObjects, 0.5*timeStep, accX, accY, accZ);
}


/**
 * Function: writeCSVRow
 * -----------------------
 * Write the status of each object into a CSV file.
 *
 * Parameters:
 *   objects     - Object[] - The array of objects.
 *   numObjects  - int      - The number of objects in the array of objects.
 *   totalEnergy - double   - The total of energy in the array of objects.
 *   time        - double   - Current time of the simulation.
 *   filePointer - *FILE    - The file pointer to writ the CSV file to.
 *   momentum    - Vector3D - Array of dimensions for momentum in 3D.
 */
void writeCSVRow(Object objects[], int numObjects, double totalEnergy, double time, FILE *filePointer, Vector3D momentum){
    // For each object, write one row
    for (int i = 0; i < numObjects; i++) {
        fprintf(filePointer, WRITE_FORMAT,
                time,            // simulation time
                i,               // object index
                objects[i].posX,
                objects[i].posY,
                objects[i].posZ,
                objects[i].velX,
                objects[i].velY,
                objects[i].velZ,
                objects[i].mass,
                objects[i].radius,
                objects[i].kinEnergy + objects[i].potEnergy,
                totalEnergy,
                momentum.posX,
                momentum.posY,
                momentum.posZ);    // same value for every row this timestep
    }
}


/**
 * Function: displayObjects2D
 * -----------------------
 * Display the simulation in a 2D screen using Raylib library.
 *
 * Parameters:
 *   objects     - Object[] - The array of objects.
 *   numObjects  - int      - The number of objects in the array of objects.
 *   energy      - double   - The total of energy in the array of objects.
 *   momentum    - Vector3D - The momentum of array of objects.
 *   step        - int      - The current time step.
 *   title       - char     - The title of the simulation.
 */
void displayObjects2D(Object objects[], int numObjects, double energy, Vector3D momentum, int step, const char *title){
    static Vector3D currCom;
    Vector3D com = getCenterOfMass(objects, numObjects);

    BeginDrawing();
    ClearBackground(SPACE_COLOR);

    char buffer[256];

    int size = 30, start = 0, increment = 40;

    // Title
    snprintf(buffer, sizeof(buffer), "Title: %s", title);
    DrawText(buffer, 10, start += increment, size, YELLOW);

    // Step
    snprintf(buffer, sizeof(buffer), "Step Number: %d", step);
    DrawText(buffer, 10, start += increment, size, WHITE);

    // Total Energy
    snprintf(buffer, sizeof(buffer), "Total Energy: %.6e", energy);
    DrawText(buffer, 10, start += increment, size, WHITE);

    // Total Energy
    snprintf(buffer, sizeof(buffer), "2KE/PE: %.6e", (get2KEtoPERatio(objects, numObjects)));
    DrawText(buffer, 10, start += increment, size, WHITE);

    if (step % UPDATE_COM == 0){
        // Center of Mass
        snprintf(buffer, sizeof(buffer), "Center of Mass (X, Y, Z): (%f, %f, %f)", com.posX, com.posY, com.posZ);
        DrawText(buffer, 10, start += increment, size, WHITE);
        currCom = com;
    } else {
        snprintf(buffer, sizeof(buffer), "Center of Mass (X, Y, Z): (%f, %f, %f)", currCom.posX, currCom.posY, currCom.posZ);
        DrawText(buffer, 10, start += increment, size, WHITE);
    }

    // Momentum X
    snprintf(buffer, sizeof(buffer), "Momentum X: %.6f", momentum.posX);
    DrawText(buffer, 10, start += increment, size, WHITE);

    // Momentum Y
    snprintf(buffer, sizeof(buffer), "Momentum Y: %.6f", momentum.posY);
    DrawText(buffer, 10, start += increment, size, WHITE);

    // Momentum Z
    snprintf(buffer, sizeof(buffer), "Momentum Z: %.6f", momentum.posZ);
    DrawText(buffer, 10, start += increment, size, WHITE);

    for (int i = 0; i < numObjects; i++) {

        // Position relative to center of mass
        double relX = objects[i].posX - com.posX;
        double relY = objects[i].posY - com.posY;

        int x = (int)(RESOLUTION_WIDTH  / 2 + relX * SCALE);
        int y = (int)(RESOLUTION_HEIGHT / 2 + relY * SCALE);

        int radius = (int)(objects[i].radius * SCALE);
        DrawCircle(x, y, radius, objects[i].color);
        
        // Choose diff color if object escape

        Color color = WHITE;
        if (objects[i].kinEnergy + objects[i].potEnergy > 0){
            color = GRAY;
        }

        if (numObjects <= NUM_OBJECT_DISPLAY_ENERGY_LIMIT){
            snprintf(buffer, sizeof(buffer), "M = %.0f object | ME = %.2f", 
                                        objects[i].mass, objects[i].kinEnergy + objects[i].potEnergy);
            DrawText(buffer, 10, start += increment, size, color);
        }

        if (MAX_TRAIL > 0){
            // Add current position to trail
            objects[i].trail[objects[i].trailIndex] = (Vector3D){objects[i].posX, objects[i].posY, objects[i].posZ};
            objects[i].trailIndex = (objects[i].trailIndex + 1) % MAX_TRAIL;
            if (objects[i].trailCount < MAX_TRAIL) objects[i].trailCount++;
            
            // Draw trails
            for (int j = 0; j < objects[i].trailCount - 1; j++) {
                int idx1 = (objects[i].trailIndex + MAX_TRAIL - objects[i].trailCount + j) % MAX_TRAIL;
                int idx2 = (objects[i].trailIndex + MAX_TRAIL - objects[i].trailCount + j + 1) % MAX_TRAIL;

                int x1 = (int)(RESOLUTION_WIDTH  / 2 + (objects[i].trail[idx1].posX - com.posX) * SCALE);
                int y1 = (int)(RESOLUTION_HEIGHT / 2 + (objects[i].trail[idx1].posY - com.posY) * SCALE);

                int x2 = (int)(RESOLUTION_WIDTH  / 2 + (objects[i].trail[idx2].posX - com.posX) * SCALE);
                int y2 = (int)(RESOLUTION_HEIGHT / 2 + (objects[i].trail[idx2].posY - com.posY) * SCALE);

                float fade = 0.2f + 0.8f * ((float)j / (objects[i].trailCount - 1));  // scale fade with actual count
                DrawLine(x1, y1, x2, y2, Fade(objects[i].color, fade));
            }
        }
    }

    // Draw center of mass
    int comX = RESOLUTION_WIDTH / 2;
    int comY = RESOLUTION_HEIGHT / 2;

    // Crosshair marker for COM
    DrawLine(comX - SCALE/10, comY, comX + SCALE/10, comY, WHITE);
    DrawLine(comX, comY - SCALE/10, comX, comY + SCALE/10, WHITE);

    EndDrawing();
}


/**
 * Function: simulate
 * -----------------------
 * Simulate the gravitating objects given a simulation time step and number of steps.
 *
 * Parameters:
 *   objects    - Object[] - The array of objects whose total energy will be obtained.
 *   numObjects - int      - The number of objects in the array of objects.
 *   startTime  - double   - The starting time in seconds.
 *   timeStep   - double   - The time step in seconds to be simulated
 *   numSteps   - int      - The number of steps to be simulated
 *   filename   - *char    - The file name in .csv to output simulation results.
 *   G          - double   - Gravitational Constant
 *   title      - *char    - Title of simulation
 */
void simulate(Object objects[], int numObjects, double startTime, double timeStep, int numSteps, const char *filename, double G, const char *title){
    FILE *filePointer = fopen(filename, OUTPUT_WRITE);
    fprintf(filePointer, OUTPUT_HEADER);
    InitWindow(RESOLUTION_WIDTH, RESOLUTION_HEIGHT, "Simulate");

    if (MAX_TRAIL > 0){
        for (int i = 0; i < numObjects; i++) {
            objects[i].trailIndex = 0;
            objects[i].trailCount = 0;
            for (int j = 0; j < MAX_TRAIL; j++) {
                objects[i].trail[j] = (Vector3D){objects[i].posX, objects[i].posY, objects[i].posZ};  // optional: clear memory
            }
        }
    }

    double currTime = startTime;
    for (int step = 0; step < numSteps; step++) {
        updateObjects(objects, numObjects, timeStep, G);
        // updateObjectsEuler(objects, numObjects, timeStep, G);
        
        // Check for the total energy of the system.
        double totEnergy = getTotalEnergy(objects, numObjects);

        // Check for the total momentum of the system.
        Vector3D momentum = getTotalMomentum(objects, numObjects);

        // Output results to the csv file.
        // writeCSVRow(objects, numObjects, totEnergy, currTime, filePointer, momentum);

        // Simulate
        displayObjects2D(objects, numObjects, totEnergy, momentum, step, title);

        // Update step
        currTime += timeStep;

        // Delay if desired
        usleep(DELAY);
    }

    fclose(filePointer);
    CloseWindow();
}


// Sample test 1
void test1000EqualBodiesSquare() {
    const double timeStep = 0.1;     
    const int numSteps = 5000;      
    const double G = 1.0;
    double mass = 1.0;
    double boxSize = 100.0;  

    const int numBodies = 1000;
    Object objects[numBodies];

    for (int i = 0; i < numBodies; i++) {
        objects[i].mass = mass;
        objects[i].radius = 1.0;

        // Random position in cube [-boxSize, boxSize] in each axis
        objects[i].posX = ((double)rand() / RAND_MAX) * 2 * boxSize - boxSize;
        objects[i].posY = ((double)rand() / RAND_MAX) * 2 * boxSize - boxSize;
        objects[i].posZ = ((double)rand() / RAND_MAX) * 2 * boxSize - boxSize;

        // Small random velocity for virialization
        double speed = 0.5;
        double alpha = ((double)rand() / RAND_MAX) * 2 * M_PI;
        double beta  = ((double)rand() / RAND_MAX) * M_PI;
        objects[i].velX = speed * sin(beta) * cos(alpha);
        objects[i].velY = speed * sin(beta) * sin(alpha);
        objects[i].velZ = speed * cos(beta);

        objects[i].color = WHITE;
    }

    simulate(objects, numBodies, 0.0, timeStep, numSteps, "thousand_equal_mass_square.csv", G, "1000EqualBodies");
}


// Sample test 2
void test3Bodies() {
    const double timeStep = 0.01;
    const int numSteps = 130000; // simulate long enough to see possible ejection
    const double G = 1.0;       // gravitational constant

    int numBodies = 3;
    Object objects[numBodies];

    // Masses: heavy, medium, light
    objects[0].mass = 100; // central heavy mass
    objects[1].mass = 10;  // medium mass
    objects[2].mass = 1;   // light mass

    // Central heavy mass at origin
    objects[0].posX = 0.0;
    objects[0].posY = 0.0;
    objects[0].posZ = 0.0;
    objects[0].velX = 0;
    objects[0].velY = 0.0;
    objects[0].velZ = 0;
    objects[0].radius = 1;
    objects[0].color = RED;

    // Medium mass
    double R2 = 20;
    double v2 = sqrt(G * objects[0].mass / R2); // circular orbit velocity
    objects[1].posX = R2;
    objects[1].posY = 0.0;
    objects[1].posZ = 0.0;
    objects[1].velX = 0;
    objects[1].velY = v2;
    objects[1].velZ = 0;
    objects[1].radius = 0.5;
    objects[1].color = YELLOW;

    // Light mass
    double R3 = 15.7;
    double v3 = sqrt(G * objects[0].mass / R3); // circular orbit velocity
    objects[2].posX = -R3;
    objects[2].posY = 0.0;
    objects[2].posZ = 0.0;
    objects[2].velX = 0;
    objects[2].velY = -v3; // opposite direction
    objects[2].velZ = 0;
    objects[2].radius = 0.25;
    objects[2].color = GREEN;

    objects[2].velY *= 0.25; // Adjustment

    // Run simulation
    simulate(objects, numBodies, 0.0000000, timeStep, numSteps, "three_body_hierarchy.csv", G, "ThreeBodiesTest");
}


/*
Do not run main, just type the following in terminal:

    gcc n_body_simulation.c -o n_body_simulation.exe -lraylib -lopengl32 -lgdi32 -lwinmm -lws2_32 -lm
    ./n_body_simulation.exe
*/
int main() {
    // test1000EqualBodiesSquare();
    test3Bodies();
    return 0;
}