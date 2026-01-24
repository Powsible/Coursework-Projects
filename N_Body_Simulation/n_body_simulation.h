#ifndef N_BODY_CODE_H 
#define N_BODY_CODE_H

#ifdef __cplusplus
    extern "C" {
#endif

// ------------------------------
// Includes
// ------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <raylib.h>
#include <unistd.h>

// ------------------------------
// Macro definitions
// ------------------------------
#define M_PI 3.14159265358979323846 // Some reason the math.h doesn't work
#define OUTPUT_WRITE "w"
#define OUTPUT_HEADER "Time(s),Object ID,posX,posY,posZ,velX,velY,velZ,Mass,Radius,Object Energy,System Total Energy,Momentum X,Momentum Y,Momentum Z\n"
#define WRITE_FORMAT "%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n"

#ifndef MAX_TRAIL
#define MAX_TRAIL 1000 // default value
#endif

// ------------------------------
// Type definitions
// ------------------------------

/**
 * Struct: Vector3D
 * ------------------
 * Object vector in 3D.
 *
 * Members:
 *   posX      - double - The x component of the position of the object.
 *   posY      - double - The y component of the position of the object.
 *   posZ      - double - The z component of the position of the object.
 */
typedef struct{
    double posX,
           posY,
           posZ;
} Vector3D;

/**
 * Struct: Object
 * ------------------
 * Representation of a gravitating object
 *
 * Members:
 *   posX       - double     - The x component of the position of the object.
 *   posY       - double     - The y component of the position of the object.
 *   posZ       - double     - The z component of the position of the object.
 *   velX       - double     - The x component of the velocity of the object.
 *   velY       - double     - The y component of the velocity of the object.
 *   velZ       - double     - The z component of the velocity of the object.
 *   mass       - double     - The mass of the object.
 *   potEnergy  - double     - The potential energy of the object.
 *   kinEnergy  - double     - The kinetic energy of the object.
 *   radius     - double     - The radius of the object.
 *   color      - Color      - Object color
 *   trail      - Vector3D[] - Trail of the object.
 *   trailIndex
 */
typedef struct {
    double posX, posY, posZ,
           velX, velY, velZ,
           mass, 
           potEnergy, 
           kinEnergy,
           radius;
    Color color;

    Vector3D trail[MAX_TRAIL]; // store previous positions
    int trailIndex, 
        trailCount; // current index
} Object;

// ------------------------------
// Function declarations
// ------------------------------

void getAcc(Object objects[], int numObjects, double accX[], double accY[], double accZ[], double G);
Vector3D getCentreOfMomentum(Object objects[], int numObjects);
void kick(Object objects[], int numObjects, double timeStep, double accX[], double accY[], double accZ[]);
void drift(Object objects[], int numObjects, double timeStep);
double getTotalEnergy(Object objects[], int numObjects);
double get2KEtoPERatio(Object objects[], int numObjects);
Vector3D getTotalMomentum(Object objects[], int numObjects);
Vector3D getCenterOfMass(Object objects[], int numObjects);
void updateObjectsEuler(Object objects[], int numObjects, double timeStep, double G);
void updateObjects(Object objects[], int numObjects, double timeStep, double G);
void writeCSVRow(Object objects[], int numObjects, double totalEnergy, double time, FILE *filePointer, Vector3D momentum);
void displayObjects2D(Object objects[], int numObjects, double energy, Vector3D momentum, int step, const char *title);
void simulate(Object objects[], int numObjects, double startTime, double timeStep, int numSteps, const char *filename, double G, const char *title);

#ifdef __cplusplus
}
#endif

#endif