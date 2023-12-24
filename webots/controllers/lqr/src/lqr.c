#include <webots/robot.h>
#include <webots/motor.h>
#include <webots/gps.h>
#include <webots/accelerometer.h>
#include <webots/position_sensor.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// global variables declaration
#define TIME_STEP 10
#define PI 3.14159
#define ARRAY_SIZE 2500 // sampling times
#define n_state 4
static double t0 = 0.0;   // start time
static double t1 = 25.0;  // end time
static double Ts = 0.010; // sampling period
static double g = 9.81;   // gravitational acceleration

/* reference: [1] jitendra singh (2023). Inverted Pendulum: Swing up & LQR Control design MATLAB Code (https://www.mathworks.com/matlabcentral/fileexchange/97182-inverted-pendulum-swing-up-lqr-control-design-matlab-code),
MATLAB Central File Exchange.*/

// saturated function
double saturated(double x, double x_max, double x_min){
    if (x > x_max){
        return x_max;
    }
    else if (x < x_min){
        return x_min;
    }
    else{
        return x;
    }
}

// symbolic function
double sign(double x){
    if (x > 0){
        return 1.0;
    }
    else if (x < 0){
        return -1.0;
    }
    else{
        return 0.0;
    }
}

struct _archive{
    double position_cart_archive[ARRAY_SIZE];
    double velocity_cart_archive[ARRAY_SIZE];
    double pendulum_angle_archive[ARRAY_SIZE];
    double pendulum_velocity_archive[ARRAY_SIZE];
    double control_archive[ARRAY_SIZE];
} archive;

struct _system_state{
    double position_cart;     // actual cart position
    double velocity_cart;     // actual cart velocity
    double pendulum_angle;    // actual pendulum angle of inverted pendulum
    double pendulum_velocity; // actual pendulum velocity of inverted pendulum
} system_state;

struct _controller{
    int controller_flag;            // controller switch flag
    double last_position_cart;      // actual cart position at last sampling moment
    double last_pendulum_angle;     // actual pendulum angle at last sampling moment
    double error_position_cart;     // error of cart position
    double pendulum_angle_measured; // measured pendulum angle of inverted pendulum
    double acceleration_cart;       // cart acceleration
    double K_position_cart;         // control gain of cart position in swing up controller
    double K_velocity_cart;         // control gain of cart velocity in swing up controller
    double K_lqr[n_state];          // control gain of LQR
    double control_swing_up;        // control input force of swing up controller, reference linear motor force
    double control_volt;            // control input voltage for linear motor
    double control_lqr;             // control input force of LQR controller
    double control;                 // control input force
    double state[n_state];          // system state variable
    double state_desired[n_state];  // desired system state variable
    double theta;                   // computed pendulum angle of inverted pendulum for LQR
    double r;                       // motor pinion radius
    double b;                       // viscous damping at pivot of Pendulum
    double c;                       // viscous friction coefficient of cart
    double L;                       // motor inductance
    double Rm;                      // motor armature resistance
    double kb;                      // motor back emf constant
    double kt;                      // motor torque constant
    double m_cart;                  // mass of cart
    double m_pendulum;              // mass of pendulum
    double moment_of_inertia;       // moment of inertia of pendulum
    double length_pendulum;         // pendulum length from pivot to centre of gravity
    double energy_desired;          // desired energy of pendulum
    double energy;                  // total energy of pendulum
    double n_saturated;             // factorization of the saturation function
    double k_swing;
    double position_cart_baseline;
} controller;

void CONTROLLER_init(){
    wb_robot_init();
    WbDeviceTag linearmotor = wb_robot_get_device("linearmotor"); // initialize linear motor
    controller.last_position_cart = 0.0;
    controller.last_pendulum_angle = PI;
    controller.r = 0.006;                     // motor pinion radius
    controller.m_cart = 0.135;                // mass of cart
    controller.m_pendulum = 0.1;              // mass of pendulum
    controller.moment_of_inertia = 0.0007176; // moment of inertia of pendulum
    controller.length_pendulum = 0.2;         // pendulum length from pivot to centre of gravity
    controller.b = 0.00007892;                // viscous damping at pivot of Pendulum
    controller.L = 0.046;                     // motor inductance
    controller.Rm = 12.5;                     // motor armature resistance
    controller.kb = 0.031;                    // motor back emf constant
    controller.kt = 0.031;                    // motor torque constant
    controller.c = 0.63;                      // viscous friction coefficient of cart
    controller.n_saturated = 3.0;             // factorization of the saturation function
    controller.k_swing = 1.2;

    // desired system state variable
    controller.state_desired[1] = PI;
    for (int j = 0; j < n_state; j++){
        printf("state_desired[%d]: %f\n", j, controller.state_desired[j]);
    }

    // numerical solution of matrix algebraic riccati equation from matlab: K = lqr(A,B,Q,R);
    double K_value[n_state] = {-75.5929, 212.5131, -63.5585, 19.3977};
    for (int j = 0; j < n_state; j++){
        controller.K_lqr[j] = 4.0 * K_value[j];
    }
}

double CONTROLLER_realize(int i){
    WbDeviceTag linearmotor = wb_robot_get_device("linearmotor");

    // obtain actual pendulum angle of inverted pendulum
    WbDeviceTag positionsensor2 = wb_robot_get_device("positionsensor2");
    wb_position_sensor_enable(positionsensor2, TIME_STEP);
    controller.pendulum_angle_measured = wb_position_sensor_get_value(positionsensor2);
    system_state.pendulum_angle = controller.pendulum_angle_measured + PI;
    printf("pendulum_angle: %lf\n", system_state.pendulum_angle);
    archive.pendulum_angle_archive[i] = system_state.pendulum_angle;

    system_state.pendulum_velocity = (system_state.pendulum_angle - controller.last_pendulum_angle) / Ts; // actual pendulum velocity of inverted pendulum
    printf("pendulum_velocity: %lf\n", system_state.pendulum_velocity);
    archive.pendulum_velocity_archive[i] = system_state.pendulum_velocity;

    WbDeviceTag positionsensor1 = wb_robot_get_device("positionsensor1");
    wb_position_sensor_enable(positionsensor1, TIME_STEP);
    system_state.position_cart = wb_position_sensor_get_value(positionsensor1); // actual cart position
    printf("position_cart: %lf\n", system_state.position_cart);
    archive.position_cart_archive[i] = system_state.position_cart;

    system_state.velocity_cart = (system_state.position_cart - controller.last_position_cart) / Ts; // actual cart velocity
    printf("velocity_cart: %lf\n", system_state.velocity_cart);
    archive.velocity_cart_archive[i] = system_state.velocity_cart;

    controller.state[0] = system_state.position_cart;     // first state variable is cart position
    controller.state[1] = system_state.pendulum_angle;    // second state variable is pendulum angle of inverted pendulum
    controller.state[2] = system_state.velocity_cart;     // third state variable is cart velocity
    controller.state[3] = system_state.pendulum_velocity; // fourth state variable is pendulum velocity of inverted pendulum

    for (int j = 0; j < n_state; j++){
        printf("controller.state[%d]: %f\n", j, controller.state[j]);
    }

    // LQR state feedback control law
    controller.control_volt = 0.0;
    for (int j = 0; j < n_state; j++){
        controller.control_volt += - controller.K_lqr[j] * (controller.state[j] - controller.state_desired[j]);
    }
    printf("control_volt: %lf\n", controller.control_volt);

    // control input voltage with saturation constraint
    controller.control_volt = saturated(controller.control_volt, 12.0, -12.0);
    // control input torque of LQR
    controller.control_lqr =  (controller.kt * controller.control_volt * controller.r - controller.kt * controller.kb * system_state.velocity_cart) / (controller.Rm * pow(controller.r, 2));
    // printf("control_lqr: %f\n", controller.control_lqr);

    controller.control = controller.control_lqr; // lqr controller
    printf("control_lqr: %f\n", controller.control_lqr);
    archive.control_archive[i] = controller.control;
    wb_motor_set_force(linearmotor, - controller.control);
    controller.last_position_cart = system_state.position_cart;
    controller.last_pendulum_angle = system_state.pendulum_angle;
}

void saveArchiveToTxt(double *archive, int size, const char *filename){

    FILE *file = fopen(filename, "w+");

    if (file == NULL){
        perror("Failed to open file");
        exit(1);
    }
    else{
        for (int i = 0; i < size; i++){
            fprintf(file, "%lf\n", archive[i]);
        }
        fclose(file);
        printf("Saved to file %s\n", filename);
    }
}

void saveArchive(){
    saveArchiveToTxt(archive.position_cart_archive, ARRAY_SIZE, "../../../report/position_cart.txt");
    saveArchiveToTxt(archive.velocity_cart_archive, ARRAY_SIZE, "../../../report/velocity_cart.txt");
    saveArchiveToTxt(archive.pendulum_angle_archive, ARRAY_SIZE, "../../../report/pendulum_angle.txt");
    saveArchiveToTxt(archive.pendulum_velocity_archive, ARRAY_SIZE, "../../../report/pendulum_velocity.txt");
    saveArchiveToTxt(archive.control_archive, ARRAY_SIZE, "../../../report/control.txt");
}

int main(int argc, char **argv){

    CONTROLLER_init(); // initialize controller parameter
    // PLANT_init();      // initialize plant parameter

    int i = 0;
    while (wb_robot_step(TIME_STEP) != -1){
        for (int j = 0; j < 30; j++){
            printf("*");
        }
        printf("\n");
        double time = i * TIME_STEP / 1000.0 + t0;
        printf("time at step %d: %f\n", i, time);

        if (time > t1){
            break;
        }

        CONTROLLER_realize(i);
        // PLANT_realize(i);
        i++;
    }

    saveArchive();

    wb_robot_cleanup();

    return 0;
}
