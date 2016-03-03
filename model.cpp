#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

enum StreamMode
{
    STREAM_NONE,
    STREAM_POSITIONS,
    STREAM_ANGLES
};

StreamMode streamMode = STREAM_ANGLES;
int streamCycle = 1000;
bool printResult = true;

// State Struct {{{

struct State
{
    double angle1;
    double angle2;
    double dangle1;
    double dangle2;
};

State operator+(State a, State b)
{
    State c;
    c.angle1 = a.angle1 + b.angle1;
    c.angle2 = a.angle2 + b.angle2;
    c.dangle1 = a.dangle1 + b.dangle1;
    c.dangle2 = a.dangle2 + b.dangle2;
    return c;
}

State operator*(State a, double b)
{
    State c;
    c.angle1 = a.angle1 * b;
    c.angle2 = a.angle2 * b;
    c.dangle1 = a.dangle1 * b;
    c.dangle2 = a.dangle2 * b;
    return c;
}

State operator*(double a, State b)
{
    return b * a;
}

// }}}

// Typedefs {{{

typedef double (*Fn)(double t);

// }}}

// Param Struct {{{

struct Param
{
    double mass;
    double g;
    double length;
    double damp1;
    double damp2;
    Fn accelx;
    Fn accely;
    Fn accelz;
};

// }}}

// Numerical Differentiation {{{

double diff(Fn fn, double x)
{
    double dx = 1.0e-6;
    double x1 = x - dx;
    double x2 = x + dx;
    double y1 = fn(x1);
    double y2 = fn(x2);
    return (y2 - y1) / (x2 - x1);
}

// }}}

// The Model {{{

State getDerivative(State &state, Param &param, double time)
{
    State d;

    double t = time;
    double th1 = state.angle1;
    double th2 = state.angle2;
    double dth1 = state.dangle1;
    double dth2 = state.dangle2;

    double m = param.mass;
    double g = param.g;
    double l = param.length;
    double ax = param.accelx(t);
    double ay = param.accely(t);
    double az = param.accelz(t);

    double c1 = cos(th1);
    double c2 = cos(th2);
    double s1 = sin(th1);
    double s2 = sin(th2);

    d.dangle1 = dth2*dth2*s1*c1 - (ax*c1*c2 + ay*c1*s2 + (az + g)*s1) / l;
    d.dangle2 = (ax*s2 - ay*c2) / (l * s1) - 2*dth1*dth2 * c1/s1;
    d.dangle1 -= param.damp1 * dth1;
    d.dangle2 -= param.damp2 * dth2;
    d.angle1 = dth1;
    d.angle2 = dth2;

    return d;
}

// }}}

// The Solver {{{

State rk4(State &state, Param &param, double step, double time)
{
    State k1 = getDerivative(state, param, time);

    State s2 = state + 0.5 * step * k1;
    State k2 = getDerivative(s2, param, time + 0.5 * step);

    State s3 = state + 0.5 * step * k2;
    State k3 = getDerivative(s3, param, time + 0.5 * step);

    State s4 = state + step * k3;
    State k4 = getDerivative(s4, param, time + step);

    return state + step / 6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4);
}

// }}}

// Test cases for convenience
int testNumber = 8;
double testAngle[] =
{
    1.5,
    1.2,
    1.0,
    0.5,
    -0.5,
    -0.6,
    -0.8,
    -1.0,
    -1.2,
    -1.4,
    -1.5,
    -2.7,
    -2.9,
    -3.0
};

double testRadPerSec[] =
{
    9.77591445809,
    4.27289896361,
    3.43948436183,
    2.40128546859,
    13.4954970057,
    7.65126437832,
    6.14655625872,
    6.24720671618,
    7.23505204194,
    10.3064690187,
    15.8955607605,
    14.3183359035,
    3.38668751069,
    2.12773706908
};
double getAngVel(double angle1)
{
    angle1 *= -1;
    return sqrt(9.81 * tan(angle1) / (sin(angle1) - 0.45));
}

double bobRadius = 0.15;
double pivotRadius = 0.45;
double pivotRadPerSecInitial = testRadPerSec[testNumber];
double pivotRadPerSecChange = 0.0;

// Pivot Motion {{{

double pivotRadPerSec(double time)
{
    return pivotRadPerSecInitial - pivotRadPerSecChange * time;
}
double pivotPhase(double time)
{
    return pivotRadPerSecInitial * time - 0.5 * pivotRadPerSecChange * time * time;
}
double pivotX(double time)
{
    return pivotRadius * cos(pivotPhase(time));
}
double pivotY(double time)
{
    return pivotRadius * sin(pivotPhase(time));
}
double pivotZ(double time)
{
    return 0.0;
}
double pivotVelX(double time)
{
    return diff(pivotX, time);
}
double pivotVelY(double time)
{
    return diff(pivotY, time);
}
double pivotVelZ(double time)
{
    return diff(pivotZ, time);
}
double accelPivotX(double time)
{
    return diff(pivotVelX, time);
}
double accelPivotY(double time)
{
    return diff(pivotVelY, time);
}
double accelPivotZ(double time)
{
    return diff(pivotVelZ, time);
}

// }}}

// Print util {{{

void printStatePlot(double time, State state, Param param)
{
    double px = pivotX(time);
    double py = pivotY(time);
    double pz = pivotZ(time);
    double l = param.length;
    double x = px + l * sin(state.angle1) * cos(state.angle2);
    double y = py + l * sin(state.angle1) * sin(state.angle2);
    double z = pz + l * cos(state.angle1);
    printf("%f %f %f %f %f %f %f\n", time, px, py, pz, x, y, z);
}

void printAngles(double time, State state, Param param)
{
    double anglediff = pivotPhase(time) - state.angle2;
    double cad = cos(anglediff);
    double sad = sin(anglediff);
    double dangle1term1 = state.dangle2 * state.dangle2 * sin(state.angle1) * cos(state.angle1);
    double dangle2term2 = -2.0 * state.dangle1 * state.dangle2 / tan(state.angle1);
    double dangle1termg = -param.g / param.length * sin(state.angle1);
    printf("%f %f %f %f %f %f %f %f %f %f\n", time, anglediff, state.dangle2, state.angle1, state.dangle1, cad, sad, dangle1term1, dangle2term2, dangle1termg);
}

// }}}

// TEST: Angvel vs Steady angle1 {{{

void testAngVelVSSteadyAngle(State initial, Param param, double step)
{
    for (double angVel = 8; angVel > 4; angVel -= 0.1)
    {
        pivotRadPerSecInitial = angVel;
        initial.dangle2 = angVel;
        State state = initial;
        double angle1Prev = 0.0;
        int sameCounter = 0;
        int counter = 0;
        for (double time = 0.0; sameCounter < 10000; time += step)
        {
            state = rk4(state, param, step, time);

            if (streamMode != STREAM_NONE)
            {
                counter++;
                if (counter >= streamCycle)
                {
                    counter = 0;
                    switch (streamMode)
                    {
                        case STREAM_POSITIONS:
                            printStatePlot(time, state, param);
                            break;
                        case STREAM_ANGLES:
                            printAngles(time, state, param);
                            break;
                    }
                }
            }

            double difference = state.angle1 - angle1Prev;
            if (difference < 0) difference *= -1;
            if (difference < 0.000001)
            {
                sameCounter++;
            }
            else
            {
                sameCounter = 0;
                angle1Prev = state.angle1;
            }
        }
        if (printResult) printf("%f %f\n", angVel, angle1Prev);
    }
}

// }}}

// TEST Single Case {{{

void testSingleCase(State initial, Param param, double step)
{
    State state = initial;
    int counter = 0;
    for (double time = 0.0; true; time += step)
    {
        state = rk4(state, param, step, time);

        if (streamMode != STREAM_NONE)
        {
            counter++;
            if (counter >= streamCycle)
            {
                counter = 0;
                switch (streamMode)
                {
                    case STREAM_POSITIONS:
                        printStatePlot(time, state, param);
                        break;
                    case STREAM_ANGLES:
                        printAngles(time, state, param);
                        break;
                }
            }
        }
    }
}

// }}}

int main(int argc, char ** argv)
{
    Param param;
    param.mass = 1.0;
    param.g = 9.81;
    param.length = 1.0;
    param.damp1 = 0.0;
    param.damp2 = 0.0;
    param.accelx = accelPivotX;
    param.accely = accelPivotY;
    param.accelz = accelPivotZ;

    State initial;
    initial.angle1 = testAngle[testNumber];//-0.65;//expectedAngle1;//-45 / 180.0 * 3.141592653;
    initial.angle2 = 0.0;//0.005;
    initial.dangle1 = 0.0;
    initial.dangle2 = pivotRadPerSec(0); //5.0;

    double step = 0.00002;
    void (*testRoutine)(State, Param, double step) = testSingleCase;
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "sp") == 0) streamMode = STREAM_POSITIONS;
        else if (strcmp(argv[i], "sa") == 0) streamMode = STREAM_ANGLES;
        else if (strcmp(argv[i], "sn") == 0) streamMode = STREAM_NONE;
        else if (strcmp(argv[i], "pr") == 0) printResult = true;
        else if (strcmp(argv[i], "npr") == 0) printResult = false;
        else if (strcmp(argv[i], "tsc") == 0) testRoutine = testSingleCase;
        else if (strcmp(argv[i], "tangvel") == 0) testRoutine = testAngVelVSSteadyAngle;
        else if (strcmp(argv[i], "-damp2") == 0 && i + 1 < argc)
        {
            i++;
            param.damp2 = atof(argv[i]);
        }
        else if (strcmp(argv[i], "-damp1") == 0 && i + 1 < argc)
        {
            i++;
            param.damp1 = atof(argv[i]);
        }
        else if (strcmp(argv[i], "-testnum") == 0 && i + 1 < argc)
        {
            i++;
            testNumber = atoi(argv[i]);
            pivotRadPerSecInitial = testRadPerSec[testNumber];
            initial.angle1 = testAngle[testNumber];
        }
        else if (strcmp(argv[i], "-step") == 0 && i + 1 < argc)
        {
            i++;
            step = atof(argv[i]);
        }
        else if (strcmp(argv[i], "-angle1") == 0 && i + 1 < argc)
        {
            i++;
            initial.angle1 = atof(argv[i]);
            pivotRadPerSecInitial = getAngVel(initial.angle1);
            initial.dangle2 = pivotRadPerSecInitial;
        }
        else if (strcmp(argv[i], "-angvel") == 0 && i + 1 < argc)
        {
            i++;
            pivotRadPerSecInitial = atof(argv[i]);
            initial.dangle2 = pivotRadPerSecInitial;
        }
    }

    double expectedAngle1 = -asin((pivotRadius + bobRadius) / param.length);

    testRoutine(initial, param, step);
}
