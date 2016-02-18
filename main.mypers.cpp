//http://www.zn.dmef.put.poznan.pl/content/009/009.pdf

#include <cmath>
#include <cstdio>

#define PLOT 0
#define IS_TEST 1

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

typedef double (*Fn)(double t);

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

double diff(Fn fn, double x)
{
    double dx = 1.0e-6;
    double x1 = x - dx;
    double x2 = x + dx;
    double y1 = fn(x1);
    double y2 = fn(x2);
    return (y2 - y1) / (x2 - x1);
}

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

double bobRadius = 0.1503;
double pivotRadius = 0.45;
//double pivotRPM = 61.97;
//double pivotRadPerSec = 8; // pivotRPM * 3.141592653 * 2 / 60.0;
double pivotRadPerSecInitial = 7.0;
double pivotRadPerSecChange = 0.0;

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
    // return pivotRadius * cos(pivotRadPerSec * time);
    // return 1.0 / (1.0 + pow(2.718, -10.0 * (time - 2.0))) + 1.0 / (1.0 + pow(2.718, 10.0 * (time - 10.0)));
    // return time < 10.0? 0.0 : 5.0 * (time - 10.0) * (time - 10.0);
}
double pivotY(double time)
{
    return pivotRadius * sin(pivotPhase(time));
    // return pivotRadius * sin(pivotRadPerSec * time);
    // return 0.0;
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
    //return time < 10.0? 0.0 : 10.0;
    //return -pivotRadius * pivotRadPerSec * cos(pivotRadPerSec * time);
}
double accelPivotY(double time)
{
    return diff(pivotVelY, time);
    //return 0;
    //return -pivotRadius * pivotRadPerSec * sin(pivotRadPerSec * time);
}
double accelPivotZ(double time)
{
    return diff(pivotVelZ, time);
    return 0.0;
}

int main()
{
    Param param;
    param.mass = 1.0;
    param.g = 9.81;
    param.length = 1.0;
    param.damp1 = 0.1 * 10;
    param.damp2 = 0.1 * 10;
    param.accelx = accelPivotX;
    param.accely = accelPivotY;
    param.accelz = accelPivotZ;

    double expectedAngle1 = -asin((pivotRadius + bobRadius) / param.length);

    State initial;
    initial.angle1 = expectedAngle1;//-45 / 180.0 * 3.141592653;
    initial.angle2 = 0.0;
    initial.dangle1 = 0.0;
    initial.dangle2 = pivotRadPerSec(0); //5.0;

    double time = 0.0;
    double step = 0.00002;
    int counter = 0;
    int sameCounter = 0;
    double angle1Prev = 0.0;

    if (IS_TEST)
    {
        for (double angVel = 8; angVel > 4; angVel -= 0.1)
        {
            printf("Angvel = %f\n", angVel);
            pivotRadPerSecInitial = angVel;
            double expectedAngle1 = -asin((pivotRadius + bobRadius) / param.length); 
            initial.angle1 = expectedAngle1;
            State state = initial;
            sameCounter = 0;
            while (sameCounter < 10000)
            {
                counter++;
                state = rk4(state, param, step, time);
                double difference = state.angle1 - angle1Prev;
                if (counter >= 10000)
                {
                    //printf("angle1 = %f\n", state.angle1);
                    counter = 0;
                }
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
            printf("%f %f\n", angVel, angle1Prev);
        }
    }
    else
    {
        State state = initial;
        for (time = 0.0; true; time += step)
        {
            counter++;
            state = rk4(state, param, step, time);
            double difference = state.angle1 - angle1Prev;
            if (difference < 0) difference *= -1;
            if (difference < 0.000001)
            {
                sameCounter++;
                if (sameCounter > 10000)
                {
                    printf("DONE!\n");
                    return 0;
                }
            }
            else
            {
                sameCounter = 0;
                angle1Prev = state.angle1;
            }
            if (counter >= 10000)
            {
                counter = 0;
                if (!PLOT)
                {
                    printf("%.17g\n", state.angle1);
                }
                else
                {
                    double px = pivotX(time);
                    double py = pivotY(time);
                    double pz = pivotZ(time);
                    double l = param.length;
                    double x = px + l * sin(state.angle1) * cos(state.angle2);
                    double y = py + l * sin(state.angle1) * sin(state.angle2);
                    double z = pz + l * cos(state.angle1);
                    printf("%f %f %f %f %f %f %f\n", pivotRadPerSec(time), px, py, pz, x, y, z);
                }
            }
        }
    }
}
