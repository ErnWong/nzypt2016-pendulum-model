//http://www.zn.dmef.put.poznan.pl/content/009/009.pdf

#include <cmath>
#include <cstdio>

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
    double damp1;
    double damp2;
    Fn accelx;
    Fn accely;
    Fn accelz;
    Fn windx;
    Fn windz;
    Fn length;
    Fn dlength;
};

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
    double b1 = param.damp1;
    double b2 = param.damp2;
    double l = param.length(t);
    double dl = param.dlength(t);
    double ax = param.accelx(t);
    double ay = param.accely(t);
    double az = param.accelz(t);
    double fx = param.windx(t);
    double fz = param.windz(t);

    double c1 = cos(th1);
    double c2 = cos(th2);
    double s1 = sin(th1);
    double s2 = sin(th2);

    double n1 = -2.0*dth1*dth2*c2*s2*m*l*l;
    n1 += m*l*c2 * (2.0*dth1*dl*c2 + ax*c1 - ay*s1 + g*s1);
    n1 += b1*dth1;

    double n2 = dth1*dth1*c2*s2*m*l*l;
    n2 += m*l*(2.0*dth2*dl + g*c1*s2 - az*c2 - ay*c1*s2 - ax*s1*s2) + b2*dth2;

    double m1 = m*l*l*c2*c2;
    double m2 = m*l*l;

    double q1 = -fx*l*c1*c2;
    double q2 = fx*l*s1*s2 + fz*l*c2;

    d.dangle1 = (q1 - n1 - b1 * dth1) / m1;
    d.dangle2 = (q2 - n2 - b2 * dth2) / m2;
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

double accelPivotX(double time)
{
    return 10.5*sin(10.0*time);
}
double accelPivotY(double time)
{
    return 0.0;
}
double accelPivotZ(double time)
{
    return 10.5*cos(10.0*time);
}

double windX(double time)
{
    return 0.0;
}
double windZ(double time)
{
    return 0.0;
}

double stringLength(double time)
{
    return 1.0;
}
double dstringLength(double time)
{
    return 0.0;
}

int main()
{
    Param param;
    param.mass = 1.0;
    param.g = 9.81;
    param.damp1 = 0.0;
    param.damp2 = 0.0;
    param.accelx = accelPivotX;
    param.accely = accelPivotY;
    param.accelz = accelPivotZ;
    param.windx = windX;
    param.windz = windZ;
    param.length = stringLength;
    param.dlength = dstringLength;

    State state;
    state.angle1 = 0.1;
    state.angle2 = 0.0;
    state.dangle1 = 0.0;
    state.dangle2 = 0.0;

    double time = 0.0;
    double step = 0.001;

    int counter = 0;

    for (time = 0.0; true; time += step)
    {
        counter++;
        state = rk4(state, param, step, time);
        if (counter > 10)
        {
            counter = 0;
            double l = param.length(time);
            double x = l * sin(state.angle1) * cos(state.angle2);
            double y = l * sin(state.angle1) * sin(state.angle2);
            double z = l * cos(state.angle1);
            printf("%f %f %f\n", x, y, z);
        }
    }
}
