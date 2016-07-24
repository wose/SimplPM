#include "norad.h"

#include "gmock/gmock.h"

using namespace ::testing;

class SomeTestVessels : public Test
{
public:
    const double epsilon = 3e-2;

    double sgpPositions[5 * 3 ] = {
        2328.96594238,  -5995.21600342,   1719.97894287,
        2456.00610352,  -6071.94232177,   1222.95977784,
        2567.39477539,  -6112.49725342,    713.97710419,
        2663.03179932,  -6115.37414551,    195.73919105,
        2742.85470581,  -6079.13580322,   -328.86091614,
    };

    double sgp4Positions[5 * 3] = {
        2328.97048951,  -5995.22076416,   1719.97067261,
        2456.10705566,  -6071.93853760,   1222.89727783,
        2567.56195068,  -6112.50384522,    713.96397400,
        2663.09078980,  -6115.48229980,    196.39640427,
        2742.55133057,  -6079.67144775,   -326.38095856,
    };

    double sgp8Positions[5 * 3] = {
        2328.87265015,  -5995.21289063,   1720.04884338,
        2456.04577637,  -6071.90490722,   1222.84086609,
        2567.68383789,  -6112.40881348,    713.29282379,
        2663.49508667,  -6115.18182373,    194.62816810,
        2743.29238892,  -6078.90783691,   -329.73434067,
    };

    double sdp4Positions[5 * 3] = {
        7473.37066650,    428.95261765,   5828.74786377,
        -3305.22537232,  32410.86328125, -24697.1767581,
        14271.28759766,  24110.46411133,  -4725.76837158,
        -9990.05883789,  22717.35522461, -23616.89062501,
        9787.86975097,  33753.34667969, -15030.81176758,
    };

    double sdp8Positions[5 * 3] = {
        7469.47631836,    415.99390792,   5829.64318848,
        -3337.38992310,  32351.39086914, -24658.63037109,
        14226.54333496,  24236.08740234,  -4856.19744873,
        -10151.59838867,  22223.69848633, -23392.39770508,
        9420.08203125,  33847.21875000, -15391.06469727,
    };

    double sgpVelocities[5 * 3] = {
        2.91110113,     -0.98164053,      -7.09049922,
        2.67852119,     -0.44705850,      -7.22800565,
        2.43952477,      0.09884824,      -7.31899641,
        2.19531813,      0.65333930,      -7.36169147,
        1.94707947,      1.21346101,      -7.35499924,
    };

    double sgp4Velocities[5 * 3] = {
        2.91207230,     -0.98341546,      -7.09081703,
        2.67938992,     -0.44829041,      -7.22879231,
        2.44024599,      0.09810869,      -7.31995916,
        2.19611958,      0.65241995,      -7.36282432,
        1.94850229,      1.21106251,      -7.35619372,
    };

    double sgp8Velocities[5 * 3] = {
        2.91210661,     -0.98353850,      -7.09081554,
        2.67936245,     -0.44820847,      -7.22888553,
        2.43992555,      0.09893919,      -7.32018769,
        2.19525236,      0.65453661,      -7.36308974,
        1.94680957,      1.21500109,      -7.35625595,
    };

    double sdp4Velocities[5 * 3] = {
        5.10715413,      6.44468284,      -0.18613096,
        -1.30113538,     -1.15131518,      -0.28333528,
        -0.32050445,      2.67984074,      -2.08405289,
        -1.01667246,     -2.29026759,       0.72892364,
        -1.09425066,      0.92358845,      -1.52230928,
    };

    double sdp8Velocities[5 * 3] = {
        5.11402285,      6.44403201,      -0.18296110,
        -1.30200730,     -1.15603013,      -0.28164955,
        -0.33951668,      2.65315416,      -2.08114153,
        -1.00112480,     -2.33532837,       0.76987664,
        -1.11986055,      0.85410149,      -1.49506933
    };

    std::vector<std::string> vessels = {
        "1 88888U          80275.98708465  .00073094  13844-3  66816-4 0    87",
        "2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518  1058",
        "1 11801U          80230.29629788  .01431103  00000-0  14311-1       2",
        "2 11801U 46.7916 230.4354 7318036  47.4722  10.4117  2.28537848     2",
    };

    double satParams[N_SAT_PARAMS];
    double position[3];
    double velocity[3];

    tle_t tle;
    void checkResult(double* acutalPosition, double* targetPosition,
                     double* actualVelocity, double* targetVelocity)
    {
        actualVelocity[0] /= 60;
        actualVelocity[1] /= 60;
        actualVelocity[2] /= 60;

        EXPECT_THAT(acutalPosition[0], DoubleNear(*targetPosition++, epsilon));
        EXPECT_THAT(acutalPosition[1], DoubleNear(*targetPosition++, epsilon));
        EXPECT_THAT(acutalPosition[2], DoubleNear(*targetPosition++, epsilon));

        EXPECT_THAT(actualVelocity[0], DoubleNear(*targetVelocity++, epsilon));
        EXPECT_THAT(actualVelocity[1], DoubleNear(*targetVelocity++, epsilon));
        EXPECT_THAT(actualVelocity[2], DoubleNear(*targetVelocity++, epsilon));
    }
};

TEST_F(SomeTestVessels, OrbitCanBePredicted_SGP)
{
    double *targetPosition = sgpPositions;
    double *targetVelocity = sgpVelocities;

    ASSERT_THAT(parse_elements(vessels[0].c_str(), vessels[1].c_str(), &tle), Eq(0));
    SGP_init(satParams, &tle);

    for (int index = 0; index < 5; ++ index) {
        ASSERT_THAT(SGP(index * 360.0f, &tle, satParams, position, velocity), Eq(0));
        checkResult(position, &targetPosition[index * 3], velocity, &targetVelocity[index * 3]);
    }
}

TEST_F(SomeTestVessels, OrbitCanBePredicted_SGP4)
{
    double *targetPosition = sgp4Positions;
    double *targetVelocity = sgp4Velocities;

    ASSERT_THAT(parse_elements(vessels[0].c_str(), vessels[1].c_str(), &tle), Eq(0));
    SGP4_init(satParams, &tle);

    for (int index = 0; index < 5; ++ index) {
        ASSERT_THAT(SGP4(index * 360.0f, &tle, satParams, position, velocity), Eq(0));
        checkResult(position, &targetPosition[index * 3], velocity, &targetVelocity[index * 3]);
    }
}

TEST_F(SomeTestVessels, OrbitCanBePredicted_SGP8)
{
    double *targetPosition = sgp8Positions;
    double *targetVelocity = sgp8Velocities;

    ASSERT_THAT(parse_elements(vessels[0].c_str(), vessels[1].c_str(), &tle), Eq(0));
    SGP8_init(satParams, &tle);

    for (int index = 0; index < 5; ++ index) {
        ASSERT_THAT(SGP8(index * 360.0f, &tle, satParams, position, velocity), Eq(0));
        checkResult(position, &targetPosition[index * 3], velocity, &targetVelocity[index * 3]);
    }
}

TEST_F(SomeTestVessels, OrbitCanBePredicted_SDP4)
{
    double *targetPosition = sdp4Positions;
    double *targetVelocity = sdp4Velocities;

    ASSERT_THAT(parse_elements(vessels[2].c_str(), vessels[3].c_str(), &tle), Eq(0));
    SDP4_init(satParams, &tle);

    for (int index = 0; index < 5; ++ index) {
        ASSERT_THAT(SDP4(index * 360.0f, &tle, satParams, position, velocity), Eq(0));
        checkResult(position, &targetPosition[index * 3], velocity, &targetVelocity[index * 3]);
    }
}

TEST_F(SomeTestVessels, OrbitCanBePredicted_SDP8)
{
    double *targetPosition = sdp8Positions;
    double *targetVelocity = sdp8Velocities;

    ASSERT_THAT(parse_elements(vessels[2].c_str(), vessels[3].c_str(), &tle), Eq(0));
    SDP8_init(satParams, &tle);

    for (int index = 0; index < 5; ++ index) {
        ASSERT_THAT(SDP8(index * 360.0f, &tle, satParams, position, velocity), Eq(0));
        checkResult(position, &targetPosition[index * 3], velocity, &targetVelocity[index * 3]);
    }
}
