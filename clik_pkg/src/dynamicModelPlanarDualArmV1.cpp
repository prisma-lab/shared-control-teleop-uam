#include "dynamicModelPlanarDualArmV1.hpp"

#define mC 1
#define mA 0.5
#define LC 1
#define off 0.14
#define iCxx 0.2
#define iAxx 0.1
#define lCz 0.5
#define lAz 0.25
#define fvC 0.1
#define fvA 0.1
#define g0 9.8

// #define mC 1.917
// #define mA 0.799
// #define LC 1
// #define off 0.14
// #define iCxx 1.0927
// #define iAxx 0.017011890024962
// #define lCz 0.9992
// #define lAz 0.1454
// #define fvC 1.311
// #define fvA 0.05
// #define g0 9.8


void B_row1(Eigen::MatrixXd& B, Eigen::VectorXd qVec) {

  //double mC = 1;
  // double mA = 0.5;
  // double LC = 1;
  // double off = 0.14;
  // double iCxx = 0.2;
  // double iAxx = 0.1;
  // double lCz = 0.5;
  // double lAz = 0.25;
  
  double qC = qVec(0);
  double qAL = qVec(1);
  double qAR = qVec(2);

  double t2 = cos(qAL);
  double t3 = cos(qAR);
  double t4 = sin(qAL);
  double t5 = sin(qAR);
  double t6 = lAz*lAz;
  double t7 = mA*t6;
  B(0,0) = iAxx*2.0+iCxx+t7*2.0+(LC*LC)*mA*2.0+(lCz*lCz)*mC+mA*(off*off)*2.0+lAz*mA*off*t4*2.0-lAz*mA*off*t5*2.0+LC*lAz*mA*t2*2.0+LC*lAz*mA*t3*2.0;
  B(0,1) = iAxx+t7+lAz*mA*(LC*t2+off*t4);
  B(0,2) = iAxx+t7+lAz*mA*(LC*t3-off*t5);

}



void B_row2(Eigen::MatrixXd& B, Eigen::VectorXd qVec) {

  // double mC = 1;
  // double mA = 0.5;
  // double LC = 1;
  // double off = 0.14;
  // double iCxx = 0.2;
  // double iAxx = 0.1;
  // double lCz = 0.5;
  // double lAz = 0.25;

  double qC = qVec(0);
  double qAL = qVec(1);
  double qAR = qVec(2);

  double t2 = lAz*lAz;
  double t3 = mA*t2;
  B(1,0) = iAxx+t3+lAz*mA*(LC*cos(qAL)+off*sin(qAL));
  B(1,1) = iAxx+t3;
  B(1,2) = 0;

}

void B_row3(Eigen::MatrixXd& B, Eigen::VectorXd qVec) {

  // double mC = 1;
  // double mA = 0.5;
  // double LC = 1;
  // double off = 0.14;
  // double iCxx = 0.2;
  // double iAxx = 0.1;
  // double lCz = 0.5;
  // double lAz = 0.25;

  double qC = qVec(0);
  double qAL = qVec(1);
  double qAR = qVec(2);

  double t2 = lAz*lAz;
  double t3 = mA*t2;
  B(2,0) = iAxx+t3+lAz*mA*(LC*cos(qAR)-off*sin(qAR));
  B(2,1) = 0;
  B(2,2) = iAxx+t3;

}

void n_row1(Eigen::VectorXd& n, Eigen::VectorXd qVec, Eigen::VectorXd qDotVec) {

  // double mC = 1;
  // double mA = 0.5;
  // double LC = 1;
  // double off = 0.14;
  // double iCxx = 0.2;
  // double iAxx = 0.1;
  // double lCz = 0.5;
  // double lAz = 0.25;
  // double fvC = 0.1;
  // double fvA = 0.1;
  // double g0 = 9.8;

  double qC = qVec(0);
  double qAL = qVec(1);
  double qAR = qVec(2);
  double qCd = qDotVec(0);
  double qALd = qDotVec(1);
  double qARd = qDotVec(2);

  double t2 = cos(qAL);
  double t3 = cos(qAR);
  double t4 = sin(qAL);
  double t5 = sin(qAR);
  double t6 = sin(qC);
  double t7 = qALd*qALd;
  double t8 = qARd*qARd;
  n(0) = fvC*qCd+g0*lCz*mC*t6+g0*lAz*mA*sin(qAL+qC)+g0*lAz*mA*sin(qAR+qC)+LC*g0*mA*t6*2.0-LC*lAz*mA*t4*t7-LC*lAz*mA*t5*t8+lAz*mA*off*t2*t7-lAz*mA*off*t3*t8-LC*lAz*mA*qALd*qCd*t4*2.0-LC*lAz*mA*qARd*qCd*t5*2.0+lAz*mA*off*qALd*qCd*t2*2.0-lAz*mA*off*qARd*qCd*t3*2.0;

}

void n_row2(Eigen::VectorXd& n, Eigen::VectorXd qVec, Eigen::VectorXd qDotVec) {

  // double mC = 1;
  // double mA = 0.5;
  // double LC = 1;
  // double off = 0.14;
  // double iCxx = 0.2;
  // double iAxx = 0.1;
  // double lCz = 0.5;
  // double lAz = 0.25;
  // double fvC = 0.1;
  // double fvA = 0.1;
  // double g0 = 9.8;

  double qC = qVec(0);
  double qAL = qVec(1);
  double qAR = qVec(2);
  double qCd = qDotVec(0);
  double qALd = qDotVec(1);
  double qARd = qDotVec(2);

  double t2 = qCd*qCd;
  n(1) = fvA*qALd+g0*lAz*mA*sin(qAL+qC)+LC*lAz*mA*t2*sin(qAL)-lAz*mA*off*t2*cos(qAL);

}


void n_row3(Eigen::VectorXd& n, Eigen::VectorXd qVec, Eigen::VectorXd qDotVec) {

  // double mC = 1;
  // double mA = 0.5;
  // double LC = 1;
  // double off = 0.14;
  // double iCxx = 0.2;
  // double iAxx = 0.1;
  // double lCz = 0.5;
  // double lAz = 0.25;
  // double fvC = 0.1;
  // double fvA = 0.1;
  // double g0 = 9.8;

  double qC = qVec(0);
  double qAL = qVec(1);
  double qAR = qVec(2);
  double qCd = qDotVec(0);
  double qALd = qDotVec(1);
  double qARd = qDotVec(2);

  double t2 = qCd*qCd;
  n(2) = fvA*qARd+g0*lAz*mA*sin(qAR+qC)+LC*lAz*mA*t2*sin(qAR)+lAz*mA*off*t2*cos(qAR);

}


