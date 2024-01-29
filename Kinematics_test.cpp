#include<iostream>
#include<cmath>
#include<Eigen/Dense>
#include<vector>
#include "Kinematics_header.h"
#include "Constants_header.h"

using namespace std;
using namespace Eigen;

Matrix4d transformation_matrix ( double alpha, double a, double theta, double d){
    Matrix4d T_matrix;
    T_matrix << cos(theta), -sin(theta)*cos(alpha),  sin(theta)*cos(alpha), d*cos(theta),
                sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), d*sin(theta),
                         0,             sin(alpha),             cos(alpha),            d,
                         0,                      0,                      0,            1;
    return T_matrix;
}

vector<Matrix4d> calculateTransformationMatrix(VectorXd Q){
    Matrix <double, 7, 4> DH_param;

    DH_param << -pi/2, 0, Q[0], d1,
                 pi/2, 0, Q[1], 0,
                -pi/2, 0, Q[2], d3,
                 pi/2, 0, Q[3], 0,
                -pi/2, 0, Q[4], d5,
                 pi/2, 0, Q[5], 0,
                    0, 0, Q[6], d7;

    vector<Matrix4d> TMatArray(7);

    Matrix4d T_07 = Matrix4d::Identity();

    for (int i = 0; i <7; ++i){
        double alpha = DH_param(i, 0);
        double a = DH_param(i, 1);
        double theta = DH_param(i, 2);
        double d = DH_param(i, 3);

        Matrix4d TMat = transformation_matrix( alpha, a, theta, d);

        T_07 *= TMat;
        TMatArray[i] = TMat;
    }

    return TMatArray;
}

Matrix4d finalTransformationMatrix( vector<Eigen::Matrix4d> TMatArray){
    Matrix4d TMat = Matrix4d::Identity();
    for (int i = 0; i < 7; i++)
    {
        TMat *= TMatArray[i];
    }

    return TMat;
}

VectorXd calculateEndEffectorPoseAndPosition(Eigen::VectorXd Q){
    vector<Matrix4d> TMatArray = calculateTransformationMatrix(Q);
    Matrix4d T_07 = finalTransformationMatrix( TMatArray);
    Matrix3d R_07 = T_07.block<3,3>(0,0);
    Vector3d t_07 = T_07.block<3,1>(0,3);
    Vector3d Pose_act = R_07.eulerAngles(2,1,0);
    Vector3d Position_act = t_07;

    VectorXd X(6);

    X << Position_act, Pose_act;

    return X;
}

MatrixXd calculateJacobian (Eigen::VectorXd Q){
    MatrixXd J(6,7);

    vector<Matrix4d> TMatArray = calculateTransformationMatrix(Q);

    Matrix4d T_07 = finalTransformationMatrix( TMatArray);
    Matrix3d R_07 = T_07.block<3,3>(0,0);
    Vector3d t_07 = T_07.block<3,1>(0,3);
    Matrix4d T_0i = Matrix4d::Identity();
    for (int i = 0; i < 7; i++)
    {
        Matrix4d T_i = TMatArray[i];
        T_0i *= T_i;

        Matrix3d R_0i = T_0i.block<3,3>(0,0);
        Vector3d t_0i = T_0i.block<3,1>(0, 3);

        Vector3d z_0i = R_0i.col(2);

        Vector3d d_calc = t_07 - t_0i;
        Vector3d J_top = z_0i.cross(d_calc);


        J.col(i).segment(0, 3) = J_top;
        J.col(i).segment(3, 3) = z_0i;
    }
    return J;
}

int main()
{
    kinematics::joint_angles Q;
    VectorXd Q_original = VectorXd::Zero(7);
    VectorXd Q_vals = VectorXd::Zero(7);
    cout << "Initial Q values" << endl << Q_vals << endl;
    for (int i = 0; i < 7; i++)
    {
        Q_original[i] = Q.q[i];
        Q_vals[i] = Q.q[i];
    }


    vector<Matrix4d> TMatArray = calculateTransformationMatrix(Q_vals);

    Matrix4d T_07 = finalTransformationMatrix( TMatArray);
    cout << "Transformation Matrix" << endl << T_07 << endl;

    VectorXd X(6);
    X = calculateEndEffectorPoseAndPosition(Q_vals);
    cout << "End effector position" << endl << X << endl;

    MatrixXd J = calculateJacobian(Q_vals);
    cout << "Jacobian Matrix" << endl << J << endl;

    MatrixXd J_inv = J.completeOrthogonalDecomposition().pseudoInverse();

    JacobiSVD<MatrixXd> svd(J);
    VectorXd singularValues = svd.singularValues();

    double conditionNumber = singularValues.maxCoeff() / singularValues.minCoeff();
    cout << " Condition Number is " << conditionNumber << endl;

    kinematics::desired_pos EndEffector;

    Vector3d XYZ_des;
    XYZ_des(0)= 0.384483;
    XYZ_des(1)= 1.21195;
    XYZ_des(2)= 0.471398;

    Vector3d Angle_des;
    Angle_des(0)= -0.0609516;
    Angle_des(1)= -0.143048;
    Angle_des(2)=  0.0913827;

    EndEffector.Desired_Position = XYZ_des ;
    EndEffector.Desired_Pose = Angle_des;

    VectorXd Xd(6);

    Xd.segment<3>(0) = EndEffector.Desired_Position;
    Xd.segment<3>(3) = EndEffector.Desired_Pose.transpose();
    cout << "Desired End effector position" << endl << Xd << endl;

    VectorXd e(6);
    e = Xd - X;
    cout << "initial error" << endl<< e << endl;
    VectorXd Q0 = VectorXd::Zero(7);
    VectorXd delta_q = VectorXd::Zero(7);
    int iter = 0;
    VectorXd joint_limits_lower(7);
    VectorXd joint_limits_upper(7);
    joint_limits_lower << -pi, -(2*pi/3), -pi, -(5*pi/3), -pi, -pi, -pi;
    joint_limits_upper <<  pi,  (2*pi/3),  pi,  (5*pi/3),  pi,  pi,  pi;;

    cout << "lower joint lim" << endl << joint_limits_lower << endl;
    cout << "upper joint lim" << endl << joint_limits_upper << endl;
    while ((abs(e.array()) > 0.001).any() && iter < 50)
    {

        J = calculateJacobian(Q0);
        J_inv = J.completeOrthogonalDecomposition().pseudoInverse();
        delta_q = Q0 + J_inv*e;

        for (int i = 0; i < 7; ++i)
        {
            if (Q0[i] < joint_limits_lower[i])
                Q0[i] = joint_limits_lower[i];
            else if (Q0[i] > joint_limits_upper[i])
                Q0[i] = joint_limits_upper[i];
        }
        VectorXd X_new(6);
        X_new = calculateEndEffectorPoseAndPosition(delta_q);
        e = Xd - X_new;

        Q0 = delta_q;

        iter++;
        cout << "Iteration: " << iter <<endl;

    }
    VectorXd Error;

    Error = Q_original - Q0;
    cout << "Final Q is" << endl << Q0 << endl;
    cout<< "Final Error is" << endl << Error << endl;
    return 0;
}

