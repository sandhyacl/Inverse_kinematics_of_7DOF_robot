#ifndef Kinematics_7dof
#define Kinematics_7dof

#include<Eigen/Dense>
#include<vector>

class kinematics
{

public:
    struct joint_angles
    {
        //joint_angles(): q{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}{}
        joint_angles(): q{ 0.7854,0.7854, 0.7854,0.7854,0.7854,0.7854,0.7854}{}
        joint_angles(double q1, double q2, double q3, double q4, double q5, double q6, double q7) : q{q1, q2, q3, q4, q5, q6, q7}
        {}

        double q[7];
    };
    struct joint_limits
    {
        joint_limits(): q_lim{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}{}

        joint_limits(double q1_lim = 180.0, double q2_lim = 120.0, double q3_lim = 180.0, double q4_lim = 150.0, double q5_lim = 180.0, double q6_lim = 180.0, double q7_lim = 180.0)
                    : q_lim{q1_lim, q2_lim, q3_lim, q4_lim, q5_lim, q6_lim, q7_lim}
        {}

        double q_lim[7];
    };

    struct desired_pos
    {
        Eigen::Vector3d Desired_Position;
        Eigen::Vector3d Desired_Pose;
    };


    static Eigen::Matrix4d transformation_matrix ( double alpha, double a, double theta, double d);

    static Eigen::MatrixXd calculateJacobian (Eigen::VectorXd Q);

    static std::vector<Eigen::Matrix4d> calculateTransformationMatrix( Eigen::VectorXd Q);

    static Eigen::Matrix4d finalTransformationMatrix( std::vector<Eigen::Matrix4d> TMatArray);

    static Eigen::VectorXd calculateEndEffectorPoseAndPosition(Eigen::VectorXd Q);

};

#endif

