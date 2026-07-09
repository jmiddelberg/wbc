#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "robot_models/pinocchio/RobotModelPinocchio.hpp"
#include "scenes/acceleration_reduced_tsid/AccelerationSceneReducedTSID.hpp"
#include "solvers/qpoases/QPOasesSolver.hpp"
#include "tasks/SpatialAccelerationTask.hpp"

using namespace std;
using namespace wbc;

BOOST_AUTO_TEST_CASE(simple_test){

    /**
     * Check if the WBC scene computes the correct result, i.e., if the reference spatial acceleration matches the solver output, back-projected to Cartesian space
     */

    // Configure Robot model
    shared_ptr<RobotModelPinocchio> robot_model = make_shared<RobotModelPinocchio>();
    RobotModelConfig config;
    config.file_or_string = "../../../../../models/rh5/urdf/rh5_legs.urdf";
    config.floating_base = true;
    config.contact_points = {types::Contact("FL_SupportCenter",1,0.6,0.2,0.08), types::Contact("FR_SupportCenter",1,0.6,0.2,0.08)};
    BOOST_CHECK_EQUAL(robot_model->configure(config), true);

    types::JointState joint_state;
    joint_state.resize(robot_model->na());
    joint_state.position << 0,0,-0.35,0.64,0,-0.27,  0,0,-0.35,0.64,0,-0.27;
    joint_state.velocity.setZero();
    joint_state.acceleration.setZero();

    types::RigidBodyState rbs;
    rbs.pose.position = Eigen::Vector3d(-0.175,0,0.876);
    rbs.pose.orientation.setIdentity();
    rbs.twist.setZero();
    rbs.acceleration.setZero();

    BOOST_CHECK_NO_THROW(robot_model->update(joint_state.position,
                                             joint_state.velocity,
                                             joint_state.acceleration,
                                             rbs.pose,
                                             rbs.twist,
                                             rbs.acceleration));

    // Configure Solver
    QPSolverPtr solver = std::make_shared<QPOASESSolver>();
    dynamic_pointer_cast<QPOASESSolver>(solver)->setMaxNoWSR(1000);
    qpOASES::Options options = dynamic_pointer_cast<QPOASESSolver>(solver)->getOptions();
    options.printLevel = qpOASES::PL_NONE;
    dynamic_pointer_cast<QPOASESSolver>(solver)->setOptions(options);

    // Configure scene
    SpatialAccelerationTaskPtr cart_task;
    cart_task = make_shared<SpatialAccelerationTask>(TaskConfig("cart_pos_ctrl",0,Eigen::VectorXd::Ones(6),1),
                                                       robot_model,
                                                       "RH5_Root_Link");
    AccelerationSceneReducedTSID wbc_scene(robot_model, solver, 1e-3);
    BOOST_CHECK_EQUAL(wbc_scene.configure({cart_task}), true);

    // Set random Reference
    types::RigidBodyState ref;
    srand (time(NULL));
    ref.acceleration.linear = Eigen::Vector3d(((double)rand())/RAND_MAX, ((double)rand())/RAND_MAX, ((double)rand())/RAND_MAX);
    ref.acceleration.angular = Eigen::Vector3d(((double)rand())/RAND_MAX, ((double)rand())/RAND_MAX, ((double)rand())/RAND_MAX);
    BOOST_CHECK_NO_THROW(cart_task->setReference(ref.acceleration));

    // Solve
    HierarchicalQP qp;
    BOOST_CHECK_NO_THROW(qp=wbc_scene.update());
    types::JointCommand solver_output = wbc_scene.solve(qp);
    Eigen::VectorXd solver_output_raw = wbc_scene.getSolverOutputRaw();

    // Check
    uint nj = robot_model->nj();

    Eigen::VectorXd qdd = solver_output_raw.head(nj);
    Eigen::VectorXd tau = solver_output.effort;
    vector<types::Wrench> contact_wrenches = wbc_scene.getContactWrenches();

    Eigen::VectorXd y_solution = robot_model->spaceJacobian(cart_task->tipFrame())*qdd + robot_model->spatialAccelerationBias(cart_task->tipFrame()).vector6d();

    for(int i = 0; i < 3; i++){
        BOOST_CHECK(fabs(ref.acceleration.linear[i] - y_solution[i]) < 1e-3);
        BOOST_CHECK(fabs(ref.acceleration.angular[i] - y_solution[i]) < 1e3);
    }

    // Check if torques respect equation of motions
    const auto& contacts = robot_model->getContacts();
    Eigen::VectorXd eq_motion_left = robot_model->jointSpaceInertiaMatrix() * qdd + robot_model->biasForces();
    Eigen::VectorXd eq_motion_right(nj);
    eq_motion_right = robot_model->selectionMatrix().transpose()*tau;
    for(uint i=0; i < contacts.size(); ++i)
        eq_motion_right += robot_model->spaceJacobian(contacts[i].frame_id).topRows(3).transpose() * contact_wrenches[i].force;


    BOOST_CHECK((eq_motion_left - eq_motion_right).cwiseAbs().maxCoeff() < 1e-3);
}

BOOST_AUTO_TEST_CASE(delta_penalty){

    /**
     * Check the penalty on the difference between consecutive solver outputs: A scene with delta penalty should show a
     * smaller jump between consecutive solutions than a scene without, when the task reference changes abruptly
     */

    // Configure Robot model
    shared_ptr<RobotModelPinocchio> robot_model = make_shared<RobotModelPinocchio>();
    RobotModelConfig config;
    config.file_or_string = "../../../../../models/rh5/urdf/rh5_legs.urdf";
    config.floating_base = true;
    config.contact_points = {types::Contact("FL_SupportCenter",1,0.6,0.2,0.08), types::Contact("FR_SupportCenter",1,0.6,0.2,0.08)};
    BOOST_CHECK_EQUAL(robot_model->configure(config), true);

    types::JointState joint_state;
    joint_state.resize(robot_model->na());
    joint_state.position << 0,0,-0.35,0.64,0,-0.27,  0,0,-0.35,0.64,0,-0.27;
    joint_state.velocity.setZero();
    joint_state.acceleration.setZero();

    types::RigidBodyState rbs;
    rbs.pose.position = Eigen::Vector3d(-0.175,0,0.876);
    rbs.pose.orientation.setIdentity();
    rbs.twist.setZero();
    rbs.acceleration.setZero();

    robot_model->update(joint_state.position, joint_state.velocity, joint_state.acceleration,
                        rbs.pose, rbs.twist, rbs.acceleration);

    // Two identical scenes with separate solvers, one with delta penalty
    auto makeSolver = [](){
        QPSolverPtr solver = std::make_shared<QPOASESSolver>();
        dynamic_pointer_cast<QPOASESSolver>(solver)->setMaxNoWSR(1000);
        qpOASES::Options options = dynamic_pointer_cast<QPOASESSolver>(solver)->getOptions();
        options.printLevel = qpOASES::PL_NONE;
        dynamic_pointer_cast<QPOASESSolver>(solver)->setOptions(options);
        return solver;
    };

    SpatialAccelerationTaskPtr cart_task;
    cart_task = make_shared<SpatialAccelerationTask>(TaskConfig("cart_pos_ctrl",0,Eigen::VectorXd::Ones(6),1),
                                                     robot_model,
                                                     "RH5_Root_Link");
    AccelerationSceneReducedTSID scene_plain(robot_model, makeSolver(), 1e-3);
    AccelerationSceneReducedTSID scene_smooth(robot_model, makeSolver(), 1e-3);
    BOOST_CHECK_EQUAL(scene_plain.configure({cart_task}), true);
    BOOST_CHECK_EQUAL(scene_smooth.configure({cart_task}), true);
    scene_smooth.setAccelerationDeltaPenalty(1e3);
    scene_smooth.setContactWrenchDeltaPenalty(1e3);

    // First reference
    types::SpatialAcceleration ref;
    ref.linear = Eigen::Vector3d(0.1,0.2,0.3);
    ref.angular = Eigen::Vector3d(0.05,-0.1,0.1);
    cart_task->setReference(ref);

    scene_plain.solve(scene_plain.update());
    Eigen::VectorXd x0_plain = scene_plain.getSolverOutputRaw();
    scene_smooth.solve(scene_smooth.update());
    Eigen::VectorXd x0_smooth = scene_smooth.getSolverOutputRaw();

    // First solution is not affected by the delta penalty (no previous solver output exists)
    BOOST_CHECK((x0_plain - x0_smooth).norm() < 1e-6);

    // Second, strongly different reference
    ref.linear = Eigen::Vector3d(-2.0,1.5,-1.0);
    ref.angular = Eigen::Vector3d(1.0,0.5,-1.5);
    cart_task->setReference(ref);

    scene_plain.solve(scene_plain.update());
    Eigen::VectorXd x1_plain = scene_plain.getSolverOutputRaw();
    scene_smooth.solve(scene_smooth.update());
    Eigen::VectorXd x1_smooth = scene_smooth.getSolverOutputRaw();

    // With delta penalty, the solution changes much less between consecutive cycles
    double jump_plain = (x1_plain - x0_plain).norm();
    double jump_smooth = (x1_smooth - x0_smooth).norm();
    BOOST_CHECK(jump_plain > 1e-3); // the reference change actually produces a jump
    BOOST_CHECK(jump_smooth < 0.1 * jump_plain);
}

BOOST_AUTO_TEST_CASE(friction_cone_slack){

    /**
     * Check the softened friction cone constraint (one slack variable per contact): When the cone is inactive, the slack
     * variables should be (close to) zero and the solution identical to the hard constraint. When the reference acceleration
     * pushes the contact wrenches against the cone, the softened scene should track the reference better than the hard one,
     * at the cost of a small cone violation
     */

    // Configure Robot model with a low friction coefficient, so that the friction cone is easily activated
    shared_ptr<RobotModelPinocchio> robot_model = make_shared<RobotModelPinocchio>();
    RobotModelConfig config;
    config.file_or_string = "../../../../../models/rh5/urdf/rh5_legs.urdf";
    config.floating_base = true;
    config.contact_points = {types::Contact("FL_SupportCenter",1,0.05,0.2,0.08), types::Contact("FR_SupportCenter",1,0.05,0.2,0.08)};
    BOOST_CHECK_EQUAL(robot_model->configure(config), true);

    types::JointState joint_state;
    joint_state.resize(robot_model->na());
    joint_state.position << 0,0,-0.35,0.64,0,-0.27,  0,0,-0.35,0.64,0,-0.27;
    joint_state.velocity.setZero();
    joint_state.acceleration.setZero();

    types::RigidBodyState rbs;
    rbs.pose.position = Eigen::Vector3d(-0.175,0,0.876);
    rbs.pose.orientation.setIdentity();
    rbs.twist.setZero();
    rbs.acceleration.setZero();

    robot_model->update(joint_state.position, joint_state.velocity, joint_state.acceleration,
                        rbs.pose, rbs.twist, rbs.acceleration);

    auto makeSolver = [](){
        QPSolverPtr solver = std::make_shared<QPOASESSolver>();
        dynamic_pointer_cast<QPOASESSolver>(solver)->setMaxNoWSR(1000);
        qpOASES::Options options = dynamic_pointer_cast<QPOASESSolver>(solver)->getOptions();
        options.printLevel = qpOASES::PL_NONE;
        dynamic_pointer_cast<QPOASESSolver>(solver)->setOptions(options);
        return solver;
    };

    // Two identical scenes with surface contacts (dim_contact == 6), one with softened friction cone
    SpatialAccelerationTaskPtr cart_task;
    cart_task = make_shared<SpatialAccelerationTask>(TaskConfig("cart_pos_ctrl",0,Eigen::VectorXd::Ones(6),1),
                                                     robot_model,
                                                     "RH5_Root_Link");
    AccelerationSceneReducedTSID scene_hard(robot_model, makeSolver(), 1e-3, 6);
    AccelerationSceneReducedTSID scene_soft(robot_model, makeSolver(), 1e-3, 6);
    BOOST_CHECK_EQUAL(scene_hard.configure({cart_task}), true);
    BOOST_CHECK_EQUAL(scene_soft.configure({cart_task}), true);
    scene_soft.setFrictionConeSlackPenalty(1e-3);

    uint nj = robot_model->nj();
    uint nc = robot_model->nc();

    auto trackingError = [&](AccelerationSceneReducedTSID& scene, const types::SpatialAcceleration& ref){
        scene.solve(scene.update());
        Eigen::VectorXd qdd = scene.getSolverOutputRaw().head(nj);
        Eigen::VectorXd y = robot_model->spaceJacobian(cart_task->tipFrame())*qdd +
                            robot_model->spatialAccelerationBias(cart_task->tipFrame()).vector6d();
        return (y - ref.vector6d()).norm();
    };

    // Mild reference: friction cone is inactive, slack variables should stay at zero
    types::SpatialAcceleration ref;
    ref.setZero();
    cart_task->setReference(ref);

    double err_hard = trackingError(scene_hard, ref);
    double err_soft = trackingError(scene_soft, ref);

    // QP of soft scene contains one additional variable per contact
    BOOST_CHECK_EQUAL(scene_hard.getSolverOutputRaw().size(), nj + 6*nc);
    BOOST_CHECK_EQUAL(scene_soft.getSolverOutputRaw().size(), nj + 6*nc + nc);

    Eigen::VectorXd slack = scene_soft.getSolverOutputRaw().tail(nc);
    BOOST_CHECK(slack.minCoeff() > -1e-6); // slack variables are positive
    BOOST_CHECK(slack.maxCoeff() < 1e-3);  // ... and (close to) zero, since the cone is inactive
    BOOST_CHECK(fabs(err_hard - err_soft) < 1e-3); // solution matches the hard constraint

    // Aggressive lateral reference: required tangential forces exceed the friction cone
    ref.linear = Eigen::Vector3d(3.0,0,0);
    ref.angular.setZero();
    cart_task->setReference(ref);

    err_hard = trackingError(scene_hard, ref);
    err_soft = trackingError(scene_soft, ref);
    slack = scene_soft.getSolverOutputRaw().tail(nc);

    BOOST_CHECK(err_hard > 1e-2);          // the hard cone constraint prevents tracking the reference
    BOOST_CHECK(slack.maxCoeff() > 1e-4);  // the softened cone is violated ...
    BOOST_CHECK(err_soft < 0.5 * err_hard);  // ... which allows better tracking of the reference
}
