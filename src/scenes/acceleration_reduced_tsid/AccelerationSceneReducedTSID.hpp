#ifndef WBCACCELERATIONSCENEREDUCEDTSID_HPP
#define WBCACCELERATIONSCENEREDUCEDTSID_HPP

#include "../../core/Scene.hpp"
#include "../../types/Wrench.hpp"

namespace wbc{

class ContactsFrictionSurfaceConstraint;

/**
 * @brief Acceleration-based implementation of the WBC Scene. It sets up and solves the following problem:
 *  \f[
 *        \begin{array}{ccc}
 *        minimize &  \| \mathbf{J}_w\ddot{\mathbf{q}} - \dot{\mathbf{v}}_d + \dot{\mathbf{J}}\dot{\mathbf{q}}\|_2\\
 *        \mathbf{\ddot{q}},\mathbf{\tau},\mathbf{f} & & \\
 *           s.t.  & \mathbf{H}\mathbf{\ddot{q}} - \mathbf{S}^T\mathbf{\tau} - \mathbf{J}_c^T\mathbf{f} = -\mathbf{h} & \\
 *                 & \mathbf{J}_{c,i}\mathbf{\ddot{q}} = -\dot{\mathbf{J}}_{c,i}\dot{\mathbf{q}}, \, \forall i& \\
 *                 & \mathbf{\tau}_m \leq \mathbf{\tau} \leq \mathbf{\tau}_M& \\
 *        \end{array}
 *  \f]
 * \f$\ddot{\mathbf{q}}\f$ - Vector of robot joint accelerations<br>
 * \f$\mathbf{v}_{d}\f$ - Desired spatial accelerations of all tasks stacked in a vector<br>
 * \f$\mathbf{J}\f$ - Task Jacobians of all tasks stacked in a single matrix<br>
 * \f$\mathbf{J}_w = \mathbf{W}\mathbf{J}\f$ - Weighted task Jacobians<br>
 * \f$\mathbf{W}\f$ - Diagonal task weight matrix<br>
 * \f$\mathbf{H}\f$ - Joint space inertia matrix<br>
 * \f$\mathbf{S}\f$ - Selection matrix<br>
 * \f$\mathbf{\tau}\f$ - actuation forces/torques<br>
 * \f$\mathbf{h}\f$ - bias forces/torques<br>
 * \f$\mathbf{f}\f$ - external forces<br>
 * \f$\mathbf{J}_{c,i}\f$ - Contact Jacobian of i-th contact point<br>
 * \f$\dot{\mathbf{J}}\dot{\mathbf{q}}\f$ - Acceleration bias<br>
 * \f$\mathbf{\tau}_m,\mathbf{\tau}_M\f$ - Joint force/torque limits<br>
 *
 * The implementation is close to the task-space-inverse dynamics (TSID) method: https://andreadelprete.github.io/teaching/tsid/1_tsid_theory.pdf.
 * It computes the required joint space accelerations \f$\ddot{\mathbf{q}}\f$, torques \f$\mathbf{\tau}\f$ and contact wrenches \f$\mathbf{f}\f$, required to achieve the given task space
 * accelerations \f$\mathbf{v}_{d}\f$ under consideration of the equations of motion (eom), rigid contacts and joint force/torque limits. Note that onyl a single hierarchy level is allowed here,
 * prioritization can be achieved by assigning suitable task weights \f$\mathbf{W}\f$.
 */
class AccelerationSceneReducedTSID : public Scene{
protected:
    static SceneRegistry<AccelerationSceneReducedTSID> reg;

    Eigen::VectorXd robot_acc, solver_output_acc;
    std::vector<types::Contact> contacts;
    std::vector< TaskPtr > tasks;
    std::vector< ConstraintPtr > constraints;
    HierarchicalQP hqp;
    bool configured;
    types::JointCommand solver_output_joints;
    uint dim_contact;
    bool use_spatial_acc_bias;
    double acceleration_penalty;
    double contact_wrench_penalty;
    double acceleration_delta_penalty;
    double wrench_delta_penalty;
    double friction_cone_slack_penalty;
    std::shared_ptr<ContactsFrictionSurfaceConstraint> friction_surface_constraint;

    bool contactsHaveChanged(const std::vector<types::Contact>& old_contacts, const std::vector<types::Contact>& new_contacts){
        if(old_contacts.size() != new_contacts.size())
            return true;
        for(uint i = 0; i < old_contacts.size(); i++){
            if(old_contacts[i].active != new_contacts[i].active)
                return true;
        }
        return false;
    }

public:
    AccelerationSceneReducedTSID(RobotModelPtr robot_model, QPSolverPtr solver, const double dt, uint dim_contact = 3, bool use_spatial_acc_bias = true);
    virtual ~AccelerationSceneReducedTSID(){}

    /**
     * @brief Configure the WBC scene. Create tasks and sort them by priority given the task config
     * @param tasks Tasks used in optimization function. Size has to be > 0. All tasks have to be valid. See tasks and TaskConfig.hpp for more details.
     */
    virtual bool configure(const std::vector<TaskPtr> &tasks);

    /**
     * @brief Update the wbc scene and return the (updated) optimization problem
     * @return The (updated) optimization problem
     */
    virtual const HierarchicalQP& update();

    /**
     * @brief Solve the given optimization problem
     * @return Solver output as joint acceleration command
     */
    virtual const types::JointCommand& solve(const HierarchicalQP& hqp);

    /**
     * @brief Set acceleration regularization term.
     * @param reg This value is added to the diagonal of the Hessian matrix inside the QP to reduce the risk of infeasibility. Default is 1e-8.
     */
    void setAccelerationPenalty(const double reg){acceleration_penalty=reg;}

    /**
     * @brief Set contact Wrench regularization term.
     * @param reg This value is added to the diagonal of the Hessian matrix inside the QP to reduce the risk of infeasibility. Default is 1e-8
     */
    void setContactWrenchPenalty(const double reg){contact_wrench_penalty=reg;}

    /**
     * @brief Set weight for penalizing the difference between consecutive joint accelerations (solver output), i.e., the term
     * \f$ w\|\ddot{\mathbf{q}} - \ddot{\mathbf{q}}_{prev}\|_2^2 \f$ is added to the cost function. This smoothes the solver output over time,
     * which is helpful e.g. on real robots with noisy state estimation. Higher values give smoother, but less reactive motion. Default is 0 (disabled).
     */
    void setAccelerationDeltaPenalty(const double w){acceleration_delta_penalty=w;}

    /**
     * @brief Set weight for penalizing the difference between consecutive contact wrenches (solver output), i.e., the term
     * \f$ w\|\mathbf{f} - \mathbf{f}_{prev}\|_2^2 \f$ is added to the cost function. This smoothes the contact force distribution over time,
     * which is helpful e.g. on real robots with noisy state estimation, where the force distribution may otherwise jump between the
     * contact points. Higher values give smoother, but less reactive force distributions. Default is 0 (disabled).
     */
    void setContactWrenchDeltaPenalty(const double w){wrench_delta_penalty=w;}

    /**
     * @brief Soften the contact surface friction cone constraint using a single slack variable \f$s_i\f$ per contact, i.e., the
     * hard constraint \f$\mathbf{A}_i\mathbf{f}_i \leq \mathbf{0}\f$ is replaced by \f$\mathbf{A}_i\mathbf{f}_i \leq s_i\mathbf{1}, s_i \geq 0\f$
     * and the term \f$ w\sum_i s_i^2 \f$ is added to the cost function. This avoids hard active-set switching (chattering) and
     * infeasibility due to noisy state estimates on a real robot, at the cost of allowing small friction cone violations.
     * Higher values approximate the hard constraint more closely. Only available for surface contacts (dim_contact == 6).
     * Default is 0 (hard constraint, disabled).
     */
    void setFrictionConeSlackPenalty(const double w);
};

} // namespace wbc

#endif
