#ifndef WBC_CORE_CONTACTS_FRICTION_SURFACE_CONSTRAINT_HPP
#define WBC_CORE_CONTACTS_FRICTION_SURFACE_CONSTRAINT_HPP

#include "../core/Constraint.hpp"

namespace wbc {

/**
 * @brief Polyhedral wrench friction (coulomb) cone constraint according to
 * https://scaron.info/robotics/wrench-friction-cones.html (Wrench friction cone for surface contacts)
 *
 * If use_slack is set, the cone is softened using a single slack variable s per contact: A*f <= s, s >= 0.
 * The slack variables are assumed to be appended to the qp variables, i.e., the variable order becomes
 * (qdd,f_ext,s). The bounds s >= 0 and the penalty on s have to be added by the scene.
 */
class ContactsFrictionSurfaceConstraint : public Constraint{
public:
    /** @brief Default constructor */
    explicit ContactsFrictionSurfaceConstraint(bool _reduced = false, bool _use_slack = false) :
        Constraint(Constraint::inequality), reduced(_reduced), use_slack(_use_slack) { }

    virtual ~ContactsFrictionSurfaceConstraint() = default;

    virtual void update(RobotModelPtr robot_model) override;

    /** @brief Enable/disable softening of the friction cone using one slack variable per contact */
    void setUseSlack(bool enable){use_slack = enable;}

    /** @brief Return whether the friction cone is softened using slack variables */
    bool useSlack() const {return use_slack;}

private:
    bool reduced; // if torques are removed from the qp formulation or not
    bool use_slack; // if the friction cone is softened using one slack variable per contact
};

}

#endif // WBC_CORE_CONTACTS_FRICTION_SURFACE_CONSTRAINT_HPP
