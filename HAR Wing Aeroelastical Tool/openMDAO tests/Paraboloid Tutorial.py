from __future__ import division, print_function

import openmdao.api as om


class Paraboloid(om.ExplicitComponent):
    """
    Evaluates the equation f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3.
    """

    def setup(self):
        self.add_input('x', val=0.0)
        self.add_input('y', val=0.0)

        self.add_output('f_xy', val=0.0)

        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        """
        f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3

        Minimum at: x = 6.6667; y = -7.3333
        """
        x = inputs['x']
        y = inputs['y']

        outputs['f_xy'] = (x-3.0)**2 + x*y + (y+4.0)**2 - 3.0


if __name__ == "__main__":
    # build the model
    prob = om.Problem()
    indeps = prob.model.add_subsystem('indeps', om.IndepVarComp())
    indeps.add_output('x', 3.0)
    indeps.add_output('y', -4.0)

    prob.model.add_subsystem('parab', Paraboloid())

    # define the component whose output will be constrained
    prob.model.add_subsystem('const', om.ExecComp('g = x + y'))

    prob.model.connect('indeps.x', ['parab.x', 'const.x'])
    prob.model.connect('indeps.y', ['parab.y', 'const.y'])

    # setup the optimization
    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options['optimizer'] = 'COBYLA'

    prob.model.add_design_var('indeps.x', lower=-50, upper=50)
    prob.model.add_design_var('indeps.y', lower=-50, upper=50)
    prob.model.add_objective('parab.f_xy')

    # to add the constraint to the model
    prob.model.add_constraint('const.g', lower=0, upper=10.)
    # prob.model.add_constraint('const.g', equals=0.)

    prob.setup()
    prob.run_driver()

    # minimum value
    print(prob['parab.f_xy'])

    # location of the minimum
    print(prob['indeps.x'])
    print(prob['indeps.y'])