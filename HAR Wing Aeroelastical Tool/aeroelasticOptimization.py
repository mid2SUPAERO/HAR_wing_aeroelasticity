import openmdao.api as om
import numpy as np
from aeroelasticProblemFunc import aeroelasticProblemFunc
import matplotlib.pyplot as plt

class AeroelasticProblem(om.ExplicitComponent):

    def setup(self):
        self.add_input('c', val=0.5)
        self.add_input('t', val=0.005)
        self.add_input('L', val=5.0)
        self.add_input('rootSweepAng', val=0.0)
        self.add_input('rootTwistAng', val=0.0)

        self.add_output('J', val=0.0)

        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        """
        f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3

        Minimum at: x = 6.6667; y = -7.3333
        """
        c = inputs['c']
        t = inputs['t']
        L = inputs['L']
        rootSweepAng = inputs['rootSweepAng']
        rootTwistAng = inputs['rootTwistAng']

        outputs['J'], _, _, _ = aeroelasticProblemFunc(t, c, L, rootSweepAng, rootTwistAng)


if __name__ == "__main__":
    # build the model
    prob = om.Problem()
    indeps = prob.model.add_subsystem('indeps', om.IndepVarComp())
    indeps.add_output('c', 0.5)
    indeps.add_output('t', 0.005)
    indeps.add_output('L', 5.0)
    indeps.add_output('rootSweepAng', 0.0)
    indeps.add_output('rootTwistAng', 0.0)

    prob.model.add_subsystem('aeroelastic', AeroelasticProblem())

    # define the component whose output will be constrained
    prob.model.add_subsystem('const', om.ExecComp('g = c - 10*t'))

    prob.model.connect('indeps.c', ['aeroelastic.c', 'const.c'])
    prob.model.connect('indeps.t', ['aeroelastic.t', 'const.t'])
    prob.model.connect('indeps.L', 'aeroelastic.L')
    prob.model.connect('indeps.rootSweepAng', 'aeroelastic.rootSweepAng')
    prob.model.connect('indeps.rootTwistAng', 'aeroelastic.rootTwistAng')

    # setup the optimization
    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options['optimizer'] = 'COBYLA'

    prob.model.add_design_var('indeps.c', lower=0.1, upper=1)
    prob.model.add_design_var('indeps.t', lower=0.001, upper=0.05)
    prob.model.add_design_var('indeps.L', lower=3, upper=15)
    prob.model.add_design_var('indeps.rootSweepAng', lower=0.0, upper=0.1745)
    prob.model.add_design_var('indeps.rootTwistAng', lower=0.0, upper=0.1745)
    prob.model.add_objective('aeroelastic.J')

    # to add the constraint to the model
    prob.model.add_constraint('const.g', lower=0)
    # prob.model.add_constraint('const.g', equals=0.)

    prob.setup()
    prob.run_driver()

    # minimum value
    print(prob['aeroelastic.J'])

    # location of the minimum
    print(prob['indeps.c'])
    print(prob['indeps.t'])
    print(prob['indeps.L'])
    print(prob['indeps.rootSweepAng'])
    print(prob['indeps.rootTwistAng'])