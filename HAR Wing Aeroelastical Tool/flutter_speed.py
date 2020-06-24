import numpy as np
from numpy.linalg import eig


def flutter_speed(Vwind_initial, Vwind_final, tol, ap, strain_eq, altitude, betaeq, keq, throttle, deltaflap, freeDEG, T=np.array([0, 1]), elev=np.array([0, 0])):
    unstable_speed = 1000
    unstable_eig_value = None
    unstable_eig_vec = None

    _, Aaeroelast, _ = ap.linearize(strain_eq, betaeq, keq, throttle, deltaflap, Vwind_initial, freeDEG, T, elev)

    eig_val_initial, eig_vec_initial = eig(Aaeroelast)

    eig_val = 0*np.diag(eig_val_initial)
    eig_vec = 0*eig_vec_initial

    if np.max(eig_val_initial.real) > 0:
        print('Initial speed is already unstable! \n Cant find instability speed in this interval ')
        unstable_speed = Vwind_initial
    else:
        _, Aaeroelast, _= ap.linearize(strain_eq, betaeq, keq, throttle, deltaflap, Vwind_final, freeDEG, T, elev)

        eig_val_final, eig_vec_final = eig(Aaeroelast)

        if np.max(eig_val_final.real) < 0:
            print('Final speed is stable! \n Cant find instability speed in this interval! ')
            unstable_speed = Vwind_final
        else:
            print('Calculating flutter speed')
            diff = Vwind_final - Vwind_initial
            Vnew = (Vwind_final + Vwind_initial) / 2

            while diff > tol:
                _, Aaeroelast, _ = ap.linearize(strain_eq, betaeq, keq, throttle, deltaflap, Vnew, freeDEG, T, elev)

                eig_val, eig_vec = eig(Aaeroelast)
                eig_val = np.diag(eig_val)

                real_max = np.max(eig_val.real)
                ind_max = np.where(np.isclose(eig_val.real, real_max))

                if real_max > 0:
                    Vwind_final = Vnew
                else:
                    Vwind_initial = Vnew

                diff = Vwind_final - Vwind_initial
                Vnew = (Vwind_final + Vwind_initial) / 2

            unstable_speed = Vwind_initial
            unstable_eig_value = eig_val[ind_max[0][0], ind_max[1][0]]
            unstable_eig_vec = eig_vec[:, ind_max[1][0]]

    return unstable_speed#, unstable_eig_value, unstable_eig_vec
