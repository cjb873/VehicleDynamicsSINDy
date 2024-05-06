from numpy import sin, cos, tan, arctan, pi, absolute


def accel_const(velocity, acceleration, params):

    v_switch = params['v_switch']
    v_min = params['v_min']
    v_max = params['v_max']
    a_max = params['a_max']

    if velocity > v_switch:
        pos_limit = a_max * v_switch / velocity
    else:
        pos_limit = a_max

    if(velocity <= v_min and acceleration <= 0) or \
            (velocity >= v_max and acceleration >= 0):
        acceleration = 0.
    elif acceleration <= -a_max:
        acceleration = -a_max
    elif acceleration >= pos_limit:
        acceleration = pos_limit

    return acceleration


def steering_const(steering_angle, steering_velocity, params):

    s_min = params['s_min']
    s_max = params['s_max']
    sv_min = params['sv_min']
    sv_max = params['sv_max']

    if(steering_angle <= s_min and steering_velocity <= 0) or \
            (steering_angle >= s_max and steering_velocity >= 0):
        steering_velocity = 0.
    elif steering_velocity <= sv_min:
        steering_velocity = sv_min
    elif steering_velocity >= sv_max:
        steering_velocity = sv_max

    return steering_velocity


def e_kinematic_model(x, u_initial, p):
    Iz = p["I"]
    lr = p["lr"]
    lf = p["lf"]
    C_Sr = p["C_Sr"]
    C_Sf = p["C_Sf"]
    g = 9.81
    h = p["h"]
    m = p["m"]
    mu = p["mu"]
    u =[]

    u.append(steering_const(x[3], u_initial[0], p))
    u.append(accel_const(x[4], u_initial[1], p))

    # states
    # x0 - x
    # x1 - y
    # x2 - heading
    # x3 - velocity in x direction
    # x4 - velocity in y direction
    # x5 - angular velocity (change in heading)
    # x6 - steering angle

    # inputs
    # u0 - steering angle velocity
    # u1 - acceleration
    
    Frx = m * u[1]

    f = [x[3] * cos(x[2]) - x[4] * sin(x[2]),
         x[3] * sin(x[2]) + x[4] * cos(x[2]),
         x[5],
         (1 / m) * Frx,
         lr / (lr + lf) * (u[0] * x[3] + x[6] * (Frx * m)),
         1 / (lr + lf) * (u[0] * x[3] + x[6] * (Frx * m)),
         u[0]]

    return f


def kinematic_st_model(x, u_initial, p):
    length = p['lf'] + p['lr']

    # states
    # x0 - x position (global)
    # x1 - y position (global)
    # x2 - steering angle
    # x3 - velocity
    # x4 - heading
    
    # inputs
    # u0 - steering angle velocity
    # u1 - acceleration in x direction

    u = []
    # consider constraints
    u.append(steering_const(x[2], u_initial[0], p))
    u.append(accel_const(x[3], u_initial[1], p))


    f = [x[3] * cos(x[4] % (2*pi)),
         x[3] * sin(x[4] % (2*pi)),
         u[0],
         u[1],
         (x[3] / length) * tan(x[2])]

    return f

def st_model(x, u_initial, p):
    Iz = p['I']
    lr = p['lr']
    lf = p['lf']
    C_Sr = p['C_Sr']
    C_Sf = p['C_Sf']
    g = 9.81
    h = p['h']
    m = p['m']
    mu = p['mu']

    # states
    # x0 - x position
    # x1 - y position
    # x2 - steering angle
    # x3 - velocity
    # x4 - heading
    # x5 - heading derivative
    # x6 - slip angle

    u = []
    # consider constraints
    u.append(steering_const(x[2], u_initial[0], p))
    u.append(accel_const(x[3], u_initial[1], p))

    if absolute(x[3]) < 0.1:
        # Use kinematic model with reference point at center of mass
        # wheelbase
        lwb = lf + lr

        beta = arctan(tan(x[2]) * lr/lwb)
        
        # kinematic model
        f = [x[3] * cos(beta + x[4]),
             x[3] * sin(beta + x[4]),
             u[0],
             u[1],
             x[3] * cos(beta) * tan(x[2]) / lwb]

        # derivative of slip angle and yaw rate
        d_beta = (lr * u[0]) / (lwb * cos(x[2]) ** 2 * (1 + (tan(x[2]) ** 2 * lr / lwb) ** 2))
        dd_psi = 1 / lwb * (u[1] * cos(x[6]) * tan(x[2]) -
                            x[3] * sin(x[6]) * d_beta * tan(x[2]) +
                            x[3] * cos(x[6]) * u[0] / cos(x[2]) ** 2)
        f.append(dd_psi)
        f.append(d_beta)

    else:
        f = [x[3] * cos(x[4] + x[6]),
             x[3] * sin(x[4] + x[6]),
             u[0],
             u[1],
             x[5],
             ((mu * m) / (Iz * (lr + lf))) * (lf * C_Sf * (g * lr - u[1] * h) * x[2] + (lr * C_Sr * (g * lf + u[1] * h) - lf * C_Sf * (g * lr - u[1] * h)) * x[6]\
                                             - (lf ** 2 * C_Sf * (g * lr - u[1] * h) + lr ** 2 * C_Sr * (g * lf + u[1]*h)) * (x[5]/x[3])),
             (mu / (x[3] * (lf+lr))) * (C_Sf *(g * lr - u[1] * h) * x[2] - (C_Sr * (g*lf + u[1] * h) + C_Sf * (g * lr - u[1] * h)) * x[6] \
                                       + (C_Sr * (g * lf + u[1] * h) * lr - C_Sf * (g * lr - u[1]*h) * lf) * (x[5]/x[3])) - x[5]
            ]
    f[4] = f[4] % (2 * pi)
    f[6] = f[6] % (2 * pi)
    return f
