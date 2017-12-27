import numpy as np

from run_nbody_sim import G


class body:
    def __init__(self, name, mass, radius, xy, velocity, color):
        self.name = name
        self.mass = mass                                    # kg (str)
        self.radius = radius                                # m (float)
        self.xs = [xy[0]]                                   # m (floats)
        self.ys = [xy[1]]                                   # m (floats)
        self.velocity = velocity                            # [x,y] [m/s, m/s] (floats)
        self.acceleration = [0,0]                           # [x,y] [m s-2, m s-2] (floats)
        self.color = color                                  # plot color (str)
        self.hostNode = None                                # Useful for gridding stuff later
        self.fname = name + '.txt'


    def isEqual(self, otherBody):
        if self.name == otherBody.name:
            return True
        return False


    def update_gravitational_accel(self, otherBody):
        # Calculate the gravitational acceleration on a body due to a different body
        # Zero this out elsewhere and add on.
        #self.acceleration = [0,0]
        if self.name != otherBody.name:
            # Update gravitational acceleration (a = f/m) -> self.m cancels
            # take x or y projections with x/y hats
            x = (self.xs[-1] - otherBody.xs[-1])
            y = (self.ys[-1] - otherBody.ys[-1])
            d = (np.sqrt(x**2 + y**2))
            xhat = x/d
            yhat = y/d
            self.acceleration[0] += xhat * G * otherBody.mass * d**-2
            self.acceleration[1] += yhat * G * otherBody.mass * d**-2


    def update_position_and_velocity(self, dt):

        self.velocity[0] += self.acceleration[0] * dt           # x component
        self.velocity[1] += self.acceleration[1] * dt           # y component

        self.xs.append(self.xs[-1] + self.velocity[0]*dt + 0.5*self.acceleration[0]* (dt**2))
        self.ys.append(self.ys[-1] + self.velocity[1]*dt + 0.5*self.acceleration[1]* (dt**2))


        # Now write this position change to the body's file.
        g = open(self.fname, 'a')
        outstr = str(self.xs[-1]) + " " + str(self.ys[-1]) + '\n'
        g.write(outstr)
        g.close()

