import copy
import random
import time

import numpy as np
from matplotlib import pyplot as plt

from run_nbody_sim import G, AU, mSol, bf, rSol, mEarth, rEarth, dt, \
    elapsedTimeBH
from body import body


class systemBF:
    def __init__(self):
        #self.ax = plt.subplots()                            # not convinced this is necessary
        self.bodies = []                                    # list of bodies


    def addBod(self, body):
        self.bodies.append(body)


    def removeBod(self, body):
        self.bodies.remove(body)


    def EscapeVelocity(self, body):
        mtot = self.bodies[0].mass
        vEscape = np.sqrt( 2*G*mtot/d)


    def CenterOfMass(self, bodies):
        mcomx = 0
        mcomy = 0
        mtot = 0
        for i in bodies:
            mcomx += i.mass * i.xs[-1]
            mcomy += i.mass * i.ys[-1]
            mtot += i.mass
        com = np.sqrt( (mcomx/mtot)**2 + (mcomy*mtot)**2 )
        return com


    def collisions(self, body1, body2):
        d = np.sqrt( (body1.xs[-1] - body2.xs[-1])**2 + (body1.ys[-1] - body2.ys[-1])**2 )
        # Are we colliding?
        if d < body1.radius + body2.radius:
            # Kill one body, add its mass to the other.
            body1.mass += body2.mass
            body1.radius = rocheLimit(body1.mass)
            self.removeBod(body2)


    def plotRV(self, nTimesteps, dt, starIndex):
        nTimesteps = nTimesteps/1000
        plt.ion()
        plotVerticalRange = AU
        star = self.bodies[starIndex]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.axis([0, nTimesteps, -plotVerticalRange, plotVerticalRange])

        for t in range(nTimesteps):
            self.step_BruteForce(dt*6)
            ax.plot(t, 0, '-k')

            # Find CoM so that it can be subtracted to keep the star in place
            mxcomx = 0
            mtot = 0
            for i in self.bodies:
                mxcomx += i.xs[-1]*i.mass
                mtot += i.mass
            comx = mxcomx/mtot

            #ax.plot(t, self.bodies[0].xs[-1], '.y')
            #ax.plot(t, self.bodies[0].ys[-1], '.r')
            #ax.axvline(x=t)
            #ax.plot(t, self.bodies[0].xs[-1] - comx, '.y')
            ax.plot(t, self.bodies[0].velocity[0], '.y')
            plt.show(block=False)
            plt.pause(0.0001)


    def step_BruteForce(self, dt):
        # To be sure that we're not taking updated positions, copy all the bodies' current positions and analyze from that.
        oldBods = copy.deepcopy(self.bodies)
        for i in self.bodies:
            #print "\nCurrently Brute-Force stepping ", i.name
            # Calculate the gravitaional acceleration due to all other bodies on i
            for ob in oldBods:
                i.aGrav(ob)
            # Now update the positions and velocities.
            i.positionAndVelocityUpdater(dt)


            # If a body is too far out of the FOV, just get rid of it.
            # Removing this for the time being just so no bodies get lost -> more accurate timing
            """
            if abs(np.sqrt( (i.xs[-1])**2 + (i.ys[-1])**2)) > 1.5*axlims:
                self.removeBod(i)
                print "removed", i.name, "for being too far away"
                """


    def run_BruteForce(self, nTimesteps, nSmallBods):
        # Start the timer:
        startBF = time.time()
        mStar = mSol
        # exampleBody =  body(name, mass, radius, xy, velocity, color)
        bf.addBod(body('Star', mStar, rSol, [0,0], [0,-0], 'y'))
        #bf.addBod(body('Mercury', 0.055*mEarth, 0.3829*rEarth, [-0.4481*AU, 0], [0, -55410], 'blue'))
        bf.addBod(body('Venus', 0.815*mEarth, 0.949*rEarth, [0.721*AU, 0], [0, 34910], 'orange'))
        bf.addBod(body('Earth', mEarth, rEarth, [0, AU], [-29838, -0], 'g'))
        #bf.addBod(body('Mars', 0.10745*mEarth, 0.531*rEarth, [0, -1.52*AU], [240740, 0], 'red'))
        #bf.addBod(body('Jupiter', 317*mEarth, 11*rEarth, [5.2*AU, 0], [0, 13048], 'magenta'))
        #bf.addBod(body('FastJupiter', 317*mEarth, 11*rEarth, [5.2*AU, 0], [0, 40000], 'magenta'))
        #bf.addBod(body('Saturn', 95.16*mEarth, 9.14*rEarth, [0,10.06*AU], [-10180, 0], 'orange'))
        #bf.addBod(body('Uranus', 14.53*mEarth, 3.976*rEarth, [-19.91*AU, 0], [0, -7058], 'blue'))
        #bf.addBod(body('Neptune', 17.148*mEarth, 3.86*rEarth, [0, -29.95*AU], [5413, 0], 'cyan'))

        # V2434 Ori
        #bf.addBod(body('V2434a', 3.5*mSol, 3.5*rSol, [-220*AU, 0], [0, -np.sqrt((-G*6*mSol)/(220*AU))], 'blue'))
        #bf.addBod(body('V2434b', 3.*mSol, 3*rSol, [220*AU, 0], [0, np.sqrt((-G*6*mSol)/(220*AU))], 'cyan'))
        #"""
        # Inner ring of small random bodies:
        if nSmallBods > 0:
            for i in range(0, nSmallBods):
                print "Adding small bod!", i
                r = random.uniform(4.9*AU, 5.1*AU)
                x = random.uniform(0,r)
                y = np.sqrt(r**2 - x**2)

                # velocity direction is the sign of the quadrant times the sign of the coordinate (x or y)

                # THIS IS INVERTED RELATIVE TO THE BF VERSION
                v = np.sqrt(r*AU/(-G*mStar))

                """
                s = (x*y)/abs(x*y)
                # Add some random directions
                variance = random.uniform(-1000, 1000)
                vx = variance * v * s * (x/abs(x))
                vy = variance * v * s * (y/abs(y))
                """

                variance = random.uniform(-1000, 1000)
                vx = -v * (y/abs(y)) * np.sin(x/y) + variance
                vy = v * (x/abs(x)) * np.cos(x/y) + variance

                bf.addBod(body('planetesimal_'+str(i), 0.05*mEarth, 0.02*rEarth, [x, y], [vx, vy], 'y'))



        # Make a file that has the filename of each body, initiate body files.
        fnames = []
        f = open('filenames.log', 'w')
        f.write("Name\n")
        f.close()

        print self.bodies

        # Initiate each body's file
        for i in range(0, len(self.bodies)):
            b = self.bodies[i]
            fnames.append(b.fname)
            f = open(b.fname, 'w')
            f.write("XS YS\n")
            f.close()

            f = open('filenames.log', 'a')
            f.write(b.fname + '\n')
            f.close()

        # Do the actual run:
        for t in range(nTimesteps):
            # Take your step
            self.step_BruteForce(dt)

            print "\n\nBF Timestep:", t, "finished"
            stepTime = time.time()
            elapsedTimeBH.append([t, stepTime - startBF])


        print "Final Duration with", nSmallBods, "small bodies over", nTimesteps, "steps (seconds):", time.time() - startBF
        #plotTrajectories()
        return [nSmallBods, time.time() - startBF]













