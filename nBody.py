# If you just have a stupid sense of humor, you'll never run out of things to laugh at
# "The larger our ignorance, the stronger the magnetic field." - Woltier
"""
Exoplanets Final Project: N-Body Simulation
@author: Jonas Powell
October 2017
"""


## PACKAGES ###
import os
import numpy as np
import matplotlib.pyplot as plt
import copy
import random
from random import randint
import imageio
import pandas as pd
import time
import matplotlib.animation as animation
import datetime
import astropy.constants as c
import astropy.units as u
#from twilio.rest import Client


### CONSTANTS ###
G = -c.G.value                                              # m3 kg-1 s-2
AU = c.au.value                                             # m
mSol = c.M_sun                                              # kg
rSol = c.R_sun                                              # m
mEarth = c.M_earth                                          # kg
rEarth = c.R_earth                                          # m


# This is the only place dt is defined (not passed as an argument anywhere). JSYK

# For stellar scales (v2434)
#dt = 30000
#axlims = 220 * AU

dt = 1800                                                   # duration of one timestep (s) (1800s = half hour)
axlims = 2 * AU                                             # AU
theta = 0.0                                                 # for BH Approx
steps = []                                                  # the output file

# Automatically make good new file names.
now = datetime.datetime.now()
months = ['jan', 'feb', 'march', 'april', 'may', 'june', 'july', 'aug', 'sep', 'oct', 'nov', 'dec']
today = months[now.month - 1] + str(now.day)
logfileName = today + 'run.log'
gifOutName = today + 'run.gif'
trajOutName = 'z' + today + 'run.png'

gifNames = []
elapsedTimeBF = []
durationBF = []
elapsedTimeBH = []
durationBH = []







#######################
### INDIVIDUAL BODY ###
#######################

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


    def aGrav(self, otherBody):
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


    def positionAndVelocityUpdater(self, dt):

        self.velocity[0] += self.acceleration[0] * dt           # x component
        self.velocity[1] += self.acceleration[1] * dt           # y component

        self.xs.append(self.xs[-1] + self.velocity[0]*dt + 0.5*self.acceleration[0]* (dt**2))
        self.ys.append(self.ys[-1] + self.velocity[1]*dt + 0.5*self.acceleration[1]* (dt**2))


        # Now write this position change to the body's file.
        g = open(self.fname, 'a')
        outstr = str(self.xs[-1]) + " " + str(self.ys[-1]) + '\n'
        g.write(outstr)
        g.close()






#############################
### SYSTEM (BRUTE FORCE) ####
#############################


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











###########################
### BARNES-HUT GRIDDING ###
###########################




class node:
    def __init__(self, name, xlims, ylims, parent):
        self.name = name
        self.bodies = []
        self.bodiesContained = []
        self.xlims = xlims
        self.ylims = ylims
        self.com = [0,0]
        self.totalMass = 0
        self.ne = None
        self.nw = None
        self.sw = None
        self.se = None
        self.parent = parent
        self.children = []
        self.lineage = None

        if self.parent == None:
            self.level = 0
        else:
            self.level = self.parent.level + 1




class tree:
    def __init__(self, initialNode):
        self.root = initialNode

    def addBod(self, n, body):
        # Give a body and a host node for that body
        self.sortBody(n, body)
        #self.traceback(n)



    def remove(self, body):
        # Remove a body from the structure by emptying its host node.
        # This should work since, after the sorting, no node should be holding more than one body.
        bods = self.searchTree(self.root, 'bodies')
        for i in range(len(bods)):
            # If there's a body here
            if len(bods[i]) > 0:
                # If these bodies are the same one (i.e. it's the one we're looking for)
                if bods[i][0].isEqual(body):
                    bods[i][0].hostNode.bodies = []



    def traceback(self, n):
        path = ''
        while n.name != 'root':
            path = (n.name + '-' + path)
            n = n.parent
        path = 'root' + '-' + path
        n.lineage = path



    def splitNode(self, n):
        # Divide the node into four subs
        nodeMidX = np.mean([n.xlims[0], n.xlims[1]])
        nodeMidY = np.mean([n.ylims[0], n.ylims[1]])

        n.ne = node('ne', [nodeMidX, n.xlims[1]], [nodeMidY, n.ylims[1]], n)
        n.nw = node('nw', [n.xlims[0], nodeMidX], [nodeMidY, n.ylims[1]], n)
        n.sw = node('sw', [n.xlims[0], nodeMidX], [n.ylims[0], nodeMidY], n)
        n.se = node('se', [nodeMidX, n.xlims[1]], [n.ylims[0], nodeMidY], n)

        # Add each node's lineage
        for i in [n.ne, n.nw, n.sw, n.se]:
            self.traceback(i)



    def sortBody(self, n, body):
        # Sort a body into the tree
        # Add body's mass to the total mass of each node it passes through.
        # If the body is trying to go somewhere outside the node, just drop it.
        """ Mechanics

            Given a body and a node:
                add the body to the node
                if there are other bodies in the same node:
                    split the node into four subnodes
                    determine which quadrant each body should be in
                    add each body to the appropriate subnode

                    if there are other bodies...
                    """


        # If the body is trying to go outside the node's limits, just pass it (delete it)
        if body.xs[-1] > n.xlims[1] or  body.xs[-1] < n.xlims[0] or body.ys[-1] > n.ylims[1] or  body.ys[-1] < n.ylims[0]:
            # Just return. If we've made it this far, then the body hasn't stuck anywhere above, so by not placing it in this one, it just disappears
            return


        # Want to add mass to every body as we go.
        n.totalMass += body.mass
        n.bodiesContained.append(body)

        # 1. No other bodies, no other subnodes available.
        if len(n.bodies) == 0 and n.ne == None:
            # Add a body to the current node so we can handle them all at once
            n.bodies.append(body)
            body.hostNode = n
            #print body.name, "into", n.name, ", finished."
            return


        # 2. There are other subnodes or other bodies or both
        else:
            # 2.1: There are other bodies
            if n.ne != None:
                #print "No other bodies, but other subs to hit"
                # Sort it out.
                nodeMidX = np.mean([n.xlims[0], n.xlims[1]])
                nodeMidY = np.mean([n.ylims[0], n.ylims[1]])

                if nodeMidX <= body.xs[-1]:                     # if body is on the east side:
                    if nodeMidY <= body.ys[-1]:                 # if body is on the north side:
                        #print "NE", n.ne.level, "for ", body.name
                        self.sortBody(n.ne, body)
                        #n.bodies.remove(b)
                    elif nodeMidY >= body.ys[-1]:               # it's on the south
                        #print "Nw", n.nw.level, " for ", body.name
                        self.sortBody(n.se, body)
                        #n.bodies.remove(b)
                else:                                           # body is on the west side
                    if nodeMidY <= body.ys[-1]:                 # if body is on the north side:
                        #print "SW", n.sw.level, "for ", body.name
                        self.sortBody(n.nw, body)
                        #n.bodies.remove(b)

                    elif nodeMidY >= body.ys[-1]:               # body is on the south side
                        #print "SE", n.se.level, "for ", body.name
                        self.sortBody(n.sw, body)
                        #n.bodies.remove(b)


            # 2.2. There are other bodies but no subnodes.
            else:
                # Add a body to the current node so we can handle them all at once
                n.bodies.append(body)

                #print "Too many bodies!", len(n.bodies), "in", n.level, n.name
                #print "Bodies: ", [i.name for i in n.bodies]


                # Create subnodes
                self.splitNode(n)

                # Loop over the bodies in this quadrant
                for b in n.bodies:
                    nodeMidX = np.mean([n.xlims[0], n.xlims[1]])
                    nodeMidY = np.mean([n.ylims[0], n.ylims[1]])

                    # BODY IS ON EAST
                    if nodeMidX <= b.xs[-1]:                    # if body is on the east side:
                        if nodeMidY <= b.ys[-1]:                # if body is on the north side:
                            #print "NE Xlims", n.ne.xlims
                            #print "NE", n.ne.level, "for ", b.name
                            self.sortBody(n.ne, b)
                            #n.bodies.remove(b)
                        elif nodeMidY >= b.ys[-1]:              # it's on the south
                            #print "SE Xlims", n.se.xlims
                            #print "SE", n.se.level, " for ", b.name
                            self.sortBody(n.se, b)
                            #n.bodies.remove(b)
                    # BODY IS ON WEST
                    else:                                       # body is on the west side
                        if nodeMidY <= b.ys[-1]:                # if body is on the north side:
                            #print "NW Xlims", n.nw.xlims
                            #print "NW", n.nw.level, "for ", b.name
                            self.sortBody(n.nw, b)
                            #n.bodies.remove(b)

                        elif nodeMidY >= b.ys[-1]:             # body is on the south side
                            #print "SW Xlims", n.sw.xlims
                            #print "SW", n.sw.level, "for ", b.name
                            self.sortBody(n.sw, b)
                            #n.bodies.remove(b)


                # Remove those bodies from this subnode
                for i in range(len(n.bodies)):
                    b = n.bodies[0]
                    n.bodies = n.bodies[1:]
                    #print "removing ", b.name, "from level", n.level
                    #print "  bodies remaining:", len(n.bodies)



    def searchTree(self, rootNode, param):
        # Search the tree for a param using Breadth First Search
        # BFS inspired by pseudo code on Wikipedia (Tree Traversal)
        # Make sure to call the param as a string.

        # Pull param from each body held in nInit and return as a list
        paramList = []
        q = [rootNode]
        while len(q) > 0:
            # Grab the node we're looking for from queue list
            n = q[0]
            # Remove that node from the list
            q = q[1:]
            # Grab what we're looking for from the node (somehow sort it? Dict maybe better?)
            paramList.append(getattr(n, param))

            # Add more nodes to our list
            if n.ne != None:
                q.append(n.ne)
            if n.nw != None:
                q.append(n.nw)
            if n.sw != None:
                q.append(n.sw)
            if n.se != None:
                q.append(n.se)

        return paramList



    def calculateCenterOfMass(self, n):
        Mcomx = 0
        Mcomy = 0
        for i in n.bodiesContained:
            Mcomx += i.xs[-1]*i.mass
            Mcomy += i.ys[-1]*i.mass

        n.com = [Mcomx/(n.totalMass), Mcomy/(n.totalMass)]
        #print "Total Mass in this Node (kg)", n.totalMass
        #print "Center of mass (AU)", [n.com[0]/AU, n.com[1]/AU]



    def visualizeTree(self, t):
        # GRID PLOTTING
        # Get the boundaries:
        xlims = self.searchTree(self.root, 'xlims')
        ylims = self.searchTree(self.root, 'ylims')

        # Set up the figure:
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect = 'equal')
        ax.axis([xlims[0][0], xlims[0][1], ylims[0][0], ylims[0][1]])

        for i in range(0, len(xlims)-1):
            ax.plot((xlims[i][0], xlims[i][1]), (ylims[i][0], ylims[i][0]), '-k')
            ax.plot((xlims[i][0], xlims[i][1]), (ylims[i][1], ylims[i][1]), '-k')

        for i in range(0, len(ylims)-1):
            ax.plot((xlims[i][0], xlims[i][0]), (ylims[i][0], ylims[i][1]), '-k')
            ax.plot((xlims[i][1], xlims[i][1]), (ylims[i][0], ylims[i][1]), '-k')

        # BODY PLOTTING
        bods = self.searchTree(self.root, 'bodies')

        for i in range(0, len(bods)):
            if len(bods[i]) > 0:
                b = bods[i][0]
                circle = plt.Circle((b.xs[-1], b.ys[-1]), radius=axlims/40, color=b.color, alpha=0.7)
                ax.add_artist(circle)


        # Calculate the CoM of each body
        self.calculateCenterOfMass(self.root)
        #com = plt.Circle((self.root.nw.com[0], self.root.nw.com[1]), radius=1*AU, color="black", alpha=0.7)
        #print "VISUALIZE COM: ", self.root.com
        #ax.add_artist(com)
        #plt.legend()
        plt.title('Barnes-Hut Orbital Paths')
        plt.xlabel('Distance (meters)')
        plt.ylabel('Distance (meters)')
        outstr = "graph_moment" + str(t) + ".jpeg"
        gifNames.append(outstr)
        print len(gifNames)
        plt.savefig(outstr)
        plt.close('all')
        #plt.show(block=False)



    def calculateGravBH(self, n, b):
        """ MECHANICS:
            For a given node and body:
            Calculate s/d
            If s/d < theta, calculate that node's CoM and gravitational effect.
            else, try each of those inner nodes
            """

        # If this node is a leaf (i.e. no subs, just a nice little body), just do the grav
        if len(n.bodies) > 0:
            b.aGrav(n.bodies[0])
            return

        # There are no bodies, so it's either an empty leaf or an empty internal branch
        else:
            # Mass is a better way of evaluating this than len(bodies) since an internal node won't have bodies but will have totalMass
            if n.totalMass > 0:
                #print "mass > 0"

                # Don't compare the body to its own node
                if b in n.bodiesContained:
                    self.calculateGravBH(n.ne, b)
                    self.calculateGravBH(n.nw, b)
                    self.calculateGravBH(n.sw, b)
                    self.calculateGravBH(n.se, b)

                else:
                    #print " Let's get S/D"
                    self.calculateCenterOfMass(n)
                    s = abs(n.xlims[1] - n.xlims[0])/2
                    d = np.sqrt((b.xs[-1] - n.com[0])**2 + (b.ys[-1] - n.com[1])**2)

                    if s/d > theta:
                        #print "     too big"
                        #print "Approx. too big; going back again."
                        self.calculateGravBH(n.ne, b)
                        self.calculateGravBH(n.nw, b)
                        self.calculateGravBH(n.sw, b)
                        self.calculateGravBH(n.se, b)


                    # s/d < theta, approximation works!
                    else:
                        # Add this acceleration contribution
                        b.aGrav(body('BHapproximated', n.totalMass, 1, n.com, [0,0], 'black'))
                        return



    def RVplotBH(self):
        #plt.ion()
        plotVerticalRange = AU

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.axis([0, nTimesteps, -plotVerticalRange, plotVerticalRange])

        # A little bit junky. Recall that searchTree returns a backwards list, so the Sun is at the end.

        emptybods = self.searchTree(self.root, 'bodies')
        bods = []
        for i in len(emptybods):
            if len(emptybods[i]) > 0:
                bods.append(i)

        star = bods[-1]
        print star.name


        xs = []
        ts = []

        for t in range(len(star.xs)):

            # Find CoM so that it can be subtracted to keep the star in place
            comx = star.hostNode.com[0]
            xs.append(star.xs[t])
            ts.append(t*dt)

        #ax.plot(t, self.bodies[0].xs[-1], '.y')
        #ax.plot(t, self.bodies[0].ys[-1], '.r')
        #ax.axvline(x=t)
        #ax.plot(t, self.bodies[0].xs[-1] - comx, '.y')
        ax.plot(ts, xs, '.y')
        plt.show(block=False)
        #plt.pause(0.0001)



    def step_BH(self, nOld):
        """ PURPOSE:
            Take on time step, moving all bodies in the system appropriately using the Barnes-Hut approximation.
            Begin by gathering all the bodies into a list.
            For each body:
                - From some node n (starting at rootNode), see if the BH approx (s/d < theta) works.
                    - If it does, calculate the gravity exerted by a fabricated body of mass n.totalMass and located at n.com.
                    - If it doesn't, do this whole process to each of its subnodes.

            """



        # Get the list of bodies
        # This might be pretty expensive
        bods = self.searchTree(n, 'bodies')

        positions = np.array([])
        for i in range(0, len(bods)):
            # We only care about this if there are bodies in the n.bodies list.
            if len(bods[i]) > 0:
                b = bods[i][0]
                #print "\nCurrently Barnes-Hut stepping ", b.name

                # Position update (and log the data)
                b.acceleration = [0,0]
                self.calculateGravBH(nOld, b)
                b.positionAndVelocityUpdater(dt)

                # Sort this body into the new graph
                b.parent = self.root
                self.addBod(self.root, b)



    def run_BarnesHut(self, nTimesteps, nSmallBods):
        # Start the timer for cost analysis:
        startBH = time.time()
        # exampleBody =  body(name, mass, radius, xy, velocity, color)
        mStar = mSol
        self.addBod(n, body('Star', mStar, rSol, [1, 1], [0,-0], 'y'))
        self.addBod(n, body('Mercury', 0.055*mEarth, 0.3829*rEarth, [-0.4481*AU, 0], [0, -55410], 'blue'))
        #self.addBod(n, body('Venus', 0.815*mEarth, 0.949*rEarth, [0.721*AU, 0], [0, 34910], 'orange'))
        #self.addBod(n, body('Earth', mEarth, rEarth, [0, AU], [-29838, -0], 'g'))
        #self.addBod(n, body('Mars', 0.10745*mEarth, 0.531*rEarth, [0, -1.52*AU], [24074, 0], 'red'))
        #self.addBod(n, body('Jupiter', 317*mEarth, 11*rEarth, [5.2*AU, 0], [0, 13048], 'magenta'))
        #self.addBod(n, body('Saturn', 95.16*mEarth, 9.14*rEarth, [0,10.06*AU], [-10180, 0], 'orange'))
        #self.addBod(n, body('Uranus', 14.53*mEarth, 3.976*rEarth, [-19.91*AU, 0], [0, -7058], 'blue'))
        self.addBod(n, body('Neptune', 17.148*mEarth, 3.86*rEarth, [0, -29.95*AU], [5413, 0], 'cyan'))
        #self.addBod(n, body('HotJupiter', 317*mEarth, 11*rEarth, [0.1*AU, 0], [0, 17315], 'red'))

        # V2434 Ori
        #self.addBod(n, body('V2434a', 3.5*mSol, 3.5*rSol, [-220*AU, 0], [0, -np.sqrt((-G*6*mSol)/(220*AU))], 'blue'))
        #self.addBod(n, body('V2434b', 3.*mSol, 3*rSol, [220*AU, 0], [0, np.sqrt((-G*6*mSol)/(220*AU))], 'cyan'))
        #"""
        # Inner ring of small random bodies:
        if nSmallBods > 0:
            for i in range(0, nSmallBods):
                r = random.uniform(4.9*AU, 5.1*AU)
                x = random.uniform(-r,r)
                signy = (-1)**(random.randint(0,1))
                y = signy * np.sqrt(r**2 - x**2)

                # velocity direction is the sign of the quadrant times the sign of the coordinate (x or y)
                # Should really just do this with trig
                v = np.sqrt((-G*mStar)/(r*AU))

                """
                s = (x*y)/abs(x*y)
                # Add some random direction variability
                variance = random.uniform(-1000, 1000)
                vx = -v * s * (x/abs(x)) * np.sin(x/y) + variance
                vy = v * s * (y/abs(y)) * np.cos(x/y) + variance

                (y/abs(y)) * (x*y)/abs(x*y) = (x * y2)/abs(x * y2) = 1 * (x/abs(x))
                """
                variance = random.uniform(-1000, 1000)
                vx = -v * (y/abs(y)) * np.sin(x/y) + variance
                vy = v * (x/abs(x)) * np.cos(x/y) + variance

                bh.addBod(n, body('planetesimal_'+str(i), 0.05*mEarth, 0.02*rEarth, [x, y], [vx, vy], 'k'))



        # Grab them bods
        bods = self.searchTree(n, 'bodies')

        # Initiate a file that has the filename of each body so that later we know what to be looking for.
        fnames = []
        f = open('filenames.log', 'w')
        f.write("Name\n")
        f.close()

        # Initiate each body's file and add its name to the filenames log
        for i in range(0, len(bods)):
            if len(bods[i]) > 0:
                b = bods[i][0]
                fnames.append(b.fname)
                f = open(b.fname, 'w')
                f.write("XS YS\n")
                f.close()

                f = open('filenames.log', 'a')
                f.write(b.fname + '\n')
                f.close()


        # Do the actual stepping.
        for t in range(nTimesteps):

            # Clear out subnodes under root for each iteration
            nOld = copy.deepcopy(self.root)
            self.root = node('root', [-axlims, axlims], [-axlims, axlims], None)


            self.step_BH(nOld)
            print "\n\nBH Timestep:", t, "finished"
            stepTime = time.time()
            elapsedTimeBH.append([t, stepTime - startBH])

            # Have a gander at the situation if you want:
            #self.visualizeTree(t)

        print "Final Duration with", nSmallBods, "small bodies over", nTimesteps, "steps (seconds):", time.time() - startBH
        plotTrajectories()
        return [nSmallBods, time.time() - startBH]








#################
### Run Stuff ###
#################

n = node('root', [-axlims, axlims], [-axlims, axlims], None)
bh = tree(n)
#bh.visualizeTree()
#bh.run_BarnesHut(30000, 0)

bf = systemBF()
#bf.run_BruteForce(100, 0)








################
### ANALYSIS ###
################



# Draw Trajectories from Data (xs,ys)
def plotTrajectories(outAsFile=-1):
    # Make a trajectory (smeared) plot from a some coordinate outfiles. Show or savefig().
    outname = trajOutName
    filenames = pd.read_table('filenames.log').Name

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect = 'equal')
    ax.axis([-axlims, axlims, -axlims, axlims])
    bg = plt.Circle((0, 0), radius=np.sqrt(2)*axlims, color='black', alpha=0.9)
    ax.add_artist(bg)

    for f in filenames:
        print f
        p = pd.read_table(f, delim_whitespace=True)
        plt.plot(p.XS, p.YS)

    #plt.legend()
    plt.title('Orbital Paths')
    plt.xlabel('Distance (meters)')
    plt.ylabel('Distance (meters)')

    if outAsFile == 1:
        plt.savefig(outname)
        print "Trajectories plotted in", outname
    elif outAsFile == -1:
        plt.show(block=False)
    else:
        print "Do you want the plot saved as a figure or shown on the screen? (1 = save as fig, -1 = show on screen)"





# Draw Trajectories from Data (xs,ys)
def plotLines(saveFigs = -1):
    # To mimic this: https://www.facebook.com/permalink.php?story_fbid=536038580109918&id=100011113415766&notif_id=1513793777241322&notif_t=feedback_reaction_generic_tagged

    # Obviously, do a run before doing this.

    # saveFigs: If you want to save every step as an image (to turn into a gif), use saveFigs=1. Otherwise (saveFigs=-1, default), just show the image.

    axlims = 1.5 * AU
    outname = 'orbitLines'

    # Figure setup
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect = 'equal')
    ax.axis([-axlims, axlims, -axlims, axlims])
    bg = plt.Circle((0, 0), radius=np.sqrt(2)*axlims, color='black', alpha=0.9)
    ax.add_artist(bg)


    # Janky (non-general) for now; fix later.
    sun = pd.read_table('Star.txt', delim_whitespace=True)
    earth = pd.read_table('Earth.txt', delim_whitespace=True)
    jupiter = pd.read_table('Jupiter.txt', delim_whitespace=True)


    plt.plot(sun.XS[0], sun.YS[0], '.r')

    i = 0
    # It should be redundant to check both (each should have the same number of timesteps), but what the hell.
    while i < len(earth.XS) and i < len(venus.XS):
        plt.plot([earth.XS[i], jupiter.XS[i]], [earth.YS[i], jupiter.YS[i]], color='blue', alpha=0.3, linestyle='-')
        i += 10

        if saveFigs == 1:
            plt.savefig(outname + str(i) + '.jpeg')

    if saveFigs == -1:
        plt.show(block=False)



# Make Timestepped Pics from Data (xs,ys)
def makePics(fnames, gifYN):
    # Figure out how many timesteps we'll need
    # This is junky right now bc it assumes that every body has the same number of timesteps
    # This is fucking slow. plotTrajectories is way faster.
    # The only good time to use this is if run_BarnesHut didn't actually output pictures but you really want a gif.
    print fnames[0]
    ts = len(pd.read_table(fnames[0], delim_whitespace=True).XS)
    picnames = []

    # Plot only once every 100 timesteps
    for t in range(0, ts):
        print "Making pics for step", t
        # Set up the plot
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect = 'equal')
        ax.axis([-axlims, axlims, -axlims, axlims])
        bg = plt.Circle((0, 0), radius=np.sqrt(2)*axlims, color='blue', alpha=0.8)
        ax.add_artist(bg)

        # Grab the position of each body
        for f in fnames:
            f = pd.read_table(f, delim_whitespace=True)
            # Only do this if the data is there
            if t < len(f.XS):
                circle = plt.Circle((f.XS[t], f.YS[t]), radius=0.01*axlims, color='orange', alpha=1)
                ax.add_artist(circle)

        #plt.legend()
        plt.title('Orbital Paths')
        plt.xlabel('Distance (meters)')
        plt.ylabel('Distance (meters)')
        outname = 'out' + str(t) + '.png'
        plt.savefig(outname)
        #plt.show(block=False)

        # Add this picture's name to a list to feed to the gif maker
        picnames.append(outname)


    # Want a gif?
    if gifYN == 'y':
        print "Making GIF: ", gifOutName
        makeGif(picnames, gifOutName)
        print "Finished GIF:", gifOutName

#makePics(fnames, 'y')




### GIF MAKER ###
def makeGif(nSteps):
    # This assumes that the filename log is called 'filenames.log'
    images = []

    for n in range(0,nSteps):
        filename = "graph_moment" + str(n) + ".jpeg"
        print "Adding: ", filename, "to gif", gifOutName
        images.append(imageio.imread(filename))
    imageio.mimsave(gifOutName, images)
    print "Gif created:", gifOutName




def plotRV(filename):
    # Why the fuck isn't this working
    star = pd.read_table(filename, delim_whitespace=True)
    plotVerticalRange = 1.5*rSol

    xs = np.array(star.XS)
    ys = np.array(star.YS)
    ts = np.arange(0, len(xs))
    print xs[2]
    print "Ts: ", ts[:10]
    print "Xs: ", xs[:10]
    print "Length:", len(xs)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.axis([-plotVerticalRange, plotVerticalRange, -plotVerticalRange, plotVerticalRange])


    plt.plot(xs, ys)
    """
    for t in range(len(star.XS)):
        #print star.XS[t]
        circle = plt.Circle((star.XS[t], star.YS[t]), radius=100, color='red', alpha=0.7)
        ax.add_artist(circle)
        #"""


    plt.show(block=False)

#plotRV('Star.txt')




# Define this globally for access later.
def calcTime(meth):
    # Meth stands for method
    outname = 'timesteps' + meth + '.txt'
    g = open(outname, 'w')
    outstr = 'nBodies Time' + '\n'
    g.write(outstr)
    g.close()

    nt = 1500
    max_nBods = 6766

    if meth == 'bh':
        nBods = 1
        nBodsLast = 1
        while nBods <= max_nBods:
            g = open(outname, 'a')
            outs = bh.run_BarnesHut(nt, nBods)
            nBods = outs[0]
            dt = outs[1]
            g.write(str(nBods) + " " +  str(dt) + '\n')
            g.close()

            # Use Golden ratio to step (see plotSpacing.py for a comparison of how it grows compared to exponential)
            nBods += nBodsLast
            nBodsLast = nBods - nBodsLast


    if meth == 'bf':
        nBods = 1
        nBodsLast = 1
        while nBods < max_nBods:
            g = open(outname, 'a')
            outs = bf.run_BruteForce(nt, nBods)
            nBods = outs[0]
            dt = outs[1]
            g.write(str(nBods) + " " +  str(dt) + '\n')
            g.close()
            nBods += nBodsLast
            nBodsLast = nBods - nBodsLast

    else:
        print "Choose your method better (bh or bf)"








# The End
