import copy
import random
import time

import numpy as np
from matplotlib import pyplot as plt

from analysis import plotTrajectories
from run_nbody_sim import axlims, gifNames, theta, AU, dt, n, mSol, rSol, mEarth, \
    rEarth, G, bh, elapsedTimeBH
from body import body


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








