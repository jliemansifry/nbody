## Overview
An n-body simulator with implementations for brute-force and Barnes-hut methods. Documentation is basically just the paper for this class (./projectWritings/draft.pdf). Run using the commands:

```python
bf.run_BruteForce(nTimesteps, nSmallBodies)
bh.run_BarnesHut(nTimesteps, nSmallBodies)
```

### Notes

- Theta is set at 0 for now; change (to 0.5 by convention, or whatever else floats your boat) if you actually want to do anything.

- Brute-force doesn't work well right now.

- Maybe a lot of the cost is coming from opening and closing files; maybe write data in arrays and then dump them all into a file at once every so often?

- Add collisional interactions: if d < r1 + r2, kill m2, add its mass to m1? or have them collide? idk

- plt.close('all') is a gift from above.

- Use twilio (imported below) to send myself a text when a run finishes.
    https://www.twilio.com/blog/2016/10/how-to-send-an-sms-with-python-using-twilio.html?utm_source=readthedocs&utm_medium=cpc&utm_campaign=python

### Additional Resources
#### General N-Body Stuff:
N-Body Simulations: http://physics.princeton.edu/~fpretori/Nbody/intro.htm

The Barnes-Hut Algorithm: http://arborjs.org/docs/barnes-hut

A hierarchical O(N log N) force-calculation algorithm (1986): http://adsabs.harvard.edu/abs/1986Natur.324..446B

Seminar Presentation: N-body Algorithms: http://www.cs.hut.fi/~ctl/NBody.pdf

REBOUND: An open-source multi-purpose N-body code for collisional dynamics: https://arxiv.org/pdf/1110.4876.pdf

N-body simulations (gravitational): http://www.scholarpedia.org/article/N-body_simulations_(gravitational)

Fast Parallel Tree Codes for Gravitational and Fluid Dynamical N-Body Problems: http://journals.sagepub.com/doi/pdf/10.1177/109434209400800205

N-Body / Particle Simulation Methods (lecture notes): http://www.cs.cmu.edu/afs/cs/academic/class/15850c-s96/www/nbody.html#pm

Class code: https://bitbucket.org/cs205nbody/cs-205-nbody/src/470572545febe95859c4c7685e06d98c27a9c2a7/Octree_common.py?at=master&fileviewer=file-view-default

#### Old N-Body Stuff:
N-Body Simulations (1975): http://adsbit.harvard.edu/cgi-bin/nph-iarticle_query?1975IAUS...69...57A&defaultprint=YES&filetype=.pdf
Algorithm for the Calculation of the Classical Equations of Motion of an N-Body System: http://aip.scitation.org/doi/pdf/10.1063/1.1664955

Galaxy Clustering Aarseth (1979): http://adsbit.harvard.edu/cgi-bin/nph-iarticle_query?1979ApJ...228..664A&defaultprint=YES&page_ind=0&filetype=.pdf

Gromacs (1995): http://sansan.phy.ncu.edu.tw/~hclee/ref/GROMACS95.pdf

New Horizon (2011): https://arxiv.org/abs/1112.1754

Cluster analysis of the nonlinear evolution of large-scale structure in an axion/gravitino/photino-dominated universe (1983): https://journals-aps-org.ezproxy.wesleyan.edu/prl/pdf/10.1103/PhysRevLett.51.935


#### Numerical Integrations:
Algorithm for numerical integration of the rigid-body equations of motion: https://journals.aps.org/pre/pdf/10.1103/PhysRevE.58.1169

Numerical Integration of the N-Body Problem: http://ccar.colorado.edu/asen5050/projects/projects_2013/Brown_Harrison/Code/Brown_H.pdf

A numerical integration scheme for the N-body gravitational problem: http://www.sciencedirect.com/science/article/pii/0021999173901605

WHFast: A fast and unbiased implementation of a symplectic Wisdom-Holman integrator for long term gravitational simulations (REBOUND integrator): https://arxiv.org/pdf/1506.01084.pdf

### Other Interesting Stuff:
A Fast Implementation and Performance Analysis of Collisionless N-body Code Based on GPGPU: http://www.sciencedirect.com/science/article/pii/S1877050912001329
