from run_nbody_sim import bh, bf


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

