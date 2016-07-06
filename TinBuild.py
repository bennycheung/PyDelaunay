#############################################################################
#
# Voronoi diagram calculator/ Delaunay triangulator
#
# Calculate Delaunay triangulation or the Voronoi polygons for a set of
# 2D input points.
#
#############################################################################
import math
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt

from TinAlgo import *

#------------------------------------------------------------------
def usage():
    print """
TinBuild - compute Voronoi diagram or Delaunay triangulation

TinBuild [-t -v -p -d]  [filename]

TinBuild reads from filename (or standard input if no filename given) for a set
of points in the plane and writes either the Voronoi diagram or the Delaunay
triangulation to the standard output.  Each input line should consist of two
real numbers, separated by white space.

If option -t is present, the Delaunay triangulation is produced.
Each output line is a triple i j k, which are the indices of the three points
in a Delaunay triangle. Points are numbered starting at 0.

If option -v is present, the Voronoi diagram is produced.

Other options include:

d    Output textually
p    Plot graphically

"""
#------------------------------------------------------------------
class Context( object ):
    def __init__(self):
        self.doPrint = 0
        self.debug   = 0
        self.plot    = 0
        self.triangulate = False
        self.voronoi = False
        self.vertices  = []    # list of vertex 2-tuples: (x,y)
        self.triangles = []    # 3-tuple of vertex indices
        self.xmin = self.ymin = self.xmax = self.ymax = None

    def set_bounds(self,bounds):
        if not bounds == None:
            self.xmin = bounds.xmin
            self.ymin = bounds.ymin
            self.xmax = bounds.xmax
            self.ymax = bounds.ymax
        else:
            self.xmin = self.ymin = self.xmax = self.ymax = None


# Draw the delaunay triangles
def draw_delaunay(builder, sites):
    for site in sites:
        vpList = builder.getSiteDelaunay(site)
        for vp in vpList:
            plt.plot([site.x,vp.x], [site.y,vp.y], 'b')

# Draw the voronoi polygons
def draw_voronoi(builder, sites):
    for site in sites:
        vpList = builder.getSiteVoronoi(site)
        vpList.append(vpList[0])
        plt.plot([p.x for p in vpList], [p.y for p in vpList], 'b', lw=2)


def draw_triangle(edge):
    edge._qedge.mark = 1
    # draw edge
    p0 = edge.org()
    p1 = edge.dest()
    p2 = edge.lnext().dest()
    plt.plot([p0.x,p1.x], [p0.y,p1.y], 'g--')

    # recurse to the left face edges
    ledge= edge.onext()
    if (ledge._qedge.mark == 0):
        draw_triangle(ledge)

    redge = edge.lnext().sym()
    if (redge._qedge.mark == 0):
        draw_triangle(redge)

    # take the opposite edge
    edge = edge.sym()
    # recurse to the rightface edges
    ledge= edge.onext()
    if (ledge._qedge.mark == 0):
        draw_triangle(ledge)

    redge = edge.lnext().sym()
    if (redge._qedge.mark == 0):
        draw_triangle(redge)



def draw_triangles(builder, sites):
    base = builder.locateSite(sites[0])
    draw_triangle(base)


#-----------------------------------------------------------------------------
if __name__=="__main__":
    try:
        optlist,args = getopt.getopt(sys.argv[1:],"thdpv")
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    doHelp = 0
    c = Context()
    c.doPrint = 1
    for opt in optlist:
        if opt[0] == "-d":  c.debug = 1
        if opt[0] == "-p":  c.plot  = 1
        if opt[0] == "-t":  c.triangulate = 1
        if opt[0] == "-v":  c.voronoi = 1
        if opt[0] == "-h":  doHelp = 1

    if not doHelp:
        pts = []
        fp = sys.stdin
        if len(args) > 0:
            fp = open(args[0],'r')
        for line in fp:
            fld = line.split()
            x = float(fld[0])
            y = float(fld[1])
            pts.append(Site(x,y))
        if len(args) > 0: fp.close()

    if doHelp or len(pts) == 0:
        usage()
        sys.exit(2)

    # build the TIN
    sl= SiteList(pts)
    builder = TinBuilder(sl.xmin, sl.xmax, sl.ymin, sl.ymax)
    for i in range(0, len(pts)):
        pts[i].sitenum = i
        edge = builder.insertSite(pts[i])

    # use matplotlib to debug graphically
    if c.plot:
        if c.triangulate:
            draw_triangles(builder, pts)
        if c.voronoi:
            draw_voronoi(builder, pts)

        plt.plot([p.x for p in pts], [p.y for p in pts], 'ro')
        plt.axis([sl.xmin, sl.xmax, sl.ymin, sl.ymax])
        plt.show()
    else:
        for i in range(0, len(pts)):
            print "p %d %.3f %.3f" %(pts[i].sitenum, pts[i].x, pts[i].y)

        if c.triangulate:
            tin = builder.getDelaunay(pts)
            for i in range(0, len(tin)):
                (p0, p1, p2) = tin[i]
                print "t %d %d %d %d" %(i, p0, p1, p2)

        if c.voronoi:
            vor = builder.getVoronoi(pts)
            for i in range(0, len(vor)):
                plist = vor[i]
                pstr = ' '.join(["%.3f %.3f" %(p[0],p[1]) for p in plist])
                print "v %d %d %s" %(i, len(plist), pstr)

