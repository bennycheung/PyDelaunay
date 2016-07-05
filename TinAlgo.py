#------------------------------------------------------------------
# TinAlgo
#
# An Incremental Algorithm for the Construction of Delauny Diagram
# author: bcheung 2013-05-07
#------------------------------------------------------------------

import math
import sys
import getopt
import time
TOL = 1e-6
BIG_FLOAT = 1e38

#------------------------------------------------------------------
# Computational Geometry functions
#------------------------------------------------------------------
"""
* Find the Circum-center of a circle which enclosed the given topological element.
* This algorithm copied from "A programmer geometry" by Adrian Bowyer p66.
* Computational expensive.
"""
def triCircumCenter(p0, p1, p2):
    xl = p0.x
    yl = p0.y
    xk = p1.x
    yk = p1.y
    xm = p2.x
    ym = p2.y

    xlk = xl - xk
    ylk = yl - yk
    xmk = xm - xk
    ymk = ym - yk
    det = xlk * ymk - xmk * ylk

    detinv = 0.5 / det
    rlksq = xlk * xlk + ylk * ylk
    rmksq = xmk * xmk + ymk * ymk
    xcc = detinv * (rlksq * ymk - rmksq * ylk)
    ycc = detinv * (xlk * rmksq - xmk * rlksq)
    return Site(xcc + xk, ycc + yk)

"""
* Returns twice the area of the oriented triangle (a, b, c), i.e., the
* area is positive if the triangle is oriented counterclockwise.
"""
def triArea(a, b, c):
    return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x)

"""
* Returns True if the point d is inside the circle defined by the
* points a, b, c. See Guibas and Stolfi (1985) p.107.
"""
def triInCircle(a, b, c, d):
    t1 = (a.x*a.x + a.y*a.y) * triArea(b, c, d)
    t2 = (b.x*b.x + b.y*b.y) * triArea(a, c, d)
    t3 = (c.x*c.x + c.y*c.y) * triArea(a, b, d)
    t4 = (d.x*d.x + d.y*d.y) * triArea(a, b, c)
    return (t1 - t2 + t3 - t4) > 0

"""
* Returns True if the points a, b, c are in a counterclockwise order
"""
def triCCW(a, b, c):
    return triArea(a, b, c) > 0

"""
* Point is on the right face of an Edge
"""
def rightOfEdge(p, e):
    return triCCW(p, e.dest(), e.org())

"""
* Point is on the left face of an Edge
"""
def leftOfEdge(p, e):
    return triCCW(p, e.org(), e.dest())

"""
* A predicate that determines if the point x is on the edge e.
* The point is considered on if it is in the EPS-neighborhood
* of the edge.
"""
def onEdge(p, e):
    xOrg = p - e.org()
    xDest = p - e.dest()
    orgDest = e.org() - e.dest();

    t1 = xOrg.length()
    t2 = xDest.length()
    if (t1 < TOL or t2 < TOL):
        return True
    t3 = orgDest.length();
    if (t1 > t3 or t2 > t3):
        return False
    line = LineEqn(e.org(), e.dest())
    return (math.fabs(line.eval(p)) < TOL);


#------------------------------------------------------------------
"""
* Site Class
*
* Data structure to represent a Site (x,y), with index
"""
class Site(object):
    def __init__(self,x=0.0,y=0.0,sitenum=0):
        self.x = x
        self.y = y
        self.sitenum = sitenum  # reference site index

    def __str__(self):
        return "(%.3f,%.3f)" % (self.x, self.y)

    def __repr__(self):
        return "(%d,%.3f,%.3f)" % (self.sitenum, self.x, self.y)

    def __eq__(self, other):
        return math.fabs(self.x - other.x) < TOL and math.fabs(self.y - other.y) < TOL

    def __add__(self, other):
        return Site(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Site(self.x - other.x, self.y - other.y)

    def length(self):
        return math.sqrt(self.x * self.x + self.y * self.y)

    def distance(self, other):
        dx = self.x - other.x
        dy = self.y - other.y
        return math.sqrt(dx*dx + dy*dy)


#------------------------------------------------------------------
"""
* SiteList Class
*
* Data structure to represent a list of Sites
"""
class SiteList(object):
    def __init__(self,pointList):
        self.__sites = []
        self.__sitenum = 0

        self.__xmin = pointList[0].x
        self.__ymin = pointList[0].y
        self.__xmax = pointList[0].x
        self.__ymax = pointList[0].y
        for i,pt in enumerate(pointList):
            self.__sites.append(Site(pt.x,pt.y,i))
            if pt.x < self.__xmin: self.__xmin = pt.x
            if pt.y < self.__ymin: self.__ymin = pt.y
            if pt.x > self.__xmax: self.__xmax = pt.x
            if pt.y > self.__ymax: self.__ymax = pt.y
        # self.__sites.sort()

    def setSiteNumber(self,site):
        site.sitenum = self.__sitenum
        self.__sitenum += 1

    class Iterator(object):
        def __init__(this,lst):  this.generator = (s for s in lst)
        def __iter__(this):      return this
        def next(this):
            try:
                return this.generator.next()
            except StopIteration:
                return None

    def iterator(self):
        return SiteList.Iterator(self.__sites)

    def __iter__(self):
        return SiteList.Iterator(self.__sites)

    def __len__(self):
        return len(self.__sites)

    def _getxmin(self): return self.__xmin
    def _getymin(self): return self.__ymin
    def _getxmax(self): return self.__xmax
    def _getymax(self): return self.__ymax
    xmin = property(_getxmin)
    ymin = property(_getymin)
    xmax = property(_getxmax)
    ymax = property(_getymax)


#------------------------------------------------------------------
"""
* Edge Class
*
* Represent an edge for the triangle.
* The following picture illustrates the result of the edge operators:
*
* lnext <----------  dest -----------> rprev
* dprev  ---------->     <-----------  dnext
*                    ^ sym
*                    | |
*                    | |
* Left               | |
* Face         <-----+-+----  rot
*       invrot  -----+-+---->        Right
*                    | |             Face
*                    | |
*                    | |
*                    e V
* onext <---------         ----------> oprev
* lprev  --------->  org  <----------  rnext
*
*
* Edge Operators Algebra:
*   e      = current edge, directed from org to dest
*   sym    = current edge, directed from dest to org
*   rot    = dual of current edge, directed from right to left
*   invrot = dual of current edge, directed from left to right
*   onext  = next ccw edge around orgin
*   oprev  = rot->onext->rot
*   dnext  = sym->onext->sym
*   dprev  = invrot->onext->invrot
*   lnext  = invrot->onext->rot
*   lprev  = onext->sym
*   rnext  = rot->onext->invrot
*   rprev  = sym->onext
*
"""
class Edge(object):

    """
    * Make a new Edge, including QuadEdge structure
    """
    @staticmethod
    def MakeEdge():
        ql = QuadEdge()
        return ql._e[0]

    """
    * This operator affects the two edge rings around the origins of a and b,
    * and, independently, the two edge rings around the left faces of a and b.
    * In each case, (i) if the two rings are distinct, sSplice will combine
    * them into one; (ii) if the two are the same ring, sSplice will break it
    * into two separate pieces.
    * Thus, sSplice can be used both to attach the two edges together, and
    * to break them apart. See Guibas and Stolfi (1985) p.96 for more details
    * and illustrations.
    """
    @staticmethod
    def Splice(a, b):
        alpha = a.onext().rot()
        beta = b.onext().rot()
        t1 = b.onext()
        t2 = a.onext()
        t3 = beta.onext()
        t4 = alpha.onext()
        a._next = t1
        b._next = t2
        alpha._next = t3;
        beta._next = t4

    """
    * Delete an Edge, reconnecting with neighbors
    """
    @staticmethod
    def DeleteEdge(e):
        Edge.Splice(e, e.oprev())
        Edge.Splice(e.sym(), e.sym().oprev())

    def __init__(self):
        self._num = 0
        self._data = None
        self._next = None
        self._qedge = None

    def __str__(self):
        return "[%s,%s]" % (self.org(), self.dest())

    def __repr__(self):
        return "[%s,%s]" % (self.org(), self.dest())

    # Return the dual of the current edge, directed from its right to its left.
    def rot(self):
        if (self._num < 3):
            return self._qedge._e[self._num + 1]
        else:
            return self._qedge._e[self._num - 3]

    # Return the dual of the current edge, directed from its left to its right.
    def invrot(self):
        if (self._num > 0):
            return self._qedge._e[self._num - 1]
        else:
            return self._qedge._e[self._num + 3]

    # Return the edge from the destination to the origin of the current edge.
    def sym(self):
        if (self._num < 2):
            return self._qedge._e[self._num + 2]
        else:
            return self._qedge._e[self._num - 2]

    # Return the next ccw edge around (from) the origin of the current edge.
    def onext(self):
        return self._next

    # Return the next cw edge around (from) the origin of the current edge.
    def oprev(self):
        return self.rot().onext().rot()

    # Return the next ccw edge around (into) the destination of the current edge.
    def dnext(self):
        return self.sym().onext().sym()

    # Return the next cw edge around (into) the destination of the current edge.
    def dprev(self):
        return self.invrot().onext().invrot()

    # Return the ccw edge around the left face following the current edge.
    def lnext(self):
        return self.invrot().onext().rot()

    # Return the ccw edge around the left face before the current edge.
    def lprev(self):
        return self.onext().sym()

    # Return the edge around the right face ccw following the current edge.
    def rnext(self):
        return self.rot().onext().invrot()

    # Return the edge around the right face ccw before the current edge.
    def rprev(self):
        return self.sym().onext()

    # Return the origin data point.
    def org(self):
        return self._data

    # Return the destination data point.
    def dest(self):
        return self.sym()._data

    # Assign the 2 end points value.
    def endPoints(self, org, dst):
        self._data = org
        self.sym()._data = dst

    # Return the quad edge data structure.
    def qedge(self):
        return _qedge

#------------------------------------------------------------------
"""
* QuadEdge Class
*
* The quad edge data structure.
* The following is the picture illustration of the quad edge data structure.
*
*             dst
*   Left Face  x     Right Face
*             ^ |sym
*             | |[2]
*             | |      rot
*      <------+-+----- [1]
*  invrot     | |
*    [3] -----+-+----->
*             | |
*           e | |
*          [0]| V
*              x
*             org
*
* Edge Description:
*   e      [0] = current edge, directed from org to dest
*   rot    [1] = dual of current edge, directed from right to left
*   sym    [2] = current edge, directed from dest to org
*   invrot [3] = dual of current edge, directed from left to right
*
"""
class QuadEdge(object):
    """
    * The QuadEdge e[4]
    """
    def __init__(self):
        self._e = 4*[None]
        self._e[0] = Edge()
        self._e[0]._num = 0
        self._e[1] = Edge()
        self._e[1]._num = 1
        self._e[2] = Edge()
        self._e[2]._num = 2
        self._e[3] = Edge()
        self._e[3]._num = 3

        self._e[0]._next = self._e[0]
        self._e[1]._next = self._e[3]
        self._e[2]._next = self._e[2]
        self._e[3]._next = self._e[1]

        self._e[0]._qedge = self
        self._e[1]._qedge = self
        self._e[2]._qedge = self
        self._e[3]._qedge = self

        self.mark = 0

#------------------------------------------------------------------
"""
* LineEqn Class
*
* Line equation support class
"""
class LineEqn(object):
    def __init__(self):
        self._a = 0
        self._b = 0
        self._c = 0

    """
    * Computes the normalized line equation through the points p and q.
    """
    def __init__(self, p, q):
        t = q - p
        len = t.length()
        self._a = t.y / len
        self._b = -t.x / len
        self._c = -(self._a*p.x + self._b*p.y)

    """
    * Plugs point p into the line equation.
    """
    def eval(self, p):
        return (self._a * p.x + self._b * p.y + self._c)

#------------------------------------------------------------------
"""
* TinBuilder Class
*
* The Triangulated Irregular Network construction class.
* The construction is based on the quad edge data structure.
* The data structure and algorithms are based on the paper:
*
* Leonidas Guibas and Jorge Stolfi,
*   "Primitives for the Manipulation of General Subdivisions and
*   the Computation of Voronoi Diagrams",
*   ACM Trans. on Graphics, Vol 4. no.2, April 1985, pp 75-123.
*
"""
class TinBuilder(object):

    """
    * Initialize a subdivision to the triangle defined by the bounding box
    * The bounding box is specified by (min X, max X) and (min Y, max Y)
    """
    def __init__(self, minx, maxx, miny, maxy):
        pt = 5*[None]
        pt[0] = Site()
        pt[0].x = minx - TOL
        pt[0].y = miny - TOL
        pt[1] = Site()
        pt[1].x = minx - TOL
        pt[1].y = maxy + TOL
        pt[2] = Site()
        pt[2].x = maxx + TOL
        pt[2].y = maxy + TOL
        pt[3] = Site()
        pt[3].x = maxx + TOL
        pt[3].y = miny - TOL
        pt[4] = Site()
        pt[4].x = minx - TOL
        pt[4].y = miny - TOL
        e = 4*[None]
        e[0] = Edge.MakeEdge()
        e[0].endPoints(pt[1], pt[0])
        for i in range(1, 4):
            e[i] = Edge.MakeEdge()
            e[i].endPoints(pt[i + 1], pt[i])
            Edge.Splice(e[i].sym(), e[i-1])
        Edge.Splice(e[3], e[0].sym())
        self._startEdge = e[0]

    """
    * Add a new edge e connecting the destination of a to the
    * origin of b, in such a way that all three have the same
    * left face after the connection is complete.
    * Additionally, the data pointers of the new edge are set.
    """
    def connect(self, a, b):
        e = Edge.MakeEdge()
        e.endPoints(a.dest(), b.org())
        Edge.Splice(e, a.lnext())
        Edge.Splice(e.sym(), b)
        return e

    """
    * Essentially turns edge e counterclockwise inside its enclosing
    * quadrilateral. The data pointers are modified accordingly.
    """
    def swap(self, e):
        a = e.oprev()
        b = e.sym().oprev()
        Edge.Splice(e, a)
        Edge.Splice(e.sym(), b)
        Edge.Splice(e, a.lnext())
        Edge.Splice(e.sym(), b.lnext())
        e.endPoints(a.dest(), b.dest())

    """
    * Returns an edge e, s.t. either x is on e, or e is an edge of
    * a triangle containing x. The search starts from startingEdge
    * and proceeds in the general direction of x. Based on the
    * pseudocode in Guibas and Stolfi (1985) p.121.
    """
    def locateSite(self, x):
        e = self._startEdge
        while True:
            if (x == e.org() or x == e.dest()):
                return e
            elif (rightOfEdge(x, e)):
                e = e.sym()
            elif (not rightOfEdge(x, e.onext())):
                e = e.onext()
            elif (not rightOfEdge(x, e.dprev())):
                e = e.dprev()
            else:
                return e

    """
    * An Incremental Algorithm for the Construction of Delauny Diagram
    *
    * Inserts a new point into a subdivision representing a Delaunay
    * triangulation, and fixes the affected edges so that the result
    * is still a Delaunay triangulation. This is based on the
    * pseudocode from Guibas and Stolfi (1985) p.120, with slight
    * modifications and a bug fix.
    """
    def insertSite(self, x):
        e = self.locateSite(x)
        if ((x == e.org()) or (x == e.dest())):
            # point already existed
            return None
        elif (onEdge(x, e)):
            e = e.oprev()
            Edge.DeleteEdge(e.onext())

        # Connect the new point to the vertices of the containing
        # triangle (or quadrilateral, if the new point fell on an
        # existing edge.)
        base = Edge.MakeEdge()
        newedge = base

        base.endPoints(e.org(), x)
        Edge.Splice(base, e)

        self._startEdge = base

        while True:
            base = self.connect(e, base.sym())
            e = base.oprev()
            if (e.lnext() == self._startEdge):
                break

        # Examine suspect edges to ensure that the Delaunay condition is satisfied.
        while True:
            t = e.oprev()
            if (rightOfEdge(t.dest(), e) and triInCircle(e.org(), t.dest(), e.dest(), x)):
                self.swap(e)
                e = e.oprev()
            elif (e.onext() == self._startEdge):
                # no more suspect edges
                return newedge
            else:
                # pop a suspect edge
                e = e.onext().lprev()


    def insertSiteList(self, siteList):
        for site in pts:
            edge = self.insertSite(site)


    def getSiteVoronoi(self, x):
        base = self.locateSite(x)
        if (not base.org(), x):
            base = base.sym()

        vpList = []
        t = base;
        while True:
            p0 = t.org()
            p1 = t.dest()
            p2 = t.lnext().dest()
            vp = triCircumCenter(p0, p1, p2)
            vpList.append(vp)
            t = t.onext()
            if (t == base):
                break
        return vpList

    def getSiteDelaunay(self, x):
        base = self.locateSite(x)
        if (not base.org(), x):
            base = base.sym()

        vpList = []
        t = base;
        while True:
            vpList.append(t.dest())
            t = t.onext()
            if (t == base):
                break
        return vpList

    def getTriangle(self, edge, tin, mark=1):
        edge._qedge.mark = mark
        # draw edge
        p0 = edge.org()
        p1 = edge.dest()
        p2 = edge.lnext().dest()
        tin.append([p0.sitenum, p1.sitenum, p2.sitenum])

        # recurse to the left face edges
        ledge= edge.onext()
        if (ledge._qedge.mark != mark):
            self.getTriangle(ledge, tin)

        redge = edge.lnext().sym()
        if (redge._qedge.mark != mark):
            self.getTriangle(redge, tin)

        # take the opposite edge
        edge = edge.sym()
        # recurse to the rightface edges
        ledge= edge.onext()
        if (ledge._qedge.mark != mark):
            self.getTriangle(ledge, tin)

        redge = edge.lnext().sym()
        if (redge._qedge.mark != mark):
            self.getTriangle(redge, tin)

    def getDelaunay(self, sites):
        timestamp = int(time.time())
        tin = []
        self.getTriangle(self._startEdge, tin, mark=timestamp)
        return tin

    def getVoronoi(self, sites):
        vor = []
        for i in range(0, len(sites)):
            site = sites[i]
            vpList = self.getSiteVoronoi(site)
            vor.append([(p.x,p.y) for p in vpList])
        return vor

