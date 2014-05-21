from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt

class Dielectric:
	"""
	Class that encapsulates the data pertaining to a dielectric
	cavity. The form of the boundary, given in the implicit form
	<math>r(\theta)</math>, has to be given as a callable function. 
	Subclasses will define the particular geometry.

	We will also generate the Delaunay triangulation of the given 
	area as a function of the number of points given to the constructor.
	Subclasses will define the particular distribution of points.

	Methods:
		- 
		-

	Attributes:
		-
		-
	"""


class homoCircle:
	"""
	Class that defines the geometry and refractive index distribution
	of the circular and homogeneous cavity.
	"""
	def __init__(self, N, k, nc, no, r0):
		self.N = N
		self.k = k
		self.nc = nc
		self.no = no
		self.r0 = r0
		self.meshDistribution(N)
		self.triangulate()

	def refIndex(self, r, theta):
		if r < self.r0:
			return nc
		else:
			return no

	def boundary(self, theta):
		return self.r0

	def meshDistribution(self, N):
		N = int(np.sqrt(N))
		self.points = np.zeros((0,2))
		theta = np.linspace(0.0,2.0*np.pi,N*N)
		for i in range(N*N):
			self.points = np.insert(self.points, len(self.points), [self.boundary(theta[i]), theta[i]], axis=0)

		self.points = np.insert(self.points, len(self.points), [0,0], axis=0)

		for j in range(len(self.points)):
			r = self.points[j,0]
			theta = self.points[j,1]
			self.points[j,0] = r*np.cos(theta)
			self.points[j,1] = r*np.sin(theta)

		x = np.linspace(-1,1,N)
		y = np.linspace(-1,1,N)
		for k in range(N):
			for h in range(N):
				if np.sqrt(x[k]*x[k]+y[h]*y[h])< 1.0:
					self.points = np.insert(self.points, len(self.points), [x[k], y[h]], axis=0)


	def triangulate(self):
		self.triangulation =  Delaunay(self.points)

	def plotMesh(self):
		self.nTriangles = self.triangulation.simplices.shape[0]
		plt.triplot(self.points[:,0],self.points[:,1], self.triangulation.simplices.copy())
		plt.plot(self.points[:,0], self.points[:,1], 'o')

		# -- Compute areas
		self.areas = np.zeros((0))
		self.centerPoints = np.zeros((0,2))
		for i in range(self.nTriangles):
			# -- Areas
			a = np.linalg.norm(self.points[self.triangulation.simplices][i,0,:]-self.points[self.triangulation.simplices][i,1,:])
			b = np.linalg.norm(self.points[self.triangulation.simplices][i,1,:]-self.points[self.triangulation.simplices][i,2,:])
			c = np.linalg.norm(self.points[self.triangulation.simplices][i,2,:]-self.points[self.triangulation.simplices][i,-1,:])
			s = (a+b+c)/2
			self.areas = np.insert(self.areas, len(self.areas),np.sqrt(s*(s-a)*(s-b)*(s-c)))

			# -- Center points
			a = np.linalg.norm(self.points[self.triangulation.simplices][i,2,:]-self.points[self.triangulation.simplices][i,1,:])
			b = np.linalg.norm(self.points[self.triangulation.simplices][i,2,:]-self.points[self.triangulation.simplices][i,0,:])
			c = np.linalg.norm(self.points[self.triangulation.simplices][i,0,:]-self.points[self.triangulation.simplices][i,1,:])
			x = a*(b*b+c*c-a*a)
			y = b*(c*c+a*a-b*b)
			z = c*(a*a+b*b-c*c)
			sumXYZ = a*x+b*y+c*z
			alpha = a*x/sumXYZ
			beta = b*y/sumXYZ
			gamma = c*z/sumXYZ

			self.centerPoints = np.insert(self.centerPoints, len(self.centerPoints), 
				alpha*self.points[self.triangulation.simplices][i,0,:]
				+beta*self.points[self.triangulation.simplices][i,1,:]
				+gamma*self.points[self.triangulation.simplices][i,2,:], axis=0)

		plt.plot(self.centerPoints[:,0], self.centerPoints[:,1], 'ro')
		print(self.triangulation.simplices.shape)
		print(self.points[self.triangulation.simplices])
		plt.show()

if __name__ == '__main__':
	y = homoCircle(100, 1.0, 2.0, 1.0, 1.0)
	y.plotMesh()