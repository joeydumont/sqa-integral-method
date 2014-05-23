from scipy.spatial import Delaunay
from scipy.special import jn, jvp, hankel1, hankel2, h1vp, h2vp
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
	def __init__(self, N, k, nc, no, r0, M):
		self.N = N
		self.k = k
		self.nc = nc
		self.no = no
		self.r0 = r0
		self.Mmax = M
		self.meshDistribution(N)
		self.triangulate()
		self.plotMesh()
		self.computeScatteringMatrix(self.Mmax)

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
			c = np.linalg.norm(self.points[self.triangulation.simplices][i,2,:]-self.points[self.triangulation.simplices][i,0,:])
			s = (a+b+c)/2.0
			self.areas = np.insert(self.areas, len(self.areas),np.sqrt(s*(s-a)*(s-b)*(s-c)))
			

			# -- Center points
			a = np.linalg.norm(self.points[self.triangulation.simplices][i,2,:]-self.points[self.triangulation.simplices][i,1,:])
			b = np.linalg.norm(self.points[self.triangulation.simplices][i,2,:]-self.points[self.triangulation.simplices][i,0,:])
			c = np.linalg.norm(self.points[self.triangulation.simplices][i,0,:]-self.points[self.triangulation.simplices][i,1,:])
			x = 1/a
			z = 1/b
			y = 1/c
			sumXYZ = a*x+b*y+c*z
			alpha = a*x/sumXYZ
			beta = b*y/sumXYZ
			gamma = c*z/sumXYZ

			self.centerPoints = np.insert(self.centerPoints, len(self.centerPoints), 
				alpha*self.points[self.triangulation.simplices][i,0,:]
				+beta*self.points[self.triangulation.simplices][i,1,:]
				+gamma*self.points[self.triangulation.simplices][i,2,:], axis=0)

		plt.plot(self.centerPoints[:,0], self.centerPoints[:,1], 'ro')
		plt.show()

	def computeScatteringMatrix(self,Mmax):
		# -- Prepare scattering matrix. 
		scatMat = np.zeros((2*10+1,2*10+1), dtype=np.complex)

		# -- Prepare the vector and matrix
		b = np.zeros((self.nTriangles), dtype=np.complex)
		M = np.zeros((self.nTriangles,self.nTriangles), dtype=np.complex)

		for i in range(self.nTriangles):
			b[i] = jn(Mmax, np.linalg.norm(self.centerPoints[i]))*np.exp(1j*Mmax*np.arctan2(self.centerPoints[i,0],self.centerPoints[i,1]))

			for j in range(self.nTriangles):
				if (i!=j):
					d = np.linalg.norm(self.centerPoints[i]-self.centerPoints[j])
					phi1 = np.arctan2(self.centerPoints[i,0],self.centerPoints[i,1])
					phi2 = np.arctan2(self.centerPoints[j,0],self.centerPoints[j,1])
					M[i,j] = 1j*hankel1(Mmax, d)*np.exp(1j*Mmax*(phi1-phi2))*self.areas[j]/4.0

		
		x = np.linalg.solve(np.eye(self.nTriangles,self.nTriangles, dtype=np.complex)-M,b)

		# -- We compute the corresponding line of the scattering matrix.
		for i in range(2*10+1):
			m = i-10
			Smm = 0.0
			for j in range(self.nTriangles):
				d = np.linalg.norm(self.centerPoints[j])
				phi = np.arctan2(self.centerPoints[j,0],self.centerPoints[j,1])
				Smm += jn(m,d)*np.exp(-1j*m*phi)*x[j]*self.areas[j]

			scatMat[i,0] = -1j*Smm/2.0
			scatMat[i,0] += (1.0 if m==0 else 0.0)
		
		print(scatMat)
		fig = plt.figure()
		plt.tripcolor(self.points[:,0], self.points[:,1],self.triangulation.simplices, np.abs(x))
		plt.show()


if __name__ == '__main__':
	y = homoCircle(700, 1.0, 2.0, 1.0, 1.0, 0)
	zc = 2.0
	zo = 1.0
	for i in range(2*10+1):
		m = i-10
		num = -(2.0*jvp(m,zc)*hankel2(m,zo)-jn(m,zc)*h2vp(m,zo))
		den = 2.0*jvp(m,zc)*hankel1(m,zo)-jn(m,zc)*h1vp(m,zo)
		print(num/den)