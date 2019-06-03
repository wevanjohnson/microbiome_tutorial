#!/usr/bin/env python
# Initial author: Solaiappan Manimaran
# Unit Test Functions to check the Pathoscope EM algorithm functions

#	Pathoscope - Predicts strains of genomes in Nextgen seq alignment file (sam/bl8)
#	Copyright (C) 2013  Johnson Lab - Boston University
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os,sys
pathoscopedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0,pathoscopedir) 

from pathoscope.pathoid import PathoID
import unittest

class TestPathoscopeFunctions(unittest.TestCase):

	def setUp(self):
		self.maxIter = 10
		self.emEpsilon = 1e-7
		self.verbose = False
		self.piPrior = 0
		self.thetaPrior = 0

	def test_pathoscope_em_1(self):
	### 3 unique reads: 2 reads to genome1 and 1 read to genome4
		U = {0: [0, 0.5], 1: [0, 0.5]}
	### non-unique reads: 3 total reads  1:[[genomes],[qij],[xij]]
		NU = {2: [[0, 2, 3], [0.4, 0.2, 0.4] , [0.33, 0.33, 0.33], 0.4], 
			3: [[0, 2], [0.6, 0.4] , [0.5,0.5], 0.6], 4: [[0, 2], [0.5, 0.5] , [0.5,0.5], 0.5]}
	### Genome hash
		genomes = {0:"ecoli", 1:"strep", 2:"anthrax", 3:"plague"}
		(initPi, pi, theta, NU) = PathoID.pathoscope_em(U, NU, genomes, 
			self.maxIter, self.emEpsilon, self.verbose, self.piPrior, self.thetaPrior)
		print initPi, pi, theta, NU
		expectedInitPi = [0.708, 0.0, 0.228, 0.064]
		expectedPi = [1.0, 0.0, 0.0, 0.0]
		expectedTheta = [1.0, 0.0, 0.0, 0.0]
		expectedNU2 = [1.0, 0.0, 0.0]
		expectedNU3 = [1.0, 0.0]
		expectedNU4 = [1.0, 0.0]
		diffInitPi = [abs(a - b) for a, b in zip(initPi, expectedInitPi)]
		diffPi = [abs(a - b) for a, b in zip(pi, expectedPi)]
		diffTheta = [abs(a - b) for a, b in zip(theta, expectedTheta)]
		diffNU2 = [abs(a - b) for a, b in zip(NU[2][2], expectedNU2)]
		diffNU3 = [abs(a - b) for a, b in zip(NU[3][2], expectedNU3)]
		diffNU4 = [abs(a - b) for a, b in zip(NU[4][2], expectedNU4)]
		delta = [0.0001]*4
		self.assertTrue(diffInitPi < delta, "Failed initPi Assertion")
		self.assertTrue(diffPi < delta, "Failed Pi Assertion")
		self.assertTrue(diffTheta < delta, "Failed Theta Assertion")

		delta = [0.0001]*3
		self.assertTrue(diffNU2 < delta, "Failed Non-Unique Read 2 Assertion")
		delta = [0.0001]*2
		self.assertTrue(diffNU3 < delta, "Failed Non-Unique Read 3 Assertion")
		self.assertTrue(diffNU4 < delta, "Failed Non-Unique Read 4 Assertion")

	def test_pathoscope_em_2(self):
	### 3 unique reads: 2 reads to genome1 and 1 read to genome4
		U = {}
	### non-unique reads: 3 total reads  1:[[genomes],[qij],[xij]]
		NU = {0: [[0, 3], [0.5, 0.5] , [0.5, 0.5], 0.5], 
			1: [[0, 3], [0.5, 0.5] , [0.5, 0.5], 0.5], 
			2: [[0, 2, 3], [0.4, 0.2, 0.4] , [0.33, 0.33, 0.33], 0.4], 
			3: [[0, 2, 3], [0.4, 0.4, 0.2] , [0.33,0.33, 0.33], 0.4], 
			4: [[0, 2, 3], [0.4, 0.4, 0.2] , [0.33, 0.33, 0.33], 0.4]}
	### Genome hash
		genomes = {0:"ecoli", 1:"strep", 2:"anthrax", 3:"plague"}
		(initPi, pi, theta, NU) = PathoID.pathoscope_em(U, NU, genomes, self.maxIter, 
			self.emEpsilon, self.verbose, self.piPrior, self.thetaPrior)
		print initPi, pi, theta, NU
		expectedInitPi = [0.44, 0.0, 0.2, 0.36]
		expectedPi = [1.0, 0.0, 0.0, 0.0]
		expectedTheta = [1.0, 0.0, 0.0, 0.0]
		expectedNU0 = [1.0, 0.0]
		expectedNU1 = [1.0, 0.0]
		expectedNU2 = [1.0, 0.0, 0.0]
		expectedNU3 = [1.0, 0.0, 0.0]
		expectedNU4 = [1.0, 0.0, 0.0]
		diffInitPi = [abs(a - b) for a, b in zip(initPi, expectedInitPi)]
		diffPi = [abs(a - b) for a, b in zip(pi, expectedPi)]
		diffTheta = [abs(a - b) for a, b in zip(theta, expectedTheta)]
		diffNU0 = [abs(a - b) for a, b in zip(NU[2][2], expectedNU0)]
		diffNU1 = [abs(a - b) for a, b in zip(NU[3][2], expectedNU1)]
		diffNU2 = [abs(a - b) for a, b in zip(NU[2][2], expectedNU2)]
		diffNU3 = [abs(a - b) for a, b in zip(NU[3][2], expectedNU3)]
		diffNU4 = [abs(a - b) for a, b in zip(NU[4][2], expectedNU4)]
		delta = [0.01]*4
		self.assertTrue(diffInitPi < delta, "Failed initPi Assertion")
		self.assertTrue(diffPi < delta, "Failed Pi Assertion")
		self.assertTrue(diffTheta < delta, "Failed Theta Assertion")

		delta = [0.0001]*2
		self.assertTrue(diffNU0 < delta, "Failed Non-Unique Read 3 Assertion")
		self.assertTrue(diffNU1 < delta, "Failed Non-Unique Read 4 Assertion")
		delta = [0.0001]*3
		self.assertTrue(diffNU2 < delta, "Failed Non-Unique Read 2 Assertion")
		self.assertTrue(diffNU3 < delta, "Failed Non-Unique Read 2 Assertion")
		self.assertTrue(diffNU4 < delta, "Failed Non-Unique Read 2 Assertion")

	def test_pathoscope_em_3(self):
	### 3 unique reads: 2 reads to genome1 and 1 read to genome4
		U = {0: [0, 0.33], 1: [0, 0.33], 5: [0, 0.33], 6: [3, 0.33], 7: [3, 0.33]}
	### non-unique reads: 3 total reads  1:[[genomes],[qij],[xij]]
		NU = {2: [[0, 2, 3], [0.33, 0.33, 0.33] , [0.33, 0.33, 0.33], 0.33], 
			3: [[0, 2, 3], [0.33, 0.33, 0.33] , [0.33,0.33, 0.33], 0.33], 
			4: [[0, 2, 3], [0.33, 0.33, 0.33] , [0.33, 0.33, 0.33], 0.33]}
	### Genome hash
		genomes = {0:"ecoli", 1:"strep", 2:"anthrax", 3:"plague"}
		(initPi, pi, theta, NU) = PathoID.pathoscope_em(U, NU, genomes, self.maxIter, 
			self.emEpsilon, self.verbose, self.piPrior, self.thetaPrior)
		print initPi, pi, theta, NU
		expectedInitPi = [0.5, 0.0, 0.125, 0.375]
		expectedPi = [0.75, 0.0, 0.0, 0.25]
		expectedTheta = [1.0, 0.0, 0.0, 0.0]
		expectedNU2 = [1.0, 0.0, 0.0]
		expectedNU3 = [1.0, 0.0, 0.0]
		expectedNU4 = [1.0, 0.0, 0.0]
		diffInitPi = [abs(a - b) for a, b in zip(initPi, expectedInitPi)]
		diffPi = [abs(a - b) for a, b in zip(pi, expectedPi)]
		diffTheta = [abs(a - b) for a, b in zip(theta, expectedTheta)]
		diffNU2 = [abs(a - b) for a, b in zip(NU[2][2], expectedNU2)]
		diffNU3 = [abs(a - b) for a, b in zip(NU[3][2], expectedNU3)]
		diffNU4 = [abs(a - b) for a, b in zip(NU[4][2], expectedNU4)]
		delta = [0.01]*4
		self.assertTrue(diffInitPi < delta, "Failed initPi Assertion")
		self.assertTrue(diffPi < delta, "Failed Pi Assertion")
		self.assertTrue(diffTheta < delta, "Failed Theta Assertion")

		delta = [0.01]*3
		self.assertTrue(diffNU2 < delta, "Failed Non-Unique Read 2 Assertion")
		self.assertTrue(diffNU3 < delta, "Failed Non-Unique Read 3 Assertion")
		self.assertTrue(diffNU4 < delta, "Failed Non-Unique Read 4 Assertion")


if __name__ == '__main__':
	unittest.main()

