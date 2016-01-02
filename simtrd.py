#! /usr/bin/env python
"""
simtrd.py
Simulation of a sex-limited meiotic drive system with neutral modifiers and arbitrary fitness for the drive allele.
"""

from __future__ import print_function # python3 compatibility

import os
import sys
import numpy as np
import copy
from collections import Counter

import ms

## a constant-size, panmictic diploid population with mutation, recombination, selection and drift
class Population:

	def __init__(self, n = 100, theta = 0, base_tr = 0.7, s = [1.0, 1.0, 0.8, 1.0], chrlen = [100.0], responder_loc = 50.0):

		## initialize population
		self.size = n
		self.sex = np.tile([0,1], n/2)

		## set position of drive locus
		self.responder_loc = responder_loc
		self.modifiers = None
		self.base_tr = base_tr
		self.tr = np.tile(0.5, self.size)

		## set selection coefficients for drive allele
		# s is 4-vector of fitnesses for: [ aa, Aa-TRD, Aa+TRD, AA ]
		self.s = s

		## set (population-scaled) mutation rate
		self.theta = theta

		## initialize chromosomes:
		## list of 1 element per chromosome; each element is a list of arrays containing mutated sites
		self.chrlen = chrlen
		self.chroms = [ [ np.ndarray((0,), dtype = np.float32) for i in range(0, 2*n) ] for j in range(0, len(chrlen)) ]

	def __repr__(self):

		txt = "\n"
		txt += ("Population with {} individuals, each with {} chromosomes.\n\n".format(self.size, len(self.chrlen)))
		txt += ("\tMeiotic drive acts on locus at {} cM on chr1.\n".format(self.responder_loc))
		txt += ("\tChromosome lengths:\n")
		for j in range(0, len(self.chrlen)):
			txt += ("\tchr{}: {} cM\n".format(j+1, self.chrlen[j]))
			if self.modifiers is not None and len(self.modifiers):
				for m in self.modifiers:
					if m[0] == j:
						txt += ("\t-- @ {} cM: beta = {}\n".format(m[1], m[2]))

		txt += ("\n")
		f = self.get_freq(0, self.responder_loc)
		txt += ("Responder frequency: {}\n".format(f))
		if self.modifiers is not None and len(self.modifiers):
			txt += ("Modifier frequencies:\n")
			for m in self.modifiers:
				f = self.get_freq(m[0], m[1])
				txt += ("\t-- chr{} @{} cM: {}\n".format(m[0]+1, m[1], f))

		txt += ("\n")
		return txt

	def print_chroms(self):

		for j in range(0, len(self.chroms)):
			print("\n\nchr{}: ".format(j+1))
			for i in range(0, self.size):
				print("({:03d})\tA: {}".format(i, sorted(self.chroms[j][2*i])))
				print("\tB: {}".format(sorted(self.chroms[j][2*i+1])))
				print("\tsex: {}".format(self.sex[i]))
				if j == 0:
					try:
						print("\tTR: {}".format(self.tr[i]))
					except:
						pass
				print("\n")

	def as_ms(self, header = True, header_only = False):
		'''
		Format chromosomes in the output format of Hudson's ms, to facilitate summaries
		with code written expecting ms-style output.
		'''
		sites = self.get_seg_sites()
		pos = sorted(sites.keys())
		breaks = np.append([0], np.cumsum(self.chrlen))

		txt = ""
		if header or header_only:
			# fake an ms-style header
			txt = "ms {} 0 -t {}\n".format(int(2*self.size), self.theta)
			txt += "0 0 0\n\n"

		if not header_only:
			# count of segregating sites and their positions
			txt += "//\nsegsites: {}\n".format(len(sites.keys()))
			txt += "positions: {}\n".format(" ".join([ str(x) for x in pos ]))

			# now the chromsomes ('samples' in ms terminology)
			samples = [ [0]*len(pos) for i in range(0, 2*self.size) ]
			for j in range(0, len(self.chrlen)):
				on_chrom = np.nonzero(np.logical_and(pos > breaks[j], pos <= breaks[j+1]))
				lo, hi = np.min(on_chrom), np.max(on_chrom)
				for i in range(0, 2*self.size):
					for p in range(lo, hi+1):
						if pos[p] in self.chroms[j][i]:
							samples[i][p] = 1

			for s in samples:
				txt += "".join([ str(x) for x in s ]) + "\n"

		return txt

	def init_from_ms(self, ms_run):
		'''
		Initialize the population from the results of an ms run to avoid having to do
		lengthy 'burn-in' phase at beginning of each simulation.
		'''

		if not isinstance(ms_run, ms.MsReader):
			raise TypeError("Argument 'ms_sample' must be an MsRun object.")

		sys.stderr.write("Initializing population from ms sample with header:\n{}\n".format(ms_run.header))

		# read in a simulation
		ms_sample = ms_run.next()

		# initialize position of mutations
		chrlen = self.chrlen[0]
		pos = np.array(ms_sample.positions, dtype = np.float32)*chrlen # scale by chromosome length

		assert(len(ms_sample.samples) >= 2*self.size)

		for i in range(0, self.size):
			alleles = np.array([ int(x) for x in ms_sample.samples[i] ])
			derived = np.nonzero(alleles)[0]
			if len(derived):
				self.chroms[0][i] = pos[ np.nonzero(alleles)[0] ]
			else:
				self.chroms[0][i] = np.ndarray((0,), dtype = np.float32)

	def init_responder(self, freq = 0.0):
		if freq > 0.0:
			if freq < 1.0:
				ncopies = round(freq*2*self.size)
			else:
				ncopies = freq
			sys.stderr.write("There will be {} copies of the responder allele.\n".format(ncopies))
			if ncopies > 0.0:
				idx = np.random.choice(range(0, 2*self.size), replace = False, size = ncopies)
				for i in idx:
					self.chroms[0][i] = np.append(self.chroms[0][i], self.responder_loc)
			else:
				pass
		else:
			pass

	def init_modifiers(self, mods):
		# mods is a list of 4-tuples: freq, chrom, position, coefficient
		self.modifiers = [ m[1:] for m in mods ]
		for m in mods:
			freq, chrom, pos, beta = m
			if freq < 1.0:
				ncopies = round(freq*2*self.size)
			else:
				ncopies = freq
			sys.stderr.write("There will be {} copies of the modifier allele.\n".format(ncopies))
			if ncopies > 0.0:
				idx = np.random.choice(range(0, 2*self.size), replace = False, size = ncopies)
				for i in idx:
					self.chroms[chrom][i] = np.append(self.chroms[chrom][i], pos)

	def mutate(self, theta = None):

		if not theta:
			if not self.theta:
				pass
			else:
				theta = self.theta

		nsites = np.random.poisson(theta)
		#print("Doing {} mutations.".format(nsites))
		indiv = np.random.choice(range(0, 2*self.size), size = nsites)
		pos = np.random.ranf(nsites)*sum(self.chrlen)
		breaks = np.cumsum(self.chrlen)
		chroms = np.searchsorted(breaks, pos)
		for i in range(0, nsites):
			self.chroms[ chroms[i] ][ indiv[i] ] = np.append( self.chroms[ chroms[i] ][ indiv[i] ], pos[i] )
		return nsites

	def recombine(self):
		# loop on chromosomes
		for j in range(0, len(self.chrlen)):
			l = self.chrlen[j]
			d = l/100
			## loop on individuals
			# set starting strand for crossovers; random equiprobable
			strands = np.random.binomial(1, 0.5, size = self.size)
			for i in range(0, self.size):
				#print("recombining individual {}".format(i))
				# initialize crossed-over chromosomes
				n = [ np.ndarray((0,), dtype = np.float32), np.ndarray((0,), dtype = np.float32) ]
				s = [ self.chroms[j][ 2*i ], self.chroms[j][ 2*i+1 ] ]
				# draw number of crossovers
				nco = np.random.poisson(d)
				#print("there will be {} crossovers".format(nco))
				if (nco > 0):
					# draw breakpoints
					breaks = np.random.ranf(nco)*l
					breaks = np.append(np.insert(breaks, 0, 0), l)
					breaks.sort()
					#print(breaks)
					# set starting strand
					curr_strand = strands[i]
					#print("starting with strand {}".format(curr_strand))
					for b in range(1, len(breaks)):
						# each breakpoint is the right-end of an interval
						oth_strand = int(not curr_strand)
						#print("strands are [ {} ] / [ {} ]".format(curr_strand, oth_strand))
						# copy from 'current strand' to 'top chromosome', and 'other strand' to 'bottom chromosome'
						n[0] = np.append(n[0], s[curr_strand][ np.logical_and(s[curr_strand] <= breaks[b], s[curr_strand] > breaks[b-1]) ])
						n[1] = np.append(n[1], s[oth_strand][ np.logical_and(s[oth_strand] <= breaks[b], s[oth_strand] > breaks[b-1]) ])
						curr_strand = int(not curr_strand)
				else:
					# if no crossovers, just recopy parental chromosomes
					n = [ np.copy(x) for x in s ]
				# finally, update the population
				self.chroms[j][ 2*i ] = n[0]
				self.chroms[j][ 2*i+1 ] = n[1]


	def geno_at_locus(self, chrom, pos):
		geno = np.zeros(self.size)
		for i in range(0, self.size):
			geno[i] = np.sum(self.chroms[chrom][2*i] == pos) + np.sum(self.chroms[chrom][2*i+1] == pos)
		return geno

	def haplo_at_locus(self, i, chrom, pos):
		return np.array([ np.sum(self.chroms[chrom][2*i] == pos), np.sum(self.chroms[chrom][2*i+1] == pos) ])


	def calc_fitness(self):

		# start with Mendelian TR
		tr = np.tile(0.5, self.size)

		# get genotype at driver locus
		driver = self.geno_at_locus(0, self.responder_loc)

		if self.modifiers is not None and len(self.modifiers):
			# get effect sizes for modifier loci
			beta = np.array([ m[2] for m in self.modifiers ])

			# get genotypes at modifier loci
			geno = np.zeros((self.size, len(self.modifiers)))
			for i in range(0, len(self.modifiers)):
				geno[:,i] = self.geno_at_locus(self.modifiers[i][0], self.modifiers[i][1])

			# calculate TRs based on modifier loci
			tr = self.base_tr + np.dot(geno, beta)

		else:
			tr[ driver == 1 ] = self.base_tr

		# TRs can be non-Mendelian only for het females
		tr[ np.logical_or(self.sex == 0, driver != 1) ] = 0.5
		self.tr = tr

		# finally, calculate fitness
		w = np.ones(self.size)
		for i in range(0, self.size):
			if driver[i] == 0:
				w[i] = self.s[0]
			elif driver[i] == 2:
				w[i] = self.s[3]
			elif driver[i] == 1:
				if tr[i] > 0.5:
					w[i] = self.s[2]
				else:
					w[i] = self.s[1]

		self.fitness = w/np.sum(w)
		return self.fitness


	def make_offspring(self):

		offspring = copy.deepcopy(self)
		pairs = self.pick_parents()
		npairs = pairs.size/2
		#print("There are {} mating pairs.".format(npairs))

		try:

			# individuals alternate male-female
			# constant pop size, so need to make two offspring per mating
			for k in range(0, 2):

				# do recombination
				gametes = copy.deepcopy(pairs)
				#gametes.mutate()
				gametes.recombine()
				#gametes.calc_fitness()

				# loop over pairs
				for i in range(0, npairs):

					## first chromosome is subject to drive; handle it separately
					# maternal transmissions first...
					tr_mat = pairs.tr[2*i+1]
					if tr_mat > 0.5:
						which_driving = int(np.nonzero( gametes.haplo_at_locus(2*i+1, 0, gametes.responder_loc) )[0])
						non_driving = int(not which_driving)
					else:
						which_driving, non_driving = (0,1)
					probs = np.array([tr_mat, 1-tr_mat])
					choices = [which_driving, non_driving]
					xmit_mat = np.random.choice(choices, p = probs)
					# now paternal transmissions, much simpler
					xmit_pat = np.random.choice([0,1])

					# finally, draw chromosomes
					new_pat_idx = 4*i + 2*k + 0
					new_mat_idx = 4*i + 2*k + 1
					dad_idx = 4*i + 0
					mom_idx = 4*i + 2
					offspring.chroms[0][ new_pat_idx ] = gametes.chroms[0][ dad_idx + xmit_pat ]
					offspring.chroms[0][ new_mat_idx ] = gametes.chroms[0][ mom_idx + xmit_mat ]

					## loop over remaining chromosomes assuming Mendelian transmission
					for j in range(1, len(gametes.chrlen)):
						xmit_which = np.random.choice([0,1], size = 2)
						offspring.chroms[j][ new_pat_idx ] = gametes.chroms[0][ dad_idx + xmit_which[0] ]
						offspring.chroms[j][ new_mat_idx ] = gametes.chroms[0][ mom_idx + xmit_which[1] ]

				#return offspring

		except Exception as e:
			print(e)
			#return pairs

		return offspring

	def pick_parents(self):

		# compute fitness
		w = self.calc_fitness()
		# sample parents proportional to their fitness
		males = np.array(np.nonzero(self.sex == 0)).flatten()
		females = np.array(np.nonzero(self.sex == 1)).flatten()
		w_male = w[males]/np.sum(w[males])
		w_female = w[females]/np.sum(w[females])
		males = np.random.choice(males, size = self.size/2, p = w_male)
		females = np.random.choice(females, size = self.size/2, p = w_female)
		assert(len(males) == len(females))

		idx = []
		sexes = []
		for i in range(0, len(males)):
			idx.append(males[i])
			idx.append(females[i])
			sexes.append(0)
			sexes.append(1)

		# make a copy of current population
		newpop = copy.deepcopy(self)
		newpop.sex = np.array(sexes)

		# now copy chromosomes from current population to new one
		for j in range(0, len(self.chrlen)):
			for i in range(0, self.size):
				newpop.tr[i] = self.tr[ idx[i] ]
				newpop.chroms[j][2*i] = self.chroms[j][ 2*idx[i] ]
				newpop.chroms[j][2*i+1] = self.chroms[j][ 2*idx[i]+1 ]

		return newpop

	def is_fixed(self, chrom, pos):
		return self.get_freq(chrom, pos) == 1.0 or self.get_freq(chrom, pos) == 0.0

	def get_freq(self, chrom, pos):
		return np.sum(self.geno_at_locus(chrom, pos))/(2*self.size)

	def get_driver_freq(self):
		return self.get_freq(0, self.responder_loc)

	def get_seg_sites(self):
		breaks = np.cumsum(np.append([0], self.chrlen))
		sites = np.ndarray((0,), dtype = np.float32)
		for j in range(0, len(self.chrlen)):
			for i in range(0, 2*self.size):
				sites = np.append(sites, breaks[j] + self.chroms[j][i])

		counts = Counter(sites)
		return counts

## container for a simulation run with a Population object
class Trajectory:

	DRIVER_LOST = 0
	DRIVER_FIXED = 1

	def __init__(self, pop = None):

		if not isinstance(pop, Population):
			raise TypeError("Argument 'pop' should be a Population object.")

		self.pop = pop
		self.generation = 0
		self.sfs = None
		self.responder_freq = None
		self.modifier_freqs = None
		self.result = None

	def evolve(self, ngen = np.inf, check_fixed = True, verbose = False, reporter = None, **kwargs):

		if self.responder_freq is None:
			self.responder_freq = np.array( self.pop.get_driver_freq(), dtype = np.float32 )

		start_at = self.generation
		g = 0
		while not (check_fixed and self.pop.is_fixed(0, self.pop.responder_loc)) and (g < ngen):
			#print("G:{}\tdriver AF: {}".format(start_at + g, self.pop.get_driver_freq()))
			self.pop.mutate()
			self.pop = self.pop.make_offspring()
			self.responder_freq = np.append(self.responder_freq, self.pop.get_driver_freq())
			g += 1
			if not g % 100 and verbose:
				sys.stderr.write("\t... generation {}\n".format(g))
			if callable(reporter):
				reporter(self, **kwargs)

		self.generation = start_at + g
		self.sfs = self.pop.get_seg_sites()
		if self.pop.get_driver_freq() == 1.0:
			self.result = self.DRIVER_FIXED
		else:
			self.result = self.DRIVER_LOST
		return self.result


def stem_and_leaf(d):
	l,t = np.sort(d), 10
	O = range(l[0]-l[0] % t, l[-1]+11, t)
	I = np.searchsorted(l, O)
	txt = ""
	for e,a,f in zip(I,I[1:], O):
		txt += "{:3d}|{}\n".format(f/t, "".join([ str(x) for x in l[e:a]-f ]))
	return txt


if __name__ == "__main__":

	pop = Population(100, base_tr = 0.5, theta = 5, s = [1.0,1.0,0.8,1.0])

	burnin = ms.MsReader(open("init.n100.t5.out", "r"))
	pop.init_from_ms(burnin)
	pop.init_responder(1)
	sys.stderr.write(str(pop) + "\n")
	sys.stdout.write(pop.as_ms(header_only = True))
	pop.init_modifiers([ (0.01, 0, 0.5, 0.15) ]) # freq,chrom,position(cM),beta

	freqs = open("freq.txt", "w")
	stats = open("summary.txt", "w")
	print("fixed","generations","attempts", file = stats)

	#afs = []
	#def report_modifiers(run, outfile = None):
		#afs.append( (run.pop.get_freq(0, 0.5), run.pop.get_driver_freq()) )
		#if outfile:
		#	print(run.generation, run.pop.get_freq(0, 0.5), run.pop.get_driver_freq(), file = outfile)

	fix_time = []
	nfix = 0
	ntries = 0
	while nfix < 100:
		traj = Trajectory(copy.deepcopy(pop))
		rez = traj.evolve(verbose = True)
		if rez:
			sys.stderr.write("sweep #{}: {} generations (after {} runs)\n".format(nfix+1, traj.generation, ntries))
			nfix += 1
			fix_time.append(traj.generation)
			#print(stem_and_leaf(traj.pop.get_seg_sites().values()))
			print(traj.pop.as_ms(header = False))
			print(rez, traj.generation, ntries + 1, file = stats)
			#print("generation","modifier","driver", file = freqs)
			#for i in range(0, len(afs)):
			#	print(i, afs[i][0], afs[i][1], file = freqs)
			ntries = 0
		else:
			ntries += 1

	sys.stderr.write("Mean fixation time: {} generations.\n".format(np.mean(fix_time)))
