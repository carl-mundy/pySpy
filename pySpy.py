import numpy as np
import matplotlib.pyplot as plt
import os, pprint

class tracer(object):

	def __init__(self, samprops, gavo_user=False, gavo_pwd=False, topN=20000, 
					maxz=3., minz=0., boxvol=500.**3.):
		""" Set up a SAM from GAVO in which to trace galaxy populations from

			Args:
				samprops (dict): Dictionary containing parameters for the SAM database and table
					one wishes to interact within. Must contain column name values for the following keys:
					{	'name'		:	, 	A nifty name for the SAM
						'table'		: 	,	The database..table location to look in on GAVO
						'box' 		: 	,	The simulation box this sam uses e.g. MR, MR7
						'columns'	:	,	List of galaxy property columns to get from the SAM database
						'rankby'	:	, 	Column with which rank objects by
						'rankdir'	:	, 	DESC or ASC - want to rank high to low or vice versa?
					}
				gavo_user (str): GAVO http username.
				gavo_pwd (str): GAVO http password.
				topN (int): Top N galaxies to download in the catalogues used to trace descendant
					galaxies.
		"""

		# Check all properties have been supplied
		reqcols = ['name','table','box','columns','rankby','rankdir']
		if not all([i in samprops for i in reqcols]):
			print 'SAM property dictionary incomplete - please check'

		# Store class properties
		self.samprops = samprops
		self.gavo_user = gavo_user
		self.gavo_pwd = gavo_pwd
		self.topN = topN
		self.boxvol = boxvol
		self.desc_stats_n, self.prog_stats_n = {}, {}
		self.desc_stats_c, self.prog_stats_c = {}, {}
		self.maxz, self.minz = maxz, minz
		self.name = self.samprops['name']
		self.db = self.samprops['table']
		self.props = self.samprops['columns']
		self.rankprop = self.samprops['rankby']
		self.rankdir = self.samprops['rankdir'].lower()
		self.path = './{0}/'.format(self.name.lower())
		self.dtypes = [('gid',np.int64),('did',np.int64),('lpid',np.int64)]
		[self.dtypes.append((p,np.float64)) for p in self.props] # assumes all float

		if not os.path.isdir(self.path):
			os.makedirs(self.path)
		if not os.path.isdir(self.path+'tmp/'):
			os.makedirs(self.path+'tmp/')
		if not os.path.isdir(self.path+'raw/'):
			os.makedirs(self.path+'raw/')
		if not os.path.isdir(self.path+'trace_progs/'):
			os.makedirs(self.path+'trace_progs/')

		# Generate the class cmd with/without gavo credentials
		self._mk_cmd()

		# Download the box information
		box_file = self.path+'box_{0:s}.txt'.format(self.samprops['box'].upper())
		if os.path.isfile(box_file):
			'box file'
			box_data = np.genfromtxt(box_file, names=('z','snapnum'), comments='#',)
		else:
			query = "SELECT z, snapnum FROM Snapshots..{0:s}".format(self.samprops['box'])
			cmd = self._cmd.replace('{QUERY}', query)
			os.system(cmd)
			box_data = np.genfromtxt(self.path+'tmp/query_snap.csv', names=['z','snapnum'], 
				comments='#', delimiter=',')[1:]
			np.savetxt(self.path+'box_{0:s}.txt'.format(self.samprops['box'].upper()), 
				box_data)
		self.boxdata = box_data

		# Download all the catalogues we will need
		self.get_catalogues()

	def get_catalogues(self):
		""" Download the catalogues in which to trace the progenitor and descendant
			galaxy populations
		"""
		cols2get = ','.join(self.props)
		getquery = "SELECT TOP {topn} galaxyId as gid, descendantId as did,\
				lastProgenitorId as lpid, {columns2get} FROM {db} WHERE snapnum = {snap} ORDER BY\
				{ordercol} {ordercolorder}".replace('{topn}',str(self.topN)).replace(
					'{db}',self.db).replace('{ordercol}',self.rankprop).replace(
					'{ordercolorder}', self.rankdir).replace('{columns2get}', cols2get)

		maxz_snap = self._z2snap(self.maxz)
		minz_snap = self._z2snap(self.minz)
		snaprange = np.arange(maxz_snap, minz_snap+1)

		for snap in snaprange:
			savepath = self.path+'raw/raw_{0}.npy'.format(snap)
			if os.path.isfile(savepath):	continue
			cmd = self._cmd.replace('{QUERY}', getquery.replace('{snap}',str(snap)))
			os.system(cmd)
			querydata = np.genfromtxt(self.path+'tmp/query_snap.csv', delimiter=',',
				dtype=self.dtypes)[1:]
			np.save(savepath, querydata)

	def prop_prog_c(self, prop, n, quiet=True):
		""" Calcuate the mean, median and sum of the property prop at all redshifts
			for galaxies with prop > c in the 'observed' sample.

			Args:
				prop (str): SAM property to calculate mean, median and sum.
				n (float): Initial number densities from which we derive a constant rank
					property limit which we use to observe at each redshift.
				quiet (bool): If True, no plot. If False, plot.
		"""

		if type(n) is not list:	n = [n]
		Numbers = np.array([int(self.boxvol*nn) for nn in n])

		maxz_snap = self._z2snap(self.maxz)
		minz_snap = self._z2snap(self.minz)
		snaprange = np.arange(maxz_snap, minz_snap+1)[::-1]

		means, meds, sums, recs, contams = [], [], [], [], []

		for N in Numbers:
			# Get ranking property at the initial selection
			idata = np.load(self.path+'trace_progs/prog_{0}_{1}.npy'.format(snaprange[0],N))
			c_lim = idata[prop][:N].min()

			Nmeans_o, Nmeans_t = [], []
			Nmeds_o, Nmeds_t = [],[]
			Nsums_o, Nsums_t = [], []
			Nrec, Ncontam = [], []

			for si, snap in enumerate(snaprange): # low to high z
				odata = np.load(self.path+'raw/raw_{0}.npy'.format(snap))[prop]
				tdata = np.load(self.path+'trace_progs/prog_{0}_{1}.npy'.format(snap,N))[prop]
				omask = odata >= c_lim
				tmask = tdata >= 0.

				print 'omask',np.sum(omask, dtype=float), 
				print 'tdata', np.sum(tdata >= c_lim, dtype=float)
				print 

				Nmeans_o.append(np.mean(odata[omask]))
				Nmeans_t.append(np.mean(tdata[tmask]))
				Nmeds_o.append(np.median(odata[omask]))
				Nmeds_t.append(np.median(tdata[tmask]))
				Nsums_o.append(np.sum(odata[omask]))
				Nsums_t.append(np.sum(tdata[tmask]))
				Nrec.append(np.sum(tdata >= c_lim,dtype=float) / np.sum(tdata >= 0,dtype=float))
				Ncontam.append((np.sum(omask,dtype=float) - np.sum(tdata >= c_lim,dtype=float)) / np.sum(omask,dtype=float))

			means.append([Nmeans_o, Nmeans_t])
			meds.append([Nmeds_o, Nmeds_t])
			sums.append([Nsums_o, Nsums_t])
			recs.append(Nrec)
			contams.append(Ncontam)

		self.prog_stats_c[prop] = {}
		self.prog_stats_c[prop]['z'] = self._snap2z(list(snaprange))
		self.prog_stats_c[prop]['mean'] = means
		self.prog_stats_c[prop]['median'] = meds
		self.prog_stats_c[prop]['sum'] = sums
		self.prog_stats_c[prop]['rec'] = recs
		self.prog_stats_c[prop]['contam'] = contams

		if not quiet:
			self._plot_prop(prop, n, self.prog_stats_c)

	def prop_desc_c(self, prop, n, quiet=True):
		""" Calcuate the mean, median and sum of the property prop at all redshifts
			for galaxies with prop > c in the 'observed' sample.

			Args:
				prop (str): SAM property to calculate mean, median and sum.
				n (float): Initial number densities from which we derive a constant rank
					property limit which we use to observe at each redshift.
				quiet (bool): If True, no plot. If False, plot.
		"""

		if type(n) is not list:	n = [n]
		Numbers = np.array([int(self.boxvol*nn) for nn in n])

		maxz_snap = self._z2snap(self.maxz)
		minz_snap = self._z2snap(self.minz)
		snaprange = np.arange(maxz_snap, minz_snap+1)

		means, meds, sums, recs, contams = [], [], [], [], []

		for N in Numbers:

			# Get ranking property at the initial selection
			Ndata = np.load(self.path+'traced_desc_{0}.npy'.format(N))
			idata = np.load(self.path+'raw/raw_{0}.npy'.format(snaprange[0]))
			c_lim = idata[prop][:N].min()

			Nmeans_o, Nmeans_t = [], []
			Nmeds_o, Nmeds_t = [],[]
			Nsums_o, Nsums_t = [], []
			Nrec, Ncontam = [], []

			for si, snap in enumerate(snaprange):
				odata = np.load(self.path+'raw/raw_{0}.npy'.format(snap))[prop]
				snapx = snaprange[si] - self._z2snap(self.maxz)
				tdata = Ndata[snapx][prop][:N]
				omask = odata >= c_lim
				tmask = tdata >= 0.

				Nmeans_o.append(np.mean(odata[omask]))
				Nmeans_t.append(np.mean(tdata[tmask]))
				Nmeds_o.append(np.median(odata[omask]))
				Nmeds_t.append(np.median(tdata[tmask]))
				Nsums_o.append(np.sum(odata[omask]))
				Nsums_t.append(np.sum(tdata[tmask]))
				Nrec.append(np.sum(tdata >= c_lim,dtype=float) / np.sum(tdata >= 0,dtype=float))
				Ncontam.append((np.sum(odata >= c_lim,dtype=float) - 
					np.sum(tdata >= c_lim,dtype=float)) / np.sum(odata >= c_lim,dtype=float))


			means.append([Nmeans_o, Nmeans_t])
			meds.append([Nmeds_o, Nmeds_t])
			sums.append([Nsums_o, Nsums_t])
			recs.append(Nrec)
			contams.append(Ncontam)

		self.desc_stats_c[prop] = {}
		self.desc_stats_c[prop]['z'] = self._snap2z(list(snaprange))
		self.desc_stats_c[prop]['mean'] = means
		self.desc_stats_c[prop]['median'] = meds
		self.desc_stats_c[prop]['sum'] = sums
		self.desc_stats_c[prop]['rec'] = recs
		self.desc_stats_c[prop]['contam'] = contams

		if not quiet:
			self._plot_prop(prop, n, self.desc_stats_c)

	def prop_desc_n(self, prop, n, quiet=True):
		""" Calculate the metric, as definied in Mundy et al. (2015), for a certain 
			property of descendant galaxies

			Args:
				prop (str): SAM property to calculate this metric
				n (list): List of number densities to calculate this metric at
		"""
		if type(n) is not list:	n = [n]
		Numbers = np.array([int(self.boxvol*nn) for nn in n])

		maxz_snap = self._z2snap(self.maxz)
		minz_snap = self._z2snap(self.minz)
		snaprange = np.arange(maxz_snap, minz_snap+1)

		means, meds, sums, recs, contams = [], [], [], [], []

		for N in Numbers:
			Ndata = np.load(self.path+'traced_desc_{0}.npy'.format(N))
			Nmeans_o, Nmeans_t = [], []
			Nmeds_o, Nmeds_t = [],[]
			Nsums_o, Nsums_t = [], []
			Nrec, Ncontam = [], []

			for si, snap in enumerate(snaprange[::-1]):
				odata = np.load(self.path+'raw/raw_{0}.npy'.format(snap))
				snapx = snap - self._z2snap(self.maxz)
				# Select by ranking property
				trankdata = Ndata[snapx][self.rankprop][:N]
				orankdata = odata[self.rankprop][:N]
				# Property we want to look at
				tdata = Ndata[snapx][prop][:N]
				odata = odata[prop][:N]
				omask = (odata >= 0.)
				tmask = (tdata >= 0.)

				Nmeans_o.append(np.mean(odata[omask]))
				Nmeans_t.append(np.mean(tdata[tmask]))
				Nmeds_o.append(np.median(odata[omask]))
				Nmeds_t.append(np.median(tdata[tmask]))
				Nsums_o.append(np.sum(odata[omask]))
				Nsums_t.append(np.sum(tdata[tmask]))

				# Ranking property fractions
				if self.rankdir == 
				Nrec.append(np.sum(trankdata >= orankdata.min(),dtype=float) / np.sum(tdata >= 0.,dtype=float))
				Ncontam.append((N - np.sum(trankdata >= orankdata.min(), dtype=float)) / N)

			means.append([Nmeans_o, Nmeans_t])
			meds.append([Nmeds_o, Nmeds_t])
			sums.append([Nsums_o, Nsums_t])
			recs.append(Nrec)
			contams.append(Ncontam)

		self.desc_stats_n[prop] = {}
		self.desc_stats_n[prop]['z'] = self._snap2z(list(snaprange[::-1]))
		self.desc_stats_n[prop]['mean'] = means
		self.desc_stats_n[prop]['median'] = meds
		self.desc_stats_n[prop]['sum'] = sums
		if prop == self.rankprop:
			self.desc_stats_n[prop]['rec'] = recs
			self.desc_stats_n[prop]['contam'] = contams

		if not quiet:
			self._plot_prop(prop,n,self.desc_stats_n)

	def prop_prog_n(self, prop, n, quiet=True):
		""" Calculate the mean, median and sum of a property traced via a 
			progenitor population.

			Args:
				prop (str): SAM property to trace and analyse
				n (float or list): Number densities to calculate this property at
				quiet (bool): If True, do not plot. If False, plot.
		"""
		if type(n) is not list:	n = [n]
		Numbers = np.array([int(self.boxvol*nn) for nn in n])

		maxz_snap = self._z2snap(self.maxz)
		minz_snap = self._z2snap(self.minz)
		snaprange = np.arange(maxz_snap, minz_snap+1)[::-1]

		means, meds, sums = [], [], []
		
		for N in Numbers:
			Nmeans_o, Nmeans_t = [], []
			Nmeds_o, Nmeds_t = [],[]
			Nsums_o, Nsums_t = [], []

			for si, snap in enumerate(snaprange):
				odata = np.load(self.path+'raw/raw_{0}.npy'.format(snap))[prop][:N]
				tdata = np.load(self.path+'trace_progs/prog_{0}_{1}.npy'.format(snap,N))[prop][:N]
				omask = odata >= 0
				tmask = tdata >= 0

				Nmeans_o.append(np.mean(odata[omask]))
				Nmeans_t.append(np.mean(tdata[tmask]))
				Nmeds_o.append(np.median(odata[omask]))
				Nmeds_t.append(np.median(tdata[tmask]))
				Nsums_o.append(np.sum(odata[omask]))
				Nsums_t.append(np.sum(tdata[tmask]))

			means.append([Nmeans_o, Nmeans_t])
			meds.append([Nmeds_o, Nmeds_t])
			sums.append([Nsums_o, Nsums_t])

		self.prog_stats_n[prop] = {}
		self.prog_stats_n[prop]['z'] = self._snap2z(list(snaprange))
		self.prog_stats_n[prop]['mean'] = means
		self.prog_stats_n[prop]['median'] = meds
		self.prog_stats_n[prop]['sum'] = sums 

		if not quiet:
			self._plot_prop(prop,n, self.prog_stats_n)

	def _plot_prop(self, prop, n, stats):

		# if len(n) > 1:
		fig, ax = plt.subplots(2,len(n), figsize=(3.*1.62*len(n),7),)
		# else:
		# 	# fig = plt.figure(figsize=(3.*1.62*len(n),3.5))
		# 	fig, = plt.subplots(2,1, figsize=(3.*1.62,3.5))
		# 	ax = [plt.subplot(111)]
		for ni, nd in enumerate(n):
			yo = np.array(stats[prop]['mean'][ni][0]) #obs
			yt = np.array(stats[prop]['mean'][ni][1]) #traced
			y1 = (yo-yt)/yt
			x1 = np.array(stats[prop]['z'])
			ax[0,ni].plot(x1,y1,'-',lw=3,c='indianred',label='mean')

			yo = np.array(stats[prop]['median'][ni][0]) #obs
			yt = np.array(stats[prop]['median'][ni][1]) #traced
			y1 = (yo-yt)/yt
			x1 = np.array(stats[prop]['z'])
			ax[0,ni].plot(x1,y1,'--',lw=3,c='royalblue',label='median')

			yo = np.array(stats[prop]['sum'][ni][0]) #obs
			yt = np.array(stats[prop]['sum'][ni][1]) #traced
			y1 = (yo-yt)/yt
			x1 = np.array(stats[prop]['z'])
			ax[0,ni].plot(x1,y1,'-.',lw=3,c='hotpink',label='sum')

			ax[0,ni].set_xlim(-0.1,3.1)
			ax[0,ni].set_title('{0}, {1}, n={2:1.1e}'.format(self.name.upper(),prop,nd))
			ax[0,ni].legend(loc='best',fontsize=10).draw_frame(False)
			ax[0,ni].set_xlabel('redshift')
			ax[0,ni].set_ylabel(r'$\kappa(z)$')

			if prop == self.rankprop:
				x1 = np.array(stats[prop]['z'])
				y1 = np.array(stats[prop]['rec'][ni])
				ax[1,ni].plot(x1,y1,'-',c='black',lw=3,label='recovery fraction')
				y1 = np.array(stats[prop]['contam'][ni])
				ax[1,ni].plot(x1,y1,':',c='darkred',lw=3,label='contamination fraction')
				ax[1,ni].legend(loc='best', fontsize=10).draw_frame(False)
				ax[1,ni].set_xlim(-0.1,3.1)
				ax[1,ni].set_ylim(-0.1,1.1)

		plt.tight_layout()
		plt.show()

	def trace_descendants(self, n, z=False):
		""" Trace the desendants of an initial high redshift sample to
			lower redshifts

			Args:
				n (list): list of number densities to trace at
				z (float): redshift to make the initial selection
		"""
		if not z:	z = self.maxz
		if type(n) is not list:	n = [n]
		Numbers = np.array([int(self.boxvol*nn) for nn in n])
		if (z > self.maxz) or (z < self.minz):
			print 'Initial selection redshift outside catalogue limits - please check'
		if (Numbers > self.topN).any():
			print 'Number density choices are outside limits - please check'

		Nprops = 3 + len(self.props)
		maxz_snap = self._z2snap(z)
		minz_snap = self._z2snap(self.minz)
		snaprange = np.arange(maxz_snap, minz_snap+1)

		for N in Numbers:
			savepath = self.path+'traced_desc_{0}.npy'.format(N)
			if os.path.isfile(savepath):	continue

			outdtypes = []
			for d in self.dtypes: outdtypes.append(d+tuple(([N])))
			desc_out_array = np.ones((len(snaprange)), dtype=outdtypes) # [snap][prop][gal]
			desc_out_array.fill(-99)

			for si, snap in enumerate(snaprange):
				# Go from HIGH (low) to LOW (high) REDSHIFT (snap)
				sdata = np.load(self.path+'raw/raw_{0}.npy'.format(snap))
				# sortx = np.argsort(sdata[self.rankprop])
				# if self.rankdir.lower() == 'desc':
				# 	sortx = sortx[::-1]
				# sdata = sdata[sortx]

				if snap == np.min(snaprange):
					traced_descIDs = sdata['did'][:N]
					traced_descIDs,dd,dx,mc = self._trim_dups(traced_descIDs,)
					for j in range(len(outdtypes)):
						desc_out_array[si][j][:] = sdata[outdtypes[j][0]][:N]
					skippedDescIds, skippedInds = [], []
					continue

				sortedSnapIds, setIds = sdata['gid'], {}
				for idx, gid in enumerate(sortedSnapIds):	setIds[gid] = idx
				# Look for previous snap's dids in current gids
				# ... gid if found, -2 if not found, -99 if lost/merger
				thisPos = np.array([setIds[gid] if (gid in setIds and gid >= 0) 
							else -2 if (gid not in setIds and gid >= 0)
							else -99 for gid in traced_descIDs])

				# Deal with objects that may have skipped a snapshot (typically <1%)
				skipFlag = False
				if len(skippedDescIds):
					additional_descids = []
					additional_inds = skippedInds.copy()
					skipFlag = True
					for ski, skid in enumerate(skippedDescIds):
						if (skid in setIds) and (skid >= 0):
							# if the skipped id is found in the current snap's data
							for j in range(len(outdtypes)):
								desc_out_array[si][j][skippedInds[ski]] = sdata[outdtypes[j][0]][setIds[skid]]		

							additional_descids.append(sdata['did'][setIds[skid]])
						else:
							additional_descids.append(-99)

				skippedInds = np.argwhere(thisPos == -2)[:,0]
				skippedDescIds = traced_descIDs[skippedInds]

				# Generate the new descIds for the next loop
				tmp_traced_descIDs = np.array([sdata['did'][j] if (j>=0) else -99 for j in thisPos])
				if skipFlag:	tmp_traced_descIDs[additional_inds] = additional_descids
				# Now overwrite objects that were -2 with their descIds
				if len(skippedDescIds):	tmp_traced_descIDs[skippedInds] = skippedDescIds
				for j in range(len(outdtypes)):
					mask = (thisPos >= 0)
					desc_out_array[si][outdtypes[j][0]][mask] = sdata[outdtypes[j][0]][thisPos[mask]]

				traced_descIDs = tmp_traced_descIDs
				traced_descIDs, dids,dixs,mc = self._trim_dups(traced_descIDs)

			np.save(savepath, desc_out_array)

	def trace_progenitors(self, n, z=False):
		""" Trace the progenitors of an initial low redshift sample to
			higher redshifts

			Args:
				n (list): list of number densities to trace at
				z (float): redshift to make the initil selection
		"""

		if not z:	z = self.minz
		if type(n) is not list:	n = [n]
		Numbers = np.array([int(self.boxvol*nn) for nn in n])
		if (z > self.maxz) or (z < self.minz):
			print 'Initial selection redshift outside catalogue limits - please check'
		if (Numbers > self.topN).any():
			print 'Number density choices are outside limits - please check'

		Nprops = 3 + len(self.props)
		maxz_snap = self._z2snap(self.maxz)
		minz_snap = self._z2snap(z)
		snaprange = np.arange(maxz_snap, minz_snap+1)[::-1] # low to high z

		for N in Numbers:
			outdtypes = []
			for d in self.dtypes: outdtypes.append(d+tuple(([N])))

			mindata = np.load(self.path+'raw/raw_{0}.npy'.format(minz_snap))
			topNids = mindata[['gid','lpid']][:N]
			del mindata

			for si, snap in enumerate(snaprange):
				savepath = self.path+'trace_progs/prog_{0}_{1}.npy'.format(snap,N)
				if os.path.isfile(savepath):	continue

				sdata = np.load(self.path+'raw/raw_{0}.npy'.format(snap))

				findx = np.argsort(topNids['gid'])
				findinx = np.argsort(sdata['gid'])
				sorted_top, sorted_snap = topNids[findx], sdata[findinx]

				matched_pairs = np.array(list(self._match(sorted_top, sorted_snap)))
				matched_pairs = [tuple(matched_pairs[j,:]) for j in range(len(matched_pairs))]
				matched_pairs = np.array(matched_pairs, dtype=[('gid0',np.int64)]+self.dtypes)

				# Select only the most massive progenitor galaxy in each snapshot
				unq, unq_idx	= np.unique( matched_pairs['gid0'] , return_inverse = True )
				unq_cnt			= np.bincount( unq_idx )
				cnt_mask		= (unq_cnt > 1)
				dup_ids			= np.array( unq[cnt_mask] )
				cnt_idx,		= np.nonzero( cnt_mask )
				idx_mask 		= np.in1d( unq_idx, cnt_idx )
				idx_idx,		= np.nonzero( idx_mask )
				srt_idx			= np.argsort( unq_idx[idx_mask] )
				dup_idx			= np.array( np.split( idx_idx[srt_idx], np.cumsum( unq_cnt[cnt_mask] )[:-1] ) )
				# print 'dup ids',dup_ids.shape

				keep = np.ones(len(matched_pairs), dtype=bool)
				if len(dup_idx.ravel()) > 0:
					for gidx in dup_idx:
						removeidx = gidx[np.argwhere(abs(matched_pairs[gidx][self.rankprop]) != 
								np.max(abs(matched_pairs[gidx][self.rankprop])) ).ravel()]
						keep[removeidx] = False
						del removeidx

				# Keep only the right things
				keepx = np.argwhere(keep).ravel()
				keep_progs = matched_pairs[keepx]

				np.save(savepath, keep_progs)

				del sdata, snap, si, matched_pairs, sorted_top, keep
				del keepx, keep_progs, dup_idx, unq, unq_idx, cnt_mask, cnt_idx
				del idx_mask, idx_idx, srt_idx, dup_ids

	def bootstrap_trace_desc(self, n, rankerr=0.5, z=False):
		""" Calculate the errors on things when you introduce a toy model of 
			observational errors on the ranking property.

			Args:
				n (list): list of number densities to calculate this for 
				rankerr (float): fractional error representing the sigma on a 
					gaussian from which to draw values from
				z (float): redshift to make the initial selection
		"""

	def _trim_dups(self,ids,remove=True):

		ids = np.array(ids)
		unq, unq_idx	= np.unique(ids, return_inverse=True)
		unq_cnt			= np.bincount(unq_idx)
		cnt_mask		= np.logical_and(unq_cnt > 1, unq >= 0)
		dup_ids			= np.array(unq[cnt_mask])
		cnt_idx,		= np.nonzero(cnt_mask)
		idx_mask		= np.in1d(unq_idx, cnt_idx)
		idx_idx, 		= np.nonzero(idx_mask)
		srt_idx			= np.argsort(unq_idx[idx_mask])
		dup_idx			= np.split(idx_idx[srt_idx], np.cumsum(unq_cnt[cnt_mask])[:-1])

		mergercount = 0.
		if remove:
			for i, j in enumerate(dup_ids):
				iter_inds = np.array(dup_idx[i])
				remove_inds = iter_inds[np.argwhere(iter_inds != iter_inds.min()).ravel()]
				ids[remove_inds] = -99
				mergercount += len(remove_inds)

		return (ids, dup_ids, dup_idx, mergercount)

	def _match(self, find, findin):

		# Both find and findin are required to be sorted by ascending id
		# find = array of tuples [(lo,hi),(lo1,hi2),...]
		in_iter = iter(findin)
		curr_in = next(in_iter)
		for lo, hi in find:
			while curr_in[0] < lo:
				curr_in = next(in_iter)
			while lo <= curr_in[0] <= hi:
				yield np.array([lo]+list(curr_in))
				curr_in = next(in_iter)

	def _mk_cmd(self):
		self._cmd = "wget --quiet --cookies=on --keep-session-cookies --save-cookies='{DIR}cookies.txt' --load-cookies='{DIR}cookies.txt' -O '{DIR}query_snap.csv'".replace('{DIR}',self.path+'tmp/')
		if self.gavo_user and self.gavo_pwd:
			self._cmd = self._cmd + " --http-user={0:s} --http-password={1:s}".format(self.gavo_user, self.gavo_pwd)
		self._cmd = self._cmd + " \"http://gavo.mpa-garching.mpg.de/MyMillennium?action=doQuery&SQL={QUERY}\""

	def _snap2z(self, snaps):

		if type(snaps) is not list:
			return self.boxdata['z'][np.argmin(np.abs(snaps - self.boxdata['snapnum']))]
		else:
			return [self.boxdata['z'][np.argmin(np.abs(i - self.boxdata['snapnum']))] for i in snaps]

	def _z2snap(self, redshifts):

		if type(redshifts) is not list:
			return int(self.boxdata['snapnum'][np.argmin( np.abs(redshifts - self.boxdata['z']) )])
		else:
			return [int(self.boxdata['snapnum'][np.argmin( np.abs(i - self.boxdata['z']) )]) for i in redshifts]
