import csv 
import nibabel as nib
import numpy as np
from nipype.interfaces import base

### Config class defines the study-dependent configuration ###
class Config():

	def __init__(self):
		base = '/media/tabs/Kibisis/'
		self.fsSubjects = base + 'fsSubjects/'
		self.mni        = base + 'templates/mni_icbm152_t1_tal_nlin_sym_09a_masked.nii.gz'
		self.data       = base + 'omres/data/'
		self.preproc    = base + 'omres/preproc/'
		self.koffer     = base + 'omres/koffer/'
		self.rois       = self.define_rois_config(base + 'templates/')
		self.grabconf   = self.define_grabconf()
		self.tasks      = ['omres', 'localiser']
		self.designs    = self.Designs()
		self.exclude    = self.list_exclusions()

	def define_grabconf(self):

		grabconf = {
			'T1' :             
				{'sessions': False, 'extension': '.nii.gz', 'suffix': 'T1w'},
			'fmap_phase' :      
				{'sessions': True,  'extension': '.nii.gz', 'suffix': 'phasediff',  'datatype': 'fmap'},
			'fmap_mag_1' :        
				{'sessions': True,  'extension': '.nii.gz', 'suffix': 'magnitude1', 'datatype': 'fmap'},
			'fmap_mag_2' :        
				{'sessions': True,  'extension': '.nii.gz', 'suffix': 'magnitude2', 'datatype': 'fmap'},
			'fullFoV' :        
				{'sessions': True,  'extension': '.nii.gz', 'suffix': 'bold',   'task': 'fullFoV'},
			'data_omres':       
				{'sessions': True,  'extension': '.nii.gz', 'suffix': 'bold',   'task': 'omres'},
			'events_omres':     
				{'sessions': True,  'extension': '.tsv',    'suffix': 'events', 'task': 'omres'},
			'data_localiser':   
				{'sessions': True,  'extension': '.nii.gz', 'suffix': 'bold',   'task': 'localiser'},
			'events_localiser': 
				{'sessions': True,  'extension': '.tsv',    'suffix': 'events', 'task': 'localiser'}
		}

		return grabconf

	def define_rois_config(self, template_base):

		f    = template_base + 'rois/'
		rois = {'vMGBL': {'method': 'nominal',   'path':  f + 'MGB-grad1-L.nii.gz',    'size': 37},
				'vMGBR': {'method': 'nominal',   'path':  f + 'MGB-grad1-R.nii.gz',    'size': 45},
				'dMGBL': {'method': 'nominal',   'path':  f + 'MGB-rest-L.nii.gz',     'size': 77},
				'dMGBR': {'method': 'nominal',   'path':  f + 'MGB-rest-R.nii.gz',     'size': 67},
				'ICL'  : {'method': 'localiser', 'prior': f + 'combined-IC-L.nii.gz',  'size': 146},
				'ICR'  : {'method': 'localiser', 'prior': f + 'combined-IC-R.nii.gz',  'size': 146},
				'MGBL' : {'method': 'localiser', 'prior': f + 'inflated-MGB-L.nii.gz', 'size': 152},
				'MGBR' : {'method': 'localiser', 'prior': f + 'inflated-MGB-R.nii.gz', 'size': 152},
				'aHGL' : {'method': 'freesurfer', 'code': 11133},
				'aHGR' : {'method': 'freesurfer', 'code': 12133},
				'PTL'  : {'method': 'freesurfer', 'code': 11136},
				'PTR'  : {'method': 'freesurfer', 'code': 12136},
				'PPL'  : {'method': 'freesurfer', 'code': 11135},
				'PPR'  : {'method': 'freesurfer', 'code': 12135},
				'LaL'  : {'method': 'freesurfer', 'code': 11134},
				'LaR'  : {'method': 'freesurfer', 'code': 12134}}
		
		# Freesurfer cortical areas legend
		# 11133 aHG-L  G_temp_sup-G_T_transv  Anterior transverse temporal gyrus (~A1)
		# 12133 aHG-R 
		# 11136 PT-L   G_temp_sup-Plan_tempo  Planum temporale of the superior temporal gyrus (~A2)
		# 12136 PT-R  
		# 11135 PP-L   G_temp_sup-Plan_polar  Planum polare of the superior temporal gyrus
		# 12135 PP-R  
		# 11134 TSL-L  G_temp_sup-Lateral     Lateral aspect of the superior temporal gyrus
		# 12134 TSL-R 

		return rois

	def list_exclusions(self):

		exclude = [{'sub': '07', 'sess': '01', 'task': 'omres', 'run': '04'}]

		return exclude

	class Designs():

		def __init__(self):
			self.designs = {'omres':     {'picker':    self.omres,
										  'task':      'omres',
										  'volterra':  False,
										  'hrf_dervs': [0, 0],
										  'method':    'classic',
										  'smoothing': 2},
							'localiser': {'picker':    self.localiser,
										  'task':      'localiser',
										  'volterra':  False,
										  'hrf_dervs': [0, 0],
										  'method':    'classic',
										  'smoothing': 2}}

			self.contrasts = {'localiser': [('localiser', 'T', ['sound', 'silence'], [1, -1])]}

		def omres(self, logfilepath):
		 
			sounds = []
			with open(logfilepath, 'r') as logfile:
				logTsv  = csv.reader(logfile, delimiter="\t")
				next(logTsv)
				next(logTsv)
				for line in logTsv:
					if line[2][-3:] == 'st0':
						onsets = [float(line[0]) + 0.550 * n for n in range(8)]
						omloc  = int(line[2][3])
						labels, amplitudes = ['std0'], [None]
						if omloc == 0:
							labels     += ['std1' for n in range(5)]
							amplitudes += list(range(5))
							labels     += ['std2' for n in range(6, 8)]
							amplitudes += [None  for n in range(6, 8)]
						else: 
							labels     += ['std1' for n in range(omloc - 2)]
							amplitudes += list(range(omloc - 2))
							labels     += [f'om{omloc}']
							amplitudes += [None]
							labels     += ['std2' for n in range(omloc, 8)]
							amplitudes += [None  for n in range(omloc, 8)]

						for ix in range(len(onsets)):
							sounds.append({'onset': onsets[ix], 
										   'label': labels[ix], 
										   'amplitude': amplitudes[ix], 
										   'duration': 0.05})

			models = ['std0', 'std1', 'std2', 'om4', 'om5', 'om6']
			events = dict()
			for mod in models:
				events[mod] = [sound for sound in sounds if sound['label'][0:4] == mod]
			
			conditions =  [key for key in events.keys() if len(events[key]) > 0]
			onsets     =  [[ev['onset'] for ev in events[key]] for key in conditions]
			durations  =  [[ev['duration'] for ev in events[key]] for key in conditions]
			
			pmod = []
			for key in conditions:
				if key in ['std0', 'std2', 'om4', 'om5', 'om6']:
					pmod.append(None)
				else:
					pmod.append(base.Bunch(name = ['amplitude'], 
										   poly = [1], 
										   param = [[ev['amplitude'] for ev in events[key]]]))

			design = base.Bunch(conditions = conditions,
								onsets     = onsets,
								durations  = durations,
								pmod       = pmod,
								amplitudes = None,
								tmod = None)

			return design

		def localiser(self, logfilepath):

			import csv 
			from nipype.interfaces import base

			sounds, silences, keypress = [], [], []
			with open(logfilepath, 'r') as logfile:
				logTsv  = csv.reader(logfile, delimiter="\t")
				next(logTsv)
				next(logTsv)
				for line in logTsv:
					event = {'onset': float(line[0]), 'duration': float(line[1])}
					if line[2] == 'silence':
						silences.append(event)
					elif line[2] == 'key':
						keypress.append(event)
					elif event['duration'] == 1:
						sounds.append(event)
					else:
						print('WARNING: Skipping unrecognised line "{}"'.format(line))
			
			conditions = ['sound', 'silence', 'keypress']
			onsets     = [[ev['onset'] for ev in cond] for cond in [sounds, silences, keypress]]
			durations  = [[ev['duration'] for ev in cond] for cond in [sounds, silences, keypress]]
			amplitudes, tmod, pmod = [None] * 3, [None] * 3, [None] * 3

			design = base.Bunch(conditions = conditions,
								onsets     = onsets,
								durations  = durations,
								amplitudes = amplitudes,
								tmod = tmod,
								pmod = pmod)

			return design




