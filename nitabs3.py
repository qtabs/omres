# nitabs v3.0 - pybids 0.16.4

import bids
import collections
import importlib
import json
# import modelbrown  ## necessary only for bayesian analyses
import nibabel
import numpy as np
import os
import pathlib
import random
import scipy
import shutil
import string
#from icecream import ic
from nipype.interfaces import freesurfer, spm, fsl, ants, c3, nipy
from nipype.interfaces import utility, io, base
from nipype.pipeline import engine
from nipype.algorithms import modelgen, rapidart, confounds
from nipype.interfaces.utility import Function

WORKINGDIR = '/media/tabs/Kibisis/tmp'


class Koffer():

	def __init__(self, configmod, ID):

		configlib = importlib.import_module(configmod)
		self.config      = configlib.Config()
		self.ID          = ID
		self.koffer_path = os.path.join(self.config.koffer, f'{ID}')
		self.json_path   = os.path.join(self.koffer_path, f'{ID}.json')
		self.koffer      = self.load_subject_koffer()

	### Input-Output functions ###
	def update_koffer_json(self):

		shutil.copy(self.json_path, f'{self.json_path.split(".")[0]}-bu.json')
		with open(self.json_path, 'w') as f:
			json.dump(self.koffer, f, sort_keys=True, indent=4)

	def load_subject_koffer(self):

		if os.path.isfile(self.json_path):
			with open(self.json_path, 'r') as f:
				koffer = json.load(f)
		else:
			print(f'Creating new koffer for subject {self.ID}')
			koffer = dict()
			if not os.path.exists(os.path.split(self.json_path)[0]):
				os.makedirs(os.path.split(self.json_path)[0])
			with open(self.json_path, 'w') as f:
				json.dump(koffer, f, sort_keys=True, indent=4)

		return koffer

	def read(self, key):

		return self.koffer[key] if type(self.koffer[key]) is list else [self.koffer[key]]

	def get_descs(self, key, sess=None):

		if sess is None:
			fpaths = self.read(key)
		else:
			fpaths = self.read(key)[sess]

		if type(fpaths) is not list:
			fpaths = [fpaths]
		
		fnames = [os.path.split(f)[-1] for f in fpaths]
		descs  = [f[(f.find('desc-')+5):].split('.')[0] for f in fnames]

		return descs

	def check_clean_and_update_koffer(self):

		def check_clean_files(files, description):
			if type(files) is not list:
				files = [files]
			if any([not os.path.exists(file) for file in files]):
				print(f'Deleting {description} of {self.ID}: one or more files do not exist')
				for file in files:
					if os.path.exists(file):
						os.remove(file)
				return False
			else:
				return True

		remove_keys = []

		for key in self.koffer.keys():
			key_files = self.read(key)
			if type(key_files) is dict:
				for sess in key_files:
					clean = check_clean_files(key_files[sess], f'session {sess} from key {key}')
					if not clean:
						remove_keys.append((key, sess))
			if type(key_files) is list:
				clean = check_clean_files(key_files, f'key {key}')
				if not clean:
					remove_keys.append(key)

		for key in remove_keys:
			if type(key) is tuple:
				del self.koffer[key[0]][key[1]]
			elif type(key) is str:
				del self.koffer[key]

		self.update_koffer_json()

	### Validation and visualisation ###
	def freeview(self, volumes, surfaces=[], volfiles='', surfiles=''):

		if type(volumes) is not list:
			volumes = [volumes]
		
		if type(surfaces) is not list:
			surfaces = [surfaces]

		volfiles = f'{volfiles} '

		for key in volumes:
			files = self.read(key)
			if type(files) is dict:
				ff = []
				for key in files:
					if type(files[key]) is list:
						ff += files[key]
					else:
						ff += [files[key]]
				files = ff
				print(files)
			if type(files) is not list:
				files = [files]
			for f in files:
				volfiles += f' {f}'

		if surfaces != []:
			surfiles = '-f'
			for key in surfaces:
				files = self.read(key)
				if type(files) is not list:
					files = [files]
				for f in files:
					surfiles += f' {f}'

		print(f'freeview {volfiles} {surfiles}')
		os.system(f'freeview {volfiles} {surfiles}')


class Tabokoffer(Koffer):

	def __init__(self, configmod, ID):
		super().__init__(configmod, f'sub-{ID}')
		self.datapath    = self.config.data
		self.layout      = bids.layout.BIDSLayout(self.datapath)
		self.sessions    = self.get_sessions()
		self.preprocpath = self.config.preproc
		self.tasks       = self.config.tasks
		if self.preprocpath is not None:
			self.layout.add_derivatives(self.preprocpath)

	### Input-Output functions ###
	def add(self, files, key, sources, subfold, descs, sess=None, runs=False, move=True):

		def copy_move_file(src, dst):
			if not os.path.exists(os.path.join(os.path.split(dst)[0])):
				os.makedirs(os.path.join(os.path.split(dst)[0]))
			if move:
				shutil.move(src, dst)
			else:
				shutil.copy(src, dst) 	

		def compute_fname(file_path, source, desc):
			ptt = 'sub-{subject}[_task-{task}][_ses-{session}][_run-{run}]_{suffix}[_desc-{desc}]'
			entities = source.get_entities()
			entities['desc'] = desc
			bidsname  = self.layout.build_path(entities, ptt, absolute_paths=False, validate=False)
			extension = '.'.join(file_path.split('.')[1:])
			fname = os.path.join(self.koffer_path, subfold, f'{bidsname}.{extension}')
			return fname

		if runs and type(files) is str:
			files = [files]

		if not runs and sess is None:
			if type(files) is str:
				files = [files]
			stored_files = []
			if type(descs) is str:
				descs = [descs] * len(files)
			if type(sources) is not list:
				sources = [sources] * len(files)
			for file, desc, source in zip(files, descs, sources):
				dst = compute_fname(file, source, desc)
				copy_move_file(file, dst)
				stored_files.append(dst)
		elif not runs and sess is not None:
			if key in self.koffer:
				stored_files = self.read(key)
			else:
				stored_files = {}
			if type(descs) is list:
				if type(files) is not list:
					files = [files]
				stored_files[sess] = []
				for file, desc in zip(files, descs):
					dst = compute_fname(file, sources, desc)
					copy_move_file(file, dst)
					stored_files[sess].append(dst)
			else:
				dst = compute_fname(files, sources, descs)
				copy_move_file(files, dst)
				stored_files[sess] = dst
		elif runs and sess is not None:
			if key in self.koffer:
				stored_files = self.read(key)
			else:
				stored_files = {}
			stored_files[sess] = []
			for file, source in zip(files, sources):
				dst = compute_fname(file, source, descs)
				copy_move_file(file, dst)
				stored_files[sess].append(dst)

		self.koffer[key] = stored_files
		self.update_koffer_json()

	def grab_raw(self, key):

		filters = dict([(k, v) for k, v in self.config.grabconf[key].items() if k != 'sessions'])

		return self.grab(filters, sessions=self.config.grabconf[key]['sessions'])

	def grab(self, filters, sessions=False):

		basiconf = {'subject': self.ID.replace('sub-', '')}

		if sessions:
			d = dict([(s, self.layout.get(session=s, **basiconf, **filters)) for s in self.sessions])
		else:
			d = self.layout.get(**basiconf, **filters)

		return d

	def get_sessions(self):
		potential_sess = [f'{n:02d}' for n in range(99)]
		subject_n      = self.ID.replace('sub-', '')
		sessions = [s for s in potential_sess if self.layout.get(subject=subject_n, session=s)!=[]]
		return sessions

	def get_metadata(self, key, **kwargs):

		filters = {'suffix':'bold', 'desc':None, 'extension':'nii.gz'}
		filters.update(kwargs)
		
		return self.grab(filters)[0].get_metadata()[key]

	def write_bids_report(self):
		report = bids.reports.BIDSReport(self.layout).generate()

	### Preprocessing ###
	def run_fmriprep(self):
		command  = f'fmriprep-docker {self.datapath} {self.preprocpath} ' 
		command += f'--level full --output-spaces anat fun ' 
		command += f'--fs-subjects-dir {self.config.fsSubjects} -w {WORKINGDIR} '
		command += f'--fs-license-file /usr/local/freesurfer/7.4.1/license.txt'
		os.system(command)

	def store_reconall_atlas(self):
		
		fs_grabber = Tode(interface = io.FreeSurferSource(), name='fs_grabber')
		fs_grabber.inputs.subjects_dir = self.config.fsSubjects
		fs_grabber.inputs.subject_id   = self.ID
		fs_grabber.run()

		fs_a2009s = [f for f in fs_grabber.result.outputs.aparc_aseg if 'a2009s' in f][0]
		nii_conv = Tode(freesurfer.MRIConvert(out_type = 'niigz'), name = 'nii_conv')
		nii_conv.inputs.out_file = 'atlas.nii.gz' # avoids extension errors (filename has .)
		nii_conv.inputs.in_file  = fs_a2009s
		nii_conv.run()           

		nii_conv_t1 = Tode(freesurfer.MRIConvert(out_type = 'niigz'), name = 'nii_conv_T1')
		nii_conv_t1.inputs.in_file = fs_grabber.result.outputs.brain
		nii_conv_t1.run() 

		source = self.grab({'desc':None, 'suffix':'T1w', 'extension':'nii.gz'})[0]
		atlas  = nii_conv.result.outputs.out_file
		t1     = nii_conv_t1.result.outputs.out_file
		self.add(atlas, 'atlas', source, 'anat', 'atlas', move=False)
		self.add(t1,    't1',    source, 'anat', 'brain', move=False)

		fs_grabber.cleanTempFiles()
		nii_conv.cleanTempFiles()
		nii_conv_t1.cleanTempFiles()

	def filter_artifacts(self, task, confound_keys):

		confounds_tsv = self.grab({'task':task, 'desc':'confounds', 'extension':'tsv'})
		tmp_folder = create_tmp_folder(mkdir=True)
		outliers, confounds = [], []

		for f in confounds_tsv:
			df = f.get_df()

			# Confounds
			keep_cols = []
			for key in confound_keys:
				if type(key) is tuple:
					keep_cols += [col for col in list(df.columns.values) if key[0] in col][0:key[1]]
				else:
					keep_cols += [key]

			regressors = np.zeros((df.shape[0], len(keep_cols)))
			for n, key in enumerate(keep_cols):
				regressors[:,n] = df[key].values

			confounds_file = create_tmp_filename(tmp_folder, extension='txt')
			np.savetxt(confounds_file, regressors, '%.6f', delimiter='\t')
			confounds.append(confounds_file)

			# Outliers
			out_cols = [col for col in list(df.columns.values) if 'non_steady_state_outlier' in col]
			outliers_file = create_tmp_filename(tmp_folder, extension='txt')
			if out_cols == []:
				indices = []
			else:
				indices = [ix for ix, v in enumerate(df[out_cols[0]].values) if v == 1]

			np.savetxt(outliers_file, indices, '%d', delimiter='\n')
			outliers.append(outliers_file)

		return outliers, confounds, tmp_folder

	### Estimation ###
	def estimate(self, design, confound_keys):

		# Design-specific definitions
		design_picker = self.config.designs.designs[design]['picker']
		task          = self.config.designs.designs[design]['task']
		volterra      = self.config.designs.designs[design]['volterra']
		hrf_dervs     = self.config.designs.designs[design]['hrf_dervs']
		smoothing     = self.config.designs.designs[design]['smoothing']

		# Preparing data and design
		data   = self.grab({'task':task,'space': 'T1w','desc':'preproc','extension':'nii.gz'})
		events = self.grab({'task':task, 'suffix':'events'})

		outliers, confounds, confounds_tmp_folder = self.filter_artifacts(task, confound_keys)
		design_bunch   = [design_picker(l) for l in events]
		repetiton_time = self.get_metadata('RepetitionTime', task=task)
	  
		unzip = MapTode(freesurfer.MRIConvert(out_type='nii'), 'in_file', 'unzip')
		unzip.inputs.in_file = data
		unzip.run()

		if smoothing is not None:
			smoother = Tode(spm.Smooth(), name = "smooth")
			smoother.inputs.fwhm = [smoothing, smoothing, smoothing]
			smoother.inputs.in_files = unzip.result.outputs.out_file
			smoother.run()
			preprocdata = smoother.result.outputs.smoothed_files
		else:
			preprocdata = unzip.result.outputs.out_file

		# Estimating beta images
		modeler = Tode(modelgen.SpecifySPMModel(), name = 'modeler')
		modeler.inputs.concatenate_runs        = False
		modeler.inputs.input_units             = 'secs'
		modeler.inputs.output_units            = 'secs'
		modeler.inputs.time_repetition         = repetiton_time
		modeler.inputs.high_pass_filter_cutoff = 128
		modeler.inputs.subject_info            = design_bunch
		modeler.inputs.functional_runs         = preprocdata
		modeler.inputs.outlier_files           = outliers
		modeler.inputs.realignment_parameters  = confounds
		modeler.run()

		designer = Tode(spm.Level1Design(), name = 'designer')
		designer.inputs.bases                    = {'hrf': {'derivs': hrf_dervs}}
		designer.inputs.timing_units             = 'secs'
		designer.inputs.interscan_interval       = repetiton_time
		designer.inputs.session_info             = modeler.result.outputs.session_info
		designer.inputs.volterra_expansion_order = (2 if volterra else 1)
		designer.run()
			
		estimator = Tode(spm.EstimateModel(), name = 'estimator')
		estimator.inputs.estimation_method = {'Classical': 1}
		estimator.inputs.spm_mat_file      = designer.result.outputs.spm_mat_file
		estimator.run()

		beta_images_func = estimator.result.outputs.beta_images
		residuals_func   = estimator.result.outputs.residual_image
		spmmat           = estimator.result.outputs.spm_mat_file
		filtered_betas_func, descriptions = filter_betas(spmmat, beta_images_func)


		# Contrast estimation
		if design in self.config.designs.contrasts:
			contrastor = Tode(spm.EstimateContrast(), name = 'contrastor')
			contrastor.inputs.contrasts      = self.config.designs.contrasts[design]
			contrastor.inputs.spm_mat_file   = estimator.result.outputs.spm_mat_file
			contrastor.inputs.beta_images    = estimator.result.outputs.beta_images
			contrastor.inputs.residual_image = estimator.result.outputs.residual_image
			contrastor.run()


		# Resampling to T1w voxelx size and warping to MNI
		T1w  = self.read('t1')[0]
		MNI  = self.config.mni
		regt = self.grab({'from': 'T1w', 'to': 'fsnative'})[0]
		warp = self.grab({'from':'T1w', 'to':'MNI152NLin2009cAsym', 'extension':'h5'})[0]

		warping_jobs = [
			{'key': 'betas_T1w', 'reference': T1w, 'transform': regt, 'image': filtered_betas_func},
			{'key': 'betas_MNI', 'reference': MNI, 'transform': warp, 'image': 'prev_out'},
			{'key': 'res_T1w',   'reference': T1w, 'transform': regt, 'image': residuals_func},
			{'key': 'res_MNI',   'reference': MNI, 'transform': warp, 'image': 'prev_out'}]
		 
		if design in self.config.designs.contrasts:
			con_images = contrastor.result.outputs.con_images
			t_images   = contrastor.result.outputs.spmT_images
			warping_jobs += [
				{'key': 'con_T1w', 'reference': T1w, 'transform': regt, 'image': con_images},
				{'key': 'con_MNI', 'reference': MNI, 'transform': warp, 'image': 'prev_out'},
				{'key': 't_T1w',   'reference': T1w, 'transform': regt, 'image': t_images},
				{'key': 't_MNI',   'reference': MNI, 'transform': warp, 'image': 'prev_out'}]
		
		warp_out, warper, zipper = dict(), dict(), dict()
		for job in warping_jobs:
			key = job['key']
			warper[key] = MapTode(ants.ApplyTransforms(), 'input_image', name='warper_'+key)
			warper[key].inputs.input_image_type = 3
			warper[key].inputs.interpolation = 'Linear'
			warper[key].inputs.invert_transform_flags = False 
			warper[key].inputs.args = '--float'
			warper[key].inputs.num_threads = 12
			warper[key].inputs.transforms = job['transform']
			warper[key].inputs.reference_image  = job['reference']
			if job['image'] == 'prev_out':
				warper[key].inputs.input_image = prev_out
			else:
				warper[key].inputs.input_image = job['image']
			warper[key].run()

			prev_out = warper[key].result.outputs.output_image
			zipper[key] = MapTode(freesurfer.MRIConvert(out_type='niigz'), 'in_file', 'zipper_'+key)
			zipper[key].inputs.in_file = prev_out
			zipper[key].run()
			warp_out[key] = zipper[key].result.outputs.out_file


		# Storing
		spmmat = estimator.result.outputs.spm_mat_file
		self.add(spmmat, f'{design}_spmmat', data[0], 'estimation', f'{design}_spmmat')

		beta_desc = [f'beta-{d}' for d in descriptions]
		img_storing_jobs = [
			{'file_key': 'betas_T1w', 'key': 'beta',       'space': 'T1w', 'desc': beta_desc},
			{'file_key': 'betas_MNI', 'key': 'beta',       'space': 'MNI', 'desc': beta_desc},
			{'file_key': 'res_T1w',   'key': 'residuals',  'space': 'T1w', 'desc': ['residuals']},
			{'file_key': 'res_MNI',   'key': 'residuals',  'space': 'MNI', 'desc': ['residuals']}]

		if design in self.config.designs.contrasts:
			con_desc = [c[0]+'-conimages' for c in self.config.designs.contrasts[design]]
			t_desc   = [c[0]+'-timages' for c in self.config.designs.contrasts[design]]
			img_storing_jobs += [
				{'file_key': 'con_T1w', 'key': 'conimages', 'space': 'T1w', 'desc': con_desc},
				{'file_key': 'con_MNI', 'key': 'conimages', 'space': 'MNI', 'desc': con_desc},
				{'file_key': 't_T1w',   'key': 'timages',   'space': 'T1w', 'desc': t_desc},
				{'file_key': 't_MNI',   'key': 'timages',   'space': 'MNI', 'desc': t_desc}]	

		for job in img_storing_jobs:
			add_key   = f"{job['key']}_{design}_{job['space']}"
			add_descs = [f"{design}-{d}_space-{job['space']}" for d in descriptions]
			self.add(warp_out[job['file_key']], add_key, data[0], 'estimation', add_descs)

		# Cleaning up
		shutil.rmtree(confounds_tmp_folder)
		unzip.cleanTempFiles()
		if smoothing is not None:
			smoother.cleanTempFiles()
		modeler.cleanTempFiles()
		designer.cleanTempFiles()
		estimator.cleanTempFiles()
		if design in self.config.designs.contrasts:
			contrastor.cleanTempFiles()
		for key in [job['key'] for job in warping_jobs]:
			warper[key].cleanTempFiles()
			zipper[key].cleanTempFiles()

	def estimate_bayesian(self, design, confound_keys):

		# Design-specific definitions
		design_picker = self.config.designs.designs[design]['picker']
		task          = self.config.designs.designs[design]['task']
		volterra      = self.config.designs.designs[design]['volterra']
		hrf_dervs     = self.config.designs.designs[design]['hrf_dervs']
		smoothing     = self.config.designs.designs[design]['smoothing']

		# Preparing data and design
		data   = self.grab({'task':task,'space': 'T1w','desc':'preproc','extension':'nii.gz'})
		events = self.grab({'task':task, 'suffix':'events'})

		outliers, confounds, confounds_tmp_folder = self.filter_artifacts(task, confound_keys)
		
		design_bunches = dict()
		for l in events:
			for bunch in design_picker(l):
				if bunch.conditions[0] not in design_bunches:
					design_bunches[bunch.conditions[0]] = []
				design_bunches[bunch.conditions[0]].append(bunch)

		repetiton_time = self.get_metadata('RepetitionTime', task=task)
	  
		unzip = MapTode(freesurfer.MRIConvert(out_type='nii'), 'in_file', 'unzip')
		unzip.inputs.in_file = data
		unzip.run()

		# Intertim storage of model-specific derivatives
		derivatives_tmp_folder = create_tmp_folder(mkdir=True)
		dervs = ['spmmat'] + [f'{k}_{s}' for k in ['logev','cbetas'] for s in ['T1w','MNI']]
		if smoothing is not None:
			dervs += [f'{k}_smoothed_{s}' for k in ['logev','cbetas'] for s in ['T1w','MNI']]
		storing_jobs = dict([(k, {'files': [], 'descs': []}) for k in dervs])

		for bmc_mod in design_bunches:
			# Estimating beta images
			modeler = Tode(modelgen.SpecifySPMModel(), name = 'modeler')
			modeler.inputs.concatenate_runs        = False
			modeler.inputs.input_units             = 'secs'
			modeler.inputs.output_units            = 'secs'
			modeler.inputs.time_repetition         = repetiton_time
			modeler.inputs.high_pass_filter_cutoff = 128
			modeler.inputs.subject_info            = design_bunches[bmc_mod]
			modeler.inputs.functional_runs         = unzip.result.outputs.out_file
			modeler.inputs.outlier_files           = outliers
			modeler.inputs.realignment_parameters  = confounds
			modeler.run()

			designer = Tode(spm.Level1Design(), name = 'designer')
			designer.inputs.bases                    = {'hrf': {'derivs': hrf_dervs}}
			designer.inputs.timing_units             = 'secs'
			designer.inputs.interscan_interval       = repetiton_time
			designer.inputs.session_info             = modeler.result.outputs.session_info
			designer.inputs.volterra_expansion_order = (2 if volterra else 1)
			designer.run()
				
			estimator = Tode(modelbrown.EstimateBayesian(), name = 'estimator')
			estimator.inputs.write_log_evidence = True
			estimator.inputs.spm_mat_file      = designer.result.outputs.spm_mat_file
			estimator.run()

			beta_images_func = estimator.result.outputs.Cbetas
			logev_func      = estimator.result.outputs.log_evidence
			spmmat          = estimator.result.outputs.spm_mat_file
			filtered_betas_func,descriptions = filter_betas(spmmat,beta_images_func,True)

			# Smoothing
			if smoothing is not None:
				smoother_cbetas = Tode(spm.Smooth(), name = "smooth_cbetas")
				smoother_cbetas.inputs.fwhm = [smoothing, smoothing, smoothing]
				smoother_cbetas.inputs.in_files = filtered_betas_func
				smoother_cbetas.run()
				filtered_betas_smoothed_func = smoother_cbetas.result.outputs.smoothed_files
				smoother_logev = Tode(spm.Smooth(), name = "smooth_logev")
				smoother_logev.inputs.fwhm = [smoothing, smoothing, smoothing]
				smoother_logev.inputs.in_files = logev_func
				smoother_logev.run()
				logev_smoothed_func = smoother_logev.result.outputs.smoothed_files

			# Resampling to T1w voxelx size and warping to MNI
			T1w  = self.read('t1')[0]
			MNI  = self.config.mni
			regt = self.grab({'from': 'T1w', 'to': 'fsnative'})[0]
			warp = self.grab({'from':'T1w', 'to':'MNI152NLin2009cAsym', 'extension':'h5'})[0]

			warping_jobs = [
				{'key': f'cbetas_T1w', 'reference': T1w, 'transform': regt, 'image': filtered_betas_func},
				{'key': f'cbetas_MNI', 'reference': MNI, 'transform': warp, 'image': 'prev_out'},
				{'key': f'logev_T1w',  'reference': T1w, 'transform': regt, 'image': logev_func},
				{'key': f'logev_MNI',  'reference': MNI, 'transform': warp, 'image': 'prev_out'}]
			
			if smoothing is not None:
				warping_jobs += [
					{'key': f'cbetas_smoothed_T1w', 'reference': T1w, 
										'transform': regt, 'image': filtered_betas_smoothed_func},
					{'key': f'cbetas_smoothed_MNI', 'reference': MNI, 
										'transform': warp, 'image': 'prev_out'},
					{'key': f'logev_smoothed_T1w',  'reference': T1w, 
										'transform': regt, 'image': logev_smoothed_func},
					{'key': f'logev_smoothed_MNI',  'reference': MNI,
										'transform': warp, 'image': 'prev_out'}]

			warp_out, warper, zipper = dict(), dict(), dict()
			for job in warping_jobs:
				key = job['key']
				warper[key] = MapTode(ants.ApplyTransforms(), 'input_image', name='warper_'+key)
				warper[key].inputs.input_image_type = 3
				warper[key].inputs.interpolation = 'Linear'
				warper[key].inputs.invert_transform_flags = False 
				warper[key].inputs.args = '--float'
				warper[key].inputs.num_threads = 12
				warper[key].inputs.transforms = job['transform']
				warper[key].inputs.reference_image  = job['reference']
				if job['image'] == 'prev_out':
					warper[key].inputs.input_image = prev_out
				else:
					warper[key].inputs.input_image = job['image']
				warper[key].run()

				prev_out = warper[key].result.outputs.output_image
				zipper[key] = MapTode(freesurfer.MRIConvert(out_type='niigz'), 'in_file', 'zipper_'+key)
				zipper[key].inputs.in_file = prev_out
				zipper[key].run()
				warp_out[key] = zipper[key].result.outputs.out_file

			# Storing
			preffix = f'bmc-{bmc_mod}_'
			mod_tmp_fold = os.path.join(derivatives_tmp_folder, bmc_mod)
			os.mkdir(mod_tmp_fold)

			dst = os.path.join(mod_tmp_fold, os.path.split(spmmat)[1])
			shutil.copy(spmmat, dst)
			storing_jobs['spmmat']['files'].append(dst)
			storing_jobs['spmmat']['descs'].append(f'{preffix}spmmat')

			for key in warp_out:
				files = warp_out[key] if type(warp_out[key]) is list else [warp_out[key]]
				dsts  = [os.path.join(mod_tmp_fold, os.path.split(f)[1]) for f in files]
				for src, dst in zip(files, dsts):
					shutil.copy(src, dst)

				suffix = '_smoothed' if 'smoothed' in key else '' 

				if 'logev' in key:
					descs = [f'{preffix}logev{suffix}']
				elif 'beta' in key:
					descs = [f'{preffix}cbeta-{d}{suffix}' for d in descriptions]

				storing_jobs[key]['files'] += dsts
				storing_jobs[key]['descs'] += descs

			# Clean up
			modeler.cleanTempFiles()
			designer.cleanTempFiles()
			estimator.cleanTempFiles()
			if smoothing is not None:
				smoother_logev.cleanTempFiles()
				smoother_cbetas.cleanTempFiles()
			for key in [job['key'] for job in warping_jobs]:
				warper[key].cleanTempFiles()
				zipper[key].cleanTempFiles()

		for key in storing_jobs:
			files = storing_jobs[key]['files']
			if key == 'spmmat':
				descs = storing_jobs[key]['descs']
			else:
				space = 'T1w' if 'T1w' in key else 'MNI' 
				descs = [f'{d}_space-{space}' for d in storing_jobs[key]['descs']]
			self.add(files, f'bmc_{key}', data[0], 'bmc', descs, move=False)

		shutil.rmtree(confounds_tmp_folder)
		shutil.rmtree(derivatives_tmp_folder)
		unzip.cleanTempFiles()

	def create_analysis_rois(self):

		base_dir = create_tmp_folder(mkdir=True)
		warpers = []

		T1w = self.read('t1')[0]
		warp_inv = self.grab({'from':'MNI152NLin2009cAsym', 'to':'T1w', 'extension':'h5'})[0]
		warp_dir = self.grab({'from':'T1w', 'to':'MNI152NLin2009cAsym', 'extension':'h5'})[0]

		# T1w space
		analysis_roi_paths, sources, descs = [], [], []
		for key, roi_config in self.config.rois.items():

			if roi_config['method'] in ['nominal', 'localiser']:

				# Warping prior/mask to T1 space
				path_key = 'path' if (roi_config['method'] == 'nominal') else 'prior'

				warpers.append(Tode(ants.ApplyTransforms(), 'warper'))
				warpers[-1].inputs.interpolation   = 'Linear'
				warpers[-1].inputs.num_threads     = 12
				warpers[-1].inputs.reference_image = T1w
				warpers[-1].inputs.input_image     = roi_config[path_key]
				warpers[-1].inputs.transforms      = warp_inv
				warpers[-1].run()

				roi_contrast = nibabel.load(warpers[-1].result.outputs.output_image).get_fdata()
				roi_affine   = nibabel.load(warpers[-1].result.outputs.output_image).affine

				# Masking prior mask
				if roi_config['method'] == 'localiser':
					con_images = self.read('conimages_localiser_T1w')
					summed_contrast = sum([nibabel.load(f).get_fdata() for f in con_images])
					roi_contrast = (roi_contrast > 0).astype(float) * summed_contrast

				# Thresholding
				threshold = 0
				for step in [1, 0.1, 0.01, 0.001, 0.0001, 0.00001]:
					while (roi_contrast > threshold).sum() > roi_config['size']:
						threshold += step
					threshold -= step
				threshold += step

				mask_array = (roi_contrast > threshold).astype(float)

			elif roi_config['method'] == 'freesurfer':
				atlas_array  = nibabel.load(self.read('atlas')[0]).get_fdata()
				atlas_affine = nibabel.load(self.read('atlas')[0]).affine
				mask_array   = (atlas_array==roi_config['code']).astype(float)

			if len(mask_array.shape) == 4:
				mask_array = mask_array[..., 0]

			mask_path = create_tmp_filename(base_dir, extension='nii.gz')
			nibabel.save(nibabel.Nifti1Image(mask_array, roi_affine), mask_path)

			if roi_config['method'] in ['nominal', 'freesurfer']:
				mask_source = self.grab_raw('T1')[0]
			elif roi_config['method'] in ['localiser']:
				mask_source = self.grab_raw('data_localiser')[self.sessions[0]][0]

			analysis_roi_paths.append(mask_path)
			sources.append(mask_source)
			descs.append(f'roi-{key}_space-T1w')

		self.add(analysis_roi_paths, 'rois_T1w', sources, 'anat', descs)
		
		# MNI space
		analysis_roi_paths, sources, descs = [], [], []
		for key, roi_config in self.config.rois.items():

			if roi_config['method'] == 'nominal':
				mask_path = roi_config['path']

			elif roi_config['method'] == 'localiser':
				roi_contrast = nibabel.load(roi_config['prior']).get_fdata()
				roi_affine   = nibabel.load(roi_config['prior']).affine

				# Masking prior mask
				if roi_config['method'] == 'localiser':
					con_images = self.read('conimages_localiser_MNI')
					summed_contrast = sum([nibabel.load(f).get_fdata() for f in con_images])
					roi_contrast = (roi_contrast > 0).astype(float) * summed_contrast

				# Thresholding
				threshold = 0
				for step in [1, 0.1, 0.01, 0.001, 0.0001, 0.00001]:
					while (roi_contrast > threshold).sum() > roi_config['size']:
						threshold += step
					threshold -= step
				threshold += step

				mask_array = (roi_contrast > threshold).astype(float)
				mask_path = create_tmp_filename(base_dir, extension='nii.gz')
				nibabel.save(nibabel.Nifti1Image(mask_array, roi_affine), mask_path)
			
			# ToDo: Save freesurfer parcelations in MNI by segmenting the template in surface space
			elif roi_config['method'] == 'freesurfer':
				atlas_array  = nibabel.load(self.read('atlas')[0]).get_fdata()
				atlas_affine = nibabel.load(self.read('atlas')[0]).affine
				mask_array   = (atlas_array==roi_config['code']).astype(float)
				tmp_img_path = create_tmp_filename(base_dir, extension='nii.gz')
				nibabel.save(nibabel.Nifti1Image(mask_array, atlas_affine), tmp_img_path)

				warpers.append(Tode(ants.ApplyTransforms(), 'warper'))
				warpers[-1].inputs.interpolation   = 'Linear'
				warpers[-1].inputs.num_threads     = 12
				warpers[-1].inputs.reference_image = self.config.mni
				warpers[-1].inputs.input_image     = tmp_img_path
				warpers[-1].inputs.transforms      = warp_dir
				warpers[-1].run()
				mask_path = warpers[-1].result.outputs.output_image

			if len(mask_array.shape) == 4:
				mask_array = mask_array[..., 0]

			if roi_config['method'] in ['nominal', 'freesurfer']:
				mask_source = self.grab_raw('T1')[0]
			elif roi_config['method'] in ['localiser']:
				mask_source = self.grab_raw('data_localiser')[self.sessions[0]][0]

			analysis_roi_paths.append(mask_path)
			sources.append(mask_source)
			descs.append(f'roi-{key}_space-MNI')

		self.add(analysis_roi_paths, 'rois_MNI', sources, 'anat', descs, move=False)

		shutil.rmtree(base_dir)
		for warper in warpers:
			warper.cleanTempFiles()

	def estimate_residuals_sigma(self, design, confound_keys):

		design_picker = self.config.designs.designs[design]['picker']
		task          = self.config.designs.designs[design]['task']
		volterra      = self.config.designs.designs[design]['volterra']
		hrf_dervs     = self.config.designs.designs[design]['hrf_dervs']
		smoothing     = self.config.designs.designs[design]['smoothing']

		data   = self.grab({'task':task,'space':'T1w','desc':'preproc','extension':'nii.gz'})
		events = self.grab({'task':task,'suffix':'events'})
		outliers, confound_series, confounds_tmp_folder = self.filter_artifacts('localiser', confound_keys)

		design_bunch   = [design_picker(l) for l in events]
		repetiton_time = self.get_metadata('RepetitionTime', task=task)

		unzip = MapTode(freesurfer.MRIConvert(out_type='nii'), 'in_file', 'unzip')
		unzip.inputs.in_file = data
		unzip.run()
		preprocdata = unzip.result.outputs.out_file

		# Estimating residuals
		modeler = Tode(modelgen.SpecifySPMModel(), name = 'modeler')
		modeler.inputs.concatenate_runs        = False
		modeler.inputs.input_units             = 'secs'
		modeler.inputs.output_units            = 'secs'
		modeler.inputs.time_repetition         = repetiton_time
		modeler.inputs.high_pass_filter_cutoff = 128
		modeler.inputs.subject_info            = design_bunch
		modeler.inputs.functional_runs         = preprocdata
		modeler.inputs.outlier_files           = outliers
		modeler.inputs.realignment_parameters  = confound_series
		modeler.run()

		designer = Tode(spm.Level1Design(), name = 'designer')
		designer.inputs.bases                    = {'hrf': {'derivs': hrf_dervs}}
		designer.inputs.timing_units             = 'secs'
		designer.inputs.interscan_interval       = repetiton_time
		designer.inputs.session_info             = modeler.result.outputs.session_info
		designer.inputs.volterra_expansion_order = (2 if volterra else 1)
		designer.run()

		estimator = Tode(spm.EstimateModel(), name = 'estimator')
		estimator.inputs.estimation_method = {'Classical': 1}
		estimator.inputs.spm_mat_file      = designer.result.outputs.spm_mat_file
		estimator.inputs.write_residuals   = True
		estimator.run()

		# Computing variance of the residuals
		tsrnr = Tode(confounds.TSNR(), name = 'tsnr_computer')
		tsrnr.inputs.in_file = estimator.result.outputs.residual_images
		tsrnr.run()
		res_sigma = tsrnr.result.outputs.stddev_file

		# Resampling to T1w voxelx size and warping to MNI
		T1w  = self.read('t1')[0]
		MNI  = self.config.mni
		regt = self.grab({'from': 'T1w', 'to': 'fsnative'})[0]
		warp = self.grab({'from':'T1w', 'to':'MNI152NLin2009cAsym', 'extension':'h5'})[0]

		warping_jobs = [
			{'key': 'res-sigma_T1w', 'reference': T1w, 'transform': regt, 'image': res_sigma},
			{'key': 'res-sigma_MNI', 'reference': MNI, 'transform': warp, 'image': 'prev_out'}]
				 
		warp_out, warper, zipper = dict(), dict(), dict()
		for job in warping_jobs:
			key = job['key']
			warper[key] = MapTode(ants.ApplyTransforms(), 'input_image', name='warper_'+key)
			warper[key].inputs.input_image_type = 3
			warper[key].inputs.interpolation = 'Linear'
			warper[key].inputs.invert_transform_flags = False 
			warper[key].inputs.args = '--float'
			warper[key].inputs.num_threads = 12
			warper[key].inputs.transforms = job['transform']
			warper[key].inputs.reference_image  = job['reference']
			if job['image'] == 'prev_out':
				warper[key].inputs.input_image = prev_out
			else:
				warper[key].inputs.input_image = job['image']
			warper[key].run()

			prev_out = warper[key].result.outputs.output_image
			zipper[key] = MapTode(freesurfer.MRIConvert(out_type='niigz'), 'in_file', 'zipper_'+key)
			zipper[key].inputs.in_file = prev_out
			zipper[key].run()
			warp_out[key] = zipper[key].result.outputs.out_file

		# Storing
		for key in warp_out:
			desc = f'res-sigma-{design}_space-{key.split("_")[1]}'
			self.add(warp_out[key], f'{key}-{design}', data[0], 'estimation', desc)

		# Cleaning up
		shutil.rmtree(confounds_tmp_folder)
		unzip.cleanTempFiles()
		modeler.cleanTempFiles()
		designer.cleanTempFiles()
		estimator.cleanTempFiles()
		tsrnr.cleanTempFiles()

		for key in [job['key'] for job in warping_jobs]:
			warper[key].cleanTempFiles()
			zipper[key].cleanTempFiles()

	### Export ###
	def export_fields(self, keys, folder='./export'):

		if not os.path.exists(folder):
			os.makedirs(folder)

		for k in keys:
			subfold = os.path.join(folder, self.read(k)[0].split('/')[-2]) + '/'
			if not os.path.exists(subfold):
				os.makedirs(subfold)
			for f in self.read(k):
				shutil.copy(f, subfold)

	def extract_betas(self, design, space='T1w', conditions=None, rois='all'):

		masks = self.extract_roi_masks(space, rois)
		session_betas = self.read(f'beta_{design}_{space}')
		session_names = self.get_descs(f'beta_{design}_{space}')
		
		betas = dict()
		for b, n in zip(session_betas, session_names):
			key = n[:n.find('_run')].replace('beta-', '').replace(f'{design}-', '')
			if conditions is None or key in conditions:
				if key not in betas:
					betas[key] = []
				betas[key].append(b)

		beta_dists = dict([(r, dict([(b, []) for b in betas])) for r in masks])
		for beta_name in betas:
			for b in betas[beta_name]:
				beta_array = nibabel.load(b).get_fdata()
				for roi_name in masks:
					beta_dists[roi_name][beta_name].append(beta_array[masks[roi_name]].reshape(-1))

		return beta_dists

	def extract_zscores(self, design, space='T1w', conditions=None, rois='all'):

		beta_dists = self.extract_betas(design, space, conditions, rois)
		zscores    = dict([(roi, dict([(cond, []) for cond in conditions])) for roi in beta_dists])

		for roi in beta_dists:
			for run_n in range(len(list(beta_dists['ICL'].values())[0])):
				roi_distribution = []
				for condition in beta_dists[roi]:
					roi_distribution.extend(beta_dists[roi][condition][run_n])
				mean, std = np.mean(roi_distribution), np.std(roi_distribution)
				if std > 1E-5: # Avoiding runs where the ROI was not in the slab
					for condition in conditions:
						zscores[roi][condition].append((beta_dists[roi][condition][run_n]-mean)/std)
				else:
					print(f'Excluding roi {roi} of {self.ID} (run-{run_n}))', end=' ')
					print(f'with (mu={mean:.2f}, sigma={std:.2f})')

		return zscores

	def extract_funcloc_contrasts(self, space='T1w', rois='all'):

		masks = self.extract_roi_masks(space, rois)
		
		con_images = self.read(f'conimages_localiser_{space}')
		global_contrast = sum([nibabel.load(f).get_fdata() for f in con_images])

		contrasts = dict([(r, []) for r in masks])
		for roi_name in masks:
			contrasts[roi_name] = global_contrast[masks[roi_name]]

		return contrasts

	def extract_logevs(self, design, space='T1w', rois=None):

		if rois is not None:
			masks = self.extract_roi_masks(space, rois)

		files  = self.read(f'{design}_logev_{space}')
		models = [f[f.find(f'-{design}-')+len(f'-{design}-'):f.find('_logev')] for f in files]

		if rois is None:
			logevs = dict([(m, []) for m in models])
		else:
			logevs = dict([(r, dict([(m, []) for m in models])) for r in masks])

		for m, f in zip(models, files):
			logev_array = nibabel.load(f).get_fdata()
			if rois is None:
				logevs[m] = logev_array
			else:
				for r in masks:
					logevs[r][m] = logev_array[masks[r]].reshape(-1)

		return logevs

	def extract_roi_masks(self, space='T1w', rois='all'):

		masks = dict()
		for mask_image, desc in zip(self.read(f'rois_{space}'), self.get_descs(f'rois_{space}')):
			key = [d.replace('roi-', '') for d in desc.split('_') if 'roi-' in d][0]
			if rois == 'all' or key in rois:
				masks[key] = (nibabel.load(mask_image).get_fdata() == 1)

		return masks

	def extract_residuals_sigma(self, design, space='T1w', rois='all'):

		masks     = self.extract_roi_masks(space, rois)
		sigma_map = nibabel.load(self.read(f'res-sigma_{space}-{design}')[0]).get_fdata()
		res_sigma = {roi: sigma_map[masks[roi]].reshape(-1) for roi in masks}

		return res_sigma

	### Validation and visualisation ###
	def check_logfile_headers(self):
		for sess in self.sessions:
			for task in self.tasks:
				logfiles = self.grab_raw(f'events_{task}')[sess]
				for n, logfile in enumerate(logfiles):
					with open(logfile, 'r') as f:
						print(f'{self.ID} sess {sess} {task} run {n+1} -> {next(f)}',end='')
			print('')
		print('')

	def check_analysis_rois(self):

		command = f'freeview -v {self.read("t1")[0]}'
		for roi_vol in self.read('rois_T1w'):
			command += f' -v {roi_vol}:opacity=0.5:colormap=heat'
		os.system(command)

	def check_localiser_contrast(self):

		command = f'freeview -v {self.read("t1")[0]}'
		for roi_vol in self.read('conimages_localiser_T1w'):
			command += f' -v {roi_vol}:opacity=0.5:colormap=heat'
		os.system(command)

	def list_available_confound_regressors(self):
		confounds_tsv = self.grab({'task':self.tasks[0], 'desc':'confounds', 'extension':'tsv'})
		print(list(confounds_tsv[0].get_df().columns.values))



class Godkoffer(Koffer):

	def __init__(self, configmod):
		super().__init__(configmod, 'group')
		self.subjects  = get_all_subs(configmod)
		self.configmod = configmod

	def add(self, files, key, subfold, descs, move=True):

		def copy_move_file(src, dst):
			if not os.path.exists(os.path.join(os.path.split(dst)[0])):
				os.makedirs(os.path.join(os.path.split(dst)[0]))
			if move:
				shutil.move(src, dst)
			else:
				shutil.copy(src, dst) 	

		def compute_fname(file_path, desc):
			extension = '.'.join(file_path.split('.')[1:])
			bidsname  = f'{self.ID}_desc-{desc}'
			fname = os.path.join(self.koffer_path, subfold, f'{bidsname}.{extension}')
			return fname

		if type(files) is str:
			files = [files]
		
		stored_files = []
		if type(descs) is str:
			descs = [descs] * len(files)
		for file, desc, in zip(files, descs):
			dst = compute_fname(file, desc)
			copy_move_file(file, dst)
			stored_files.append(dst)
		
		self.koffer[key] = stored_files
		self.update_koffer_json()

	def run_reconall(self):
		
		reconall = Tode(freesurfer.preprocess.ReconAll(), name='reconall')
		reconall.inputs.T1_files     = self.config.mni
		reconall.inputs.subjects_dir = self.config.fsSubjects
		reconall.inputs.subject_id   = 'mni'
		reconall.inputs.directive    = 'all'
		reconall.inputs.parallel     = True
		reconall.inputs.openmp       = 12
		reconall.run()
		reconall.cleanTempFiles()

	def store_reconall_atlas(self):

		fs_grabber = Tode(interface = io.FreeSurferSource(), name='fs_grabber')
		fs_grabber.inputs.subjects_dir = self.config.fsSubjects
		fs_grabber.inputs.subject_id   = 'mni'
		fs_grabber.run()

		fs_a2009s = [f for f in fs_grabber.result.outputs.aparc_aseg if 'a2009s' in f][0]
		nii_conv = Tode(freesurfer.MRIConvert(out_type = 'niigz'), name = 'nii_conv')
		nii_conv.inputs.out_file = 'atlas.nii.gz' # avoids extension errors (filename has .)
		nii_conv.inputs.in_file  = fs_a2009s
		nii_conv.run()

		nii_conv_t1 = Tode(freesurfer.MRIConvert(out_type = 'niigz'), name = 'nii_conv_T1')
		nii_conv_t1.inputs.in_file = fs_grabber.result.outputs.brain
		nii_conv_t1.run() 

		warp_transform = Tode(ants.Registration(), name = 'ants_warp_computation')
		warp_transform.inputs.transforms   = ['Rigid']
		warp_transform.inputs.metric = ['MI']
		warp_transform.inputs.smoothing_sigmas = [[8, 4, 2, 1, 0]]
		warp_transform.inputs.shrink_factors = [[16, 8, 4, 2, 1]]
		warp_transform.inputs.number_of_iterations = [[10000, 2000, 1000, 500, 200]]
		warp_transform.inputs.convergence_threshold = [1.e-6]
		warp_transform.inputs.transform_parameters = [(0.1,)]
		warp_transform.inputs.fixed_image  = self.config.mni
		warp_transform.inputs.moving_image = nii_conv_t1.result.outputs.out_file
		warp_transform.inputs.num_threads = 12
		warp_transform.inputs.verbose = True
		warp_transform.run()

		warper = Tode(ants.ApplyTransforms(), name = 'ants_warp_application')
		warper.inputs.input_image_type = 3
		warper.inputs.interpolation = 'NearestNeighbor'
		warper.inputs.invert_transform_flags = False 
		warper.inputs.num_threads = 12
		warper.inputs.reference_image = self.config.mni
		warper.inputs.input_image = nii_conv.result.outputs.out_file
		warper.inputs.transforms = warp_transform.result.outputs.forward_transforms
		warper.run()

		self.add(warper.result.outputs.output_image, 'atlas',  'anat', 'atlas')

		fs_grabber.cleanTempFiles()
		nii_conv.cleanTempFiles()
		warp_transform.cleanTempFiles()
		warper.cleanTempFiles()

	def create_analysis_rois(self, rois=None):

		base_dir = create_tmp_folder(mkdir=True)

		if rois is None:
			rois = list(self.config.rois.keys())
		elif type(rois) not in [list or tuple]:
			rois = [list]

		analysis_roi_paths, descs = [], []

		for key, roi_config in [(k, r) for k, r in self.config.rois.items() if k in rois]:

			# Nominal: directly take the template
			if roi_config['method'] == 'nominal':
				mask_path = roi_config['path']

			# Localiser: take the average contrast across subjects
			elif roi_config['method'] == 'localiser':

				# Load template masks in MNI
				prior_img   = nibabel.load(roi_config['prior']).get_fdata()
				mask_affine = nibabel.load(roi_config['prior']).affine

				contrast_img = np.zeros(prior_img.shape)
				for s in self.subjects:
					tk = Tabokoffer(self.configmod, s)
					subject_con_files = tk.read('conimages_localiser_MNI')
					subject_con_img   = sum([nibabel.load(f).get_fdata() for f in subject_con_files])
					contrast_img     += (prior_img > 0).astype(float) * subject_con_img

				# Thresholding
				threshold = 0
				for step in [1, 0.1, 0.01, 0.001, 0.0001, 0.00001]:
					while (contrast_img > threshold).sum() > roi_config['size']:
						threshold += step
					threshold -= step
				threshold += step

				mask_array = (contrast_img > threshold).astype(float)
				mask_path  = create_tmp_filename(base_dir, extension='nii.gz')
				nibabel.save(nibabel.Nifti1Image(mask_array, mask_affine), mask_path)
			
			elif roi_config['method'] == 'freesurfer':

				if 'atlas' not in self.koffer:
					self.run_reconall()
					self.store_reconall_atlas()

				atlas_array = nibabel.load(self.read('atlas')[0]).get_fdata()
				mask_affine = nibabel.load(self.read('atlas')[0]).affine
				mask_array  = (atlas_array==roi_config['code']).astype(float)
				mask_path   = create_tmp_filename(base_dir, extension='nii.gz')
				nibabel.save(nibabel.Nifti1Image(mask_array, mask_affine), mask_path)
				
			analysis_roi_paths.append(mask_path)
			descs.append(f'roi-{key}_space-MNI')

		self.add(analysis_roi_paths, 'rois_MNI', 'anat', descs, move=False)
		shutil.rmtree(base_dir)

	def extract_roi_masks(self, rois='all'):

		masks = dict()
		for mask_image, desc in zip(self.read(f'rois_MNI'), self.get_descs(f'rois_MNI')):
			key = [d.replace('roi-', '') for d in desc.split('_') if 'roi-' in d][0]
			if rois == 'all' or key in rois:
				this_mask = (nibabel.load(mask_image).get_fdata() == 1)
				if len(this_mask.shape) == 4:
					this_mask = this_mask[:, :, :, 0]
				masks[key] = this_mask

		return masks

	def extract_logevs(self, design, rois=None):

		if rois is not None:
			masks = self.extract_roi_masks(rois)

		files  = Tabokoffer(self.configmod, self.subjects[0]).read(f'{design}_logev_MNI')
		models = [f[f.find(f'-{design}-')+len(f'-{design}-'):f.find('_logev')] for f in files]

		if rois is None:
			logevs = dict([(m, []) for m in models])
		else:
			logevs = dict([(r, dict([(m, []) for m in models])) for r in masks])

		for s in self.subjects:
			files = Tabokoffer(self.configmod, s).read(f'{design}_logev_MNI')
			for m, f in zip(models, files):
				logev_array = nibabel.load(f).get_fdata()
				if rois is None:
					logevs[m].append(logev_array)
				else:
					for r in masks:
						logevs[r][m].append(logev_array[masks[r]].reshape(-1))

		return logevs


class Tode(engine.Node):
	
	def __init__(self, *args, **kwargs):
		super(Tode, self).__init__(*args, **kwargs)
		self.base_dir = create_tmp_folder()

	def cleanTempFiles(self):
		shutil.rmtree(self.base_dir)


class MapTode(engine.MapNode):
	
	def __init__(self,  *args, **kwargs):
		super(MapTode, self).__init__( *args, **kwargs)
		self.base_dir = create_tmp_folder()

	def cleanTempFiles(self):
		shutil.rmtree(self.base_dir)


def check_completed(configmod, printout=True):

	config = importlib.import_module(configmod).Config()
	layout = bids.layout.BIDSLayout(config.data)
	subs = layout.get_subjects()

	completed = dict([(subID, []) for subID in subs])

	for subID in subs:
		json_path = os.path.join(config.koffer, f'sub-{subID}', f'k{subID}.json')
		if os.path.isfile(json_path):
			with open(json_path, 'r') as f:
				koffer = json.load(f)
			if printout:
				print(f'Subject {subID}:')
			for key in koffer:
				completed[subID].append(key)
				if printout:
					print(f' |---- {key}')
		else:
			print(f'Subject {subID} does not have a koffer json')

	if not printout:
		return completed

def check_subjects_missing(configmod, key, printout=True):

	completed = check_completed(configmod, printout=False)
	missing = [subID for subID in completed.keys() if key not in completed[subID]]

	if printout:
		if missing == []:
			print(f'{key} exists in all subjects')
		else:
			print(f'{key} is missing in {", ".join(missing)}')

	if not printout:
		return missing

def get_all_subs(configmod, exclude=[]):

	config = importlib.import_module(configmod).Config()
	layout = bids.layout.BIDSLayout(config.data)
	subs = layout.get_subjects()

	if type(exclude) is not list:
		exclude = [exclude]

	for subID in exclude:
		if type(subID) is int:
			subID = f'{subID:02}'
		if subID in subs:
			subs.remove(subID)

	return subs

def filter_betas(spmmat, beta_images, bayesian=False):

	no_interest_regressors = ['Realign', 'Outlier', 'keypress', ')x']
	
	SPM         = scipy.io.loadmat(spmmat)
	spm_beta_ix = 15 if bayesian else 13               
	beta_names  = [r[5][0] for r in SPM['SPM'][0,0][spm_beta_ix][0]]

	filtered_betas, descriptions = [], []
	for d, b in zip(beta_names, beta_images):
		if not any([r in d for r in no_interest_regressors]):
			run_n_ix0 = d.find('Sn(') + 3
			run_n_ix1 = int(d[run_n_ix0:].find(')') + run_n_ix0)
			run_n = int(d[run_n_ix0:run_n_ix1])
			beta_name = f'{d[(run_n_ix1+2):]}_run-{run_n:02d}'
			for s0, s1 in [('*bf(1)', ''), ('*bf(2)', '-hrfder1'), ('*bf(3)', '-hrfder2'), ('^1', '')]:
				beta_name = beta_name.replace(s0, s1)
			filtered_betas.append(b)
			descriptions.append(beta_name)

	return filtered_betas, descriptions
	
def create_tmp_folder(mkdir=False):

	characters_to_use = string.ascii_uppercase + string.ascii_lowercase
	tmp_folder = ''.join(random.choice(characters_to_use) for _ in range(4))
	tmp_folder = os.path.join(WORKINGDIR, tmp_folder)

	if mkdir:
		os.mkdir(tmp_folder)

	return tmp_folder

def create_tmp_filename(base_dir='./', extension=None):

	characters_to_use = string.ascii_uppercase + string.ascii_lowercase
	tmp_filename = ''.join(random.choice(characters_to_use) for _ in range(4))
	
	if extension is not None:
		tmp_filename += f'.{extension}'

	return os.path.join(base_dir, tmp_filename)


