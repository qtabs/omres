import nitabs3
import sys
import time
import os
import pickle


# Parameters
configmod  = 'omresconfig'
subjects   = nitabs3.get_all_subs(configmod)
rois       = ['ICL', 'ICR', 'MGBL', 'MGBR', 'aHGL', 'aHGR']
conditions = ['std0', 'std1', 'std2', 'om4', 'om5', 'om6']
nuissance  = ['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'] + ['csf', 'white_matter']


# Storage initialisation
funcloc_t1w, zscores_t1w, betas_t1w, noise_t1w = dict(), dict(), dict(), dict()


# Pipeline
for s in subjects:
	tt = time.time()

	# Instantiate koffer
	tk = nitabs3.Tabokoffer(configmod, s)

	# Preprocessing
	# tk.check_logfile_headers()
	# tk.run_fmriprep()
	# tk.store_reconall_atlas()

	# Estimation
	# tk.estimate('localiser', nuissance)
	# tk.estimate_residuals_sigma('localiser', nuissance)
	# tk.create_analysis_rois()
	tk.estimate('omres', nuissance)

	# Extraction
	betas_t1w[s]   = tk.extract_betas('omres', 'T1w', conditions, rois)
	funcloc_t1w[s] = tk.extract_funcloc_contrasts('T1w', rois)
	zscores_t1w[s] = tk.extract_zscores('omres', 'T1w', conditions, rois)
	noise_t1w[s]   = tk.extract_residuals_sigma('localiser', 'T1w', rois)
	
	print(f'\n\n\nSubject{s} took {(time.time()-tt)/60:.1f} minutes\n\n')


# Saving results
if not os.path.exists('./results/'): os.makedirs('./results/')

for results in ['funcloc_t1w','zscores_t1w', 'betas_t1w', 'noise_t1w']:
	with open(f'./results/{results}.pickle', 'wb') as f:
		pickle.dump(globals()[results], f)

