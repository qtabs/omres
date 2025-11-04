import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import seaborn as sns
import glob
import pandas as pd
import pickle
import scipy
import importlib
import os


# PLoS Figures max 7.5 x 8.75 inchesm, fontsize 8
# Fontsizes (min PLoS 8)
plt.rc('font',   size = 8)      # controls default text sizes
plt.rc('axes',   titlesize = 8) # fontsize of the axes title
plt.rc('xtick',  labelsize = 8) # fontsize of the tick labels
plt.rc('ytick',  labelsize = 8) # fontsize of the tick labels
plt.rc('legend', fontsize  = 8) # legend fontsize
plt.rc('figure', titlesize = 9) # fontsize of the figure title

textoptions = {'fontsize': 7, 'horizontalalignment': 'right'}
colors      = {'L': 'tab:blue', 'R': 'tab:red', 'Left': 'tab:blue', 'Right': 'tab:red'}
alphas      = {1E10: 'ⁿˢ', 0.05: '*', 0.01: '**', 0.001: '***', 0.0001: '****'}
roinames    = {'IC': 'IC', 'MGB': 'MGB', 'aHG': 'PAC'}


def _format_pvalue_latex_(p):
	"""Format p-value in LaTeX scientific notation, capped at 1"""
	# Cap p-values at 1 (important for corrected p-values from Holm-Bonferroni)
	p = min(p, 1.0)

	if p >= 0.01:
		return f"${p:.2g}$"
	else:
		# Convert to scientific notation
		exp = int(np.floor(np.log10(p)))
		mantissa = p / (10 ** exp)
		if mantissa >= 9.95:  # Round up case
			mantissa = 1.0
			exp += 1
		return f"${mantissa:.1g} \\times 10^{{{exp}}}$"


def _merge_latex_tables_(df_d, df_p, df_p_corr, row_label_formatter=None, p_threshold=0.05):
	"""
	Merge three dataframes (Cohen's d, p-values, p-values corrected) into one LaTeX table.
	"""

	# Ensure all dataframes have the same index order
	df_d = df_d.sort_index()
	df_p = df_p.sort_index()
	df_p_corr = df_p_corr.sort_index()

	# Determine number of columns
	ncols = len(df_d.columns) + 1

	# Start building the LaTeX table
	col_format = 'l' + 'r' * (ncols - 1)
	lines = []
	lines.append(f"\\begin{{tabular}}{{{col_format}}}")
	lines.append("\\toprule")

	# Section 1: Cohen's d
	lines.append(f"\\multicolumn{{{ncols}}}{{l}}{{Cohen's d}} \\\\")

	# Column headers with proper spacing
	col_names = [str(c).replace('_', '\\_') for c in df_d.columns]
	header = f"{'':15} & {' & '.join([f'{c:20}' for c in col_names])} \\\\"
	lines.append(header)
	lines.append("\\midrule")

	# Data rows with bold formatting based on significance of corrected p-values
	for idx in df_d.index:
		row_label = row_label_formatter(idx) if row_label_formatter else str(idx)
		values = []
		for col in df_d.columns:
			val = df_d.loc[idx, col]
			p_corr_val = df_p_corr.loc[idx, col]
			# Bold if corrected p-value is significant
			if p_corr_val < p_threshold:
				values.append(f"\\textbf{{{val:.3g}}}".ljust(20))
			else:
				values.append(f"{val:.3g}".ljust(20))
		line = f"{row_label:15} & {' & '.join(values)} \\\\"
		lines.append(line)

	lines.append("\\midrule")
	lines.append("\\\\")

	# Section 2: p-values (corrected)
	lines.append(f"\\multicolumn{{{ncols}}}{{l}}{{$p$-values (corrected)}} \\\\")
	lines.append(header)
	lines.append("\\midrule")

	for idx in df_p_corr.index:
		row_label = row_label_formatter(idx) if row_label_formatter else str(idx)
		values = []
		for col in df_p_corr.columns:
			val = df_p_corr.loc[idx, col]
			values.append(_format_pvalue_latex_(val).ljust(20))
		line = f"{row_label:15} & {' & '.join(values)} \\\\"
		lines.append(line)

	lines.append("\\midrule")
	lines.append("\\\\")

	# Section 3: p-values (uncorrected)
	lines.append(f"\\multicolumn{{{ncols}}}{{l}}{{$p$-values (uncorrected)}} \\\\")
	lines.append(header)
	lines.append("\\midrule")

	for idx in df_p.index:
		row_label = row_label_formatter(idx) if row_label_formatter else str(idx)
		values = []
		for col in df_p.columns:
			val = df_p.loc[idx, col]
			values.append(_format_pvalue_latex_(val).ljust(20))
		line = f"{row_label:15} & {' & '.join(values)} \\\\"
		lines.append(line)

	lines.append("\\bottomrule")
	lines.append("\\end{tabular}")

	return "\n".join(lines) + "\n\n"


def fig2_omres(z, subjects, rois, conditions):

	levels  = list(dict.fromkeys([r[:-1] for r in rois]))
	options_violin = {'linewidth': 0.6, 'inner': 'box', 'width': 1, 'inner_kws': {'color': 'k'}}
	options_poinpl = {'errorbar': 'ci', 'markersize': 2, 'linewidth': 1.5, 'color': 'k'}

	contrasts = {'adaptation (std0-std1)': ('std0','std1'), 'omissions (om4-std1)': ('om4','std1')}

	p = {roi: {contrast_name: np.nan for contrast_name in contrasts} for roi in rois}
	d = {roi: {contrast_name: np.nan for contrast_name in contrasts} for roi in rois}

	fig, axs = plt.subplots(1, len(contrasts))

	for ax, contrast_name in zip(axs, contrasts):

		contrast_dict = {'Level':[],'Hemisphere':[],'Subject':[],contrast_name:[],'Val0':[],'Val1':[]}
		beta0, beta1  = contrasts[contrast_name]

		for roi in rois:
			for sub in subjects:
				if z[sub][roi]['std0'] != []:
					run_vals = [v0-v1 for v0, v1 in zip(z[sub][roi][beta0], z[sub][roi][beta1])]
					contrast = np.mean([run_val.mean() for run_val in run_vals])
					val0 = np.mean([v.mean() for v in z[sub][roi][beta0]])
					val1 = np.mean([v.mean() for v in z[sub][roi][beta1]])
				else:
					contrast = np.nan
					val0 = np.nan
					val1 = np.nan
				contrast_dict['Subject']     += [sub]
				contrast_dict['Level']       += [roinames[roi[:-1]]]
				contrast_dict['Hemisphere']  += [roi[-1]]
				contrast_dict[contrast_name] += [contrast]
				contrast_dict['Val0']        += [val0]
				contrast_dict['Val1']        += [val1]

		df = pd.DataFrame(contrast_dict)

		sns.violinplot(df, x='Level', y=contrast_name, hue='Hemisphere', ax=ax, split=True,
					  palette=colors, legend=('std0' in contrast_name), **options_violin)

		ax.hlines(0, -0.5, len(levels) - 0.5, color='k', linewidth=0.6)
		ax.set_xlim(-0.5, len(levels) - 0.5)
		ax.set_ylim(-0.75, 2.0)
		ax.set_ylabel(contrast_name)
		ax.set_xlabel('')

		for roi in rois:
			roisub = (df.Level==roinames[roi[:-1]]) & (df.Hemisphere==roi[-1])
			x = df.Val0[roisub & df.Val0.notna()]
			y = df.Val1[roisub & df.Val1.notna()]
			d[roi][contrast_name] = (np.mean(x) - np.mean(y)) / ((np.std(x)**2 + np.std(y)**2)/2)**.5
			p[roi][contrast_name] = scipy.stats.ranksums(x, y, 'greater').pvalue 


	# Correcting p-values and marking significance using Holm-Bonferroni
	p_corr       = {roi: {contrast_name: np.nan for contrast_name in contrasts} for roi in rois}
	significance = {roi: {contrast_name: '' for contrast_name in contrasts} for roi in rois}
	group_alphas = sorted(alphas.keys())

	all_keys  = [(roi, contrast_name) for roi in rois for contrast_name in contrasts]
	all_pvals = np.array([p[roi][contrast_name] for roi in rois for contrast_name in contrasts])

	all_pvals_corrected = ((-all_pvals).argsort().argsort() + 1) * all_pvals 

	for k in np.argsort(all_pvals):
		roi, cn = all_keys[k]
		sign_alpha_ix = np.argwhere([all_pvals_corrected[k] < a for a in group_alphas]).min()
		significance[roi][cn] = alphas[group_alphas[sign_alpha_ix]]
		significance[roi][cn] = '-' if significance[roi][cn] == 'ⁿˢ' else significance[roi][cn]
		p_corr[roi][cn] = all_pvals_corrected[k]
		group_alphas = group_alphas[sign_alpha_ix:]

	# Annotating of significance
	sign_options = {'fontsize': 6, 'horizontalalignment': 'center' ,'verticalalignment': 'center'}

	x_eps = 0.15
	sigmark_ys = {'adaptation': [1.4,  1.4, 1.6], 'omissions': [1.4, 1.8, 1.8]}
	sigmark_xs = {roi: n//2 + (-1)**(n%2+1) * x_eps for n, roi in enumerate(rois)}

	for ax, contrast_name in zip(axs, contrasts):
		for n, roi in enumerate(rois):
			x_text, y_text = sigmark_xs[roi], sigmark_ys[contrast_name.split(' ')[0]][n//2]
			ax.text(x_text, y_text, significance[roi][contrast_name], **sign_options)

	# Titles
	for ax, contrast_name in zip(axs, contrasts):
		ax.set_title(contrast_name.split(' ')[0])

	# Panel captions
	for ax, caption in zip(axs, ['A', 'B']):
		ax.annotate(caption, xy=(-0.18, 1.05),  xycoords="axes fraction", fontsize=10)

	plt.subplots_adjust(left=0.08, bottom=0.10, right=0.99, top=0.90, wspace=0.28)
	fig.set_size_inches(7.5, 2.0)
	# plt.savefig(f'./figures/fig2-omres.svg')
	plt.savefig(f'./figures/fig2-omres.pdf')

	df_d      = pd.DataFrame(d)
	df_p      = pd.DataFrame(p)
	df_p_corr = pd.DataFrame(p_corr)

	# Row label formatter: extract conditions from contrast name
	def format_row_label(contrast_name):
		# Extract conditions from "adaptation (std0-std1)" or "omissions (om4-std1)"
		import re
		match = re.search(r'\(([^-]+)-([^)]+)\)', contrast_name)
		if match:
			c1, c2 = match.groups()
			return f"{c1}  $>$ {c2}"
		return contrast_name

	# Merge and write single LaTeX file (don't transpose - original structure is correct)
	merged_latex = _merge_latex_tables_(df_d, df_p, df_p_corr,
									  row_label_formatter=format_row_label)
	with open('./figures/tabS1-omres.tex', 'w') as f:
		f.write(merged_latex)


def fig3_rois(z, subjects, rois):

	conditions     = ['om4', 'om5', 'om6']
	levels         = list(dict.fromkeys([r[:-1] for r in rois]))
	options_violin = {'linewidth': 0.6, 'inner': 'box', 'width': 1, 'inner_kws': {'color': 'k'}}
	options_poinpl = {'errorbar': 'ci', 'markersize': 2, 'linewidth': 1.5, 'color': 'k'}

	contrasts = [("om4", "om5"), ("om4", "om6"), ("om5", "om6")]
	p  = dict([(roinames[level], [np.nan for contrast in contrasts]) for level in levels])
	d  = dict([(roinames[level], [np.nan for contrast in contrasts]) for level in levels])


	z_dict = {'Subject':[], 'Level':[], 'Hemisphere':[], 'Condition':[], 'z-score': []}
	for roi in rois:
		for sub in subjects:
			if z[sub][roi]['std0'] != []:
				for condition in conditions:
					avg_z = np.mean([run.mean() for run in z[sub][roi][condition]])
					z_dict['Subject']    += [sub]
					z_dict['Level']      += [roinames[roi[:-1]]]
					z_dict['Hemisphere'] += ['Left' if roi[-1]=='L' else 'Right']
					z_dict['Condition']  += [condition]
					z_dict['z-score']    += [avg_z]

	df = pd.DataFrame(z_dict)

	fig, axs = plt.subplots(1, len(levels))

	for level, ax in zip(levels, axs):
		subset = df[df.Level==roinames[level]]
		legend = 'full' if level=='IC' else False
		sns.violinplot(subset, y='z-score', x='Condition', hue='Hemisphere', 
					   palette=colors, split=True, ax=ax, legend=legend, **options_violin)
		sns.pointplot(subset, y='z-score', x='Condition', ax=ax, **options_poinpl)
	
		ax.set_ylim([-1.35, 1.75])
		ax.set_title(roinames[level])
		ax.set_xticks(range(len(conditions)), labels=conditions)
		ax.set_xlabel('')

		if level == 'IC':
			ax.set_yticks(np.arange(-1, 1.6, 0.5), labels=['-1', '-.5', '0', '+.5', '+1', ''])
			ax.set_ylabel(f'z-score')
		else:
			ax.set_yticks(np.arange(-1, 1.6, 0.5), labels='')
			ax.set_ylabel('')

		for n, cont in enumerate(contrasts):
			x = [np.array(subset[subset['Condition']==c].groupby('Subject').mean('Hemisphere')) for c in cont]
			x = [np.array([val for val in x_ix.flatten() if not np.isnan(val)]) for x_ix in x]
			d[roinames[level]][n] = (np.mean(x[0]) - np.mean(x[1])) / ((np.std(x[0])**2 + np.std(x[1])**2)/2)**.5
			p[roinames[level]][n] = scipy.stats.ranksums(x[0], x[1], 'greater').pvalue 

	# Correcting p-values and marking significance using Holm-Bonferroni
	p_corr       = dict([(roinames[level], [np.nan for c in contrasts]) for level in levels])
	significance = dict([(roinames[level], ['' for c in contrasts]) for level in levels])
	group_alphas = sorted(alphas.keys())

	all_keys  = [(roinames[level], n) for level in levels for n in range(len(contrasts))]
	all_pvals = np.array([p[roinames[level]][n] for level in levels for n in range(len(contrasts))])

	all_pvals_corrected = ((-all_pvals).argsort().argsort() + 1) * all_pvals 

	for k in np.argsort(all_pvals):
		level, n = all_keys[k]
		sign_alpha_ix = np.argwhere([all_pvals_corrected[k] < a for a in group_alphas]).min()
		significance[level][n] = alphas[group_alphas[sign_alpha_ix]]
		significance[level][n] = '-' if significance[level][n] == 'ⁿˢ' else significance[level][n]
		p_corr[level][n] = all_pvals_corrected[k]
		group_alphas = group_alphas[sign_alpha_ix:]

	# Annotating significance
	sign_options = {'fontsize': 6, 'horizontalalignment': 'center' ,'verticalalignment': 'center'}

	sigbar_ys = {'IC':  [-0.95, -1.1, -0.8],
				 'MGB': [-0.6, 1.4, 0.9],
				 'PAC': [-0.3, 1.5, 0.7]}

	for level, ax in zip(levels, axs):
		for n in range(len(contrasts)):
			x0, x1 = conditions.index(contrasts[n][0]), conditions.index(contrasts[n][1])
			y0 = sigbar_ys[roinames[level]][n]
			y1 = sigbar_ys[roinames[level]][n] + .05*np.sign(sigbar_ys[roinames[level]][n])
			ax.plot([x0, x0], [y0, y1], 'k-', linewidth=0.7)
			ax.plot([x0, x1], [y1, y1], 'k-', linewidth=0.7)
			ax.plot([x1, x1], [y1, y0], 'k-', linewidth=0.7)
			y_text = y1 + .05 if y1>0 else y1 - 0.15
			y_text = y_text + 0.07 if significance[roinames[level]][n]=='-' and y_text < 0 else y_text
			ax.text(.5*(x0+x1), y_text, significance[roinames[level]][n], **sign_options)

	axs[0].legend(loc='upper right')

	plt.subplots_adjust(left=0.08, bottom=0.15, right=0.99, top=0.88, wspace=0.08)
	fig.set_size_inches(7.5, 2.0)
	# plt.savefig(f'./figures/fig3-rois.svg')
	plt.savefig(f'./figures/fig3-rois.pdf')

	contrast_names = [f'{c[0]} > {c[1]}' for c in contrasts]
	df_d      = pd.DataFrame(d, index=contrast_names)
	df_p      = pd.DataFrame(p, index=contrast_names)
	df_p_corr = pd.DataFrame(p_corr, index=contrast_names)

	# Row label formatter: convert ">" to "$>$"
	def format_row_label(label):
		return label.replace(' > ', '  $>$ ')

	# Merge and write single LaTeX file
	merged_latex = _merge_latex_tables_(df_d, df_p, df_p_corr,
									  row_label_formatter=format_row_label)
	with open('./figures/tabS2-negpe.tex', 'w') as f:
		f.write(merged_latex)


def fig4_emergence(z, funcloc, subjects, rois, conditions):

	levels  = list(dict.fromkeys([r[:-1] for r in rois]))
	options_violin = {'linewidth': 0.6, 'inner': 'box', 'width': 1, 'inner_kws': {'color': 'k'}}
	options_poinpl = {'errorbar': 'ci', 'markersize': 2, 'linewidth': 1.5, 'color': 'k'}

	contrasts = {'adaptation std0-std1':     ('std0', 'std1'), 
				 'omissions om4-std1':       ('om4',  'std1'), 
				 'sPE om4-om6':              ('om4',  'om6'),
				 'BOLD-sensitivity sound-silence_(localiser)': None}

	p = {cn: dict() for cn in contrasts if 'BOLD' not in cn}
	d = {cn: dict() for cn in contrasts if 'BOLD' not in cn}

	fig, axs = plt.subplots(1, len(contrasts))

	for ax, contrast_name in zip(axs, contrasts):

		contrast_dict = {'Level':[], 'Hemisphere':[], 'Subject':[], contrast_name:[]}

		for roi in rois:
			for sub in subjects:
				if z[sub][roi]['std0'] != []:
					if contrasts[contrast_name] is None:
						contrast = np.mean(funcloc[sub][roi])
					else:
						beta0, beta1 = contrasts[contrast_name]
						run_vals = [v0-v1 for v0, v1 in zip(z[sub][roi][beta0], z[sub][roi][beta1])]
						contrast = np.mean([(val/funcloc[sub][roi]).mean() for val in run_vals])
				else:
					contrast = np.nan
				contrast_dict['Subject']     += [sub]
				contrast_dict['Level']       += [roinames[roi[:-1]]]
				contrast_dict['Hemisphere']  += [roi[-1]]
				contrast_dict[contrast_name] += [contrast]

		df = pd.DataFrame(contrast_dict)

		sns.violinplot(df, x='Level', y=contrast_name, hue='Hemisphere', ax=ax, split=True,
					  palette=colors, legend=('BOLD' in contrast_name), **options_violin)

		ax.hlines(0, -0.5, len(levels) - 0.5, color='k', linewidth=0.6)
		ax.set_xlim(-0.5, len(levels) - 0.5)
		if 'BOLD' not in contrast_name:
			ax.set_ylim(-0.75, 1.3)
		title, ylabel = contrast_name.split(' ')
		ax.set_title(title)
		ax.set_ylabel(ylabel.replace('-', ' - ').replace('_', ' '))
		ax.set_xlabel('')

		if 'BOLD' not in contrast_name:
			for l1, l2 in zip(levels[:-1], levels[1:]):
				x = df[contrast_name][(df.Level==roinames[l1])][df[contrast_name].notna()]
				y = df[contrast_name][(df.Level==roinames[l2])][df[contrast_name].notna()]
				d[contrast_name][f'{l1}-{l2}'] = (x.mean()-y.mean())/((x.std()**2+y.std()**2)/2)**.5
				p[contrast_name][f'{l1}-{l2}'] = scipy.stats.ranksums(x, y, 'two-sided').pvalue 

	p_corr       = {cn: {comp: []} for cn in contrasts if 'BOLD' not in cn for comp in d[cn]}
	significance = {cn: {comp: []} for cn in contrasts if 'BOLD' not in cn for comp in d[cn]}

	# Correcting p-values and marking significance using Holm-Bonferroni
	group_alphas = sorted(alphas.keys())

	all_keys  = [(cn, comp) for cn in d for comp in d[cn]]
	all_pvals = np.array([p[cn][comp] for cn in d for comp in d[cn]])
	all_pvals_corrected = ((-all_pvals).argsort().argsort() + 1) * all_pvals 

	for k in np.argsort(all_pvals):
		cn, comp = all_keys[k]
		sign_alpha_ix = np.argwhere([all_pvals_corrected[k] < a for a in group_alphas]).min()
		significance[cn][comp] = alphas[group_alphas[sign_alpha_ix]]
		if significance[cn][comp] == 'ⁿˢ':
			significance[cn][comp] = '-'
		p_corr[cn][comp] = all_pvals_corrected[k]
		group_alphas = group_alphas[sign_alpha_ix:]


	# Annotations of significance
	sign_options = {'fontsize': 6, 'horizontalalignment': 'center' ,'verticalalignment': 'center'}

	x_eps = 0.15
	sigmark_ys = {'adaptation': [1.1, 0.5], 'omissions': [0.7, 0.6], 'sPE': [0.8, 0.7]}
	sigmark_xs = [(0, 1), (1, 2)]

	for ax, contrast_name in zip(axs, contrasts):
		if 'BOLD' not in contrast_name:
			for n, rois in enumerate(d[contrast_name]):
				x_text, y_text = n + 0.5, sigmark_ys[contrast_name.split(' ')[0]][n] + 0.1
				x0, x1 = n, n + 1
				y0, y1 = y_text - 0.1, y_text - 0.05
				ax.plot([x0, x0], [y0, y1], 'k-', linewidth=0.7)
				ax.plot([x0, x1], [y1, y1], 'k-', linewidth=0.7)
				ax.plot([x1, x1], [y1, y0], 'k-', linewidth=0.7)
				ax.text(x_text, y_text, significance[contrast_name][rois], **sign_options)

	# Panel captions
	for ax, caption in zip(axs, ['A', 'B', 'C', 'D']):
		ax.annotate(caption, xy=(-0.18, 1.05),  xycoords="axes fraction", fontsize=10)

	plt.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.88, wspace=0.55)
	fig.set_size_inches(7.5, 2)
	# plt.savefig(f'./figures/fig4-emergence.svg')
	plt.savefig(f'./figures/fig4-emergence.pdf')

	df_d      = pd.DataFrame(d)
	df_p      = pd.DataFrame(p)
	df_p_corr = pd.DataFrame(p_corr)

	# Rename columns to simplified names
	col_rename = {
		'adaptation (std0-std1)': 'adaptation',
		'omissions (om4-std1)': 'omissions',
		'sPE (om4-om6)': 'negative PE'
	}
	df_d = df_d.rename(columns=col_rename)
	df_p = df_p.rename(columns=col_rename)
	df_p_corr = df_p_corr.rename(columns=col_rename)

	# Row label formatter: convert "IC-MGB" to "IC $\neq$ MGB"
	def format_row_label(label):
		return label.replace('-', ' $\\neq$ ')

	# Merge and write single LaTeX file
	merged_latex = _merge_latex_tables_(df_d, df_p, df_p_corr,
									  row_label_formatter=format_row_label)
	with open('./figures/tabS3-emergence.tex', 'w') as f:
		f.write(merged_latex)


def figS2_allregs(z, subjects, rois, conditions):

	levels  = list(dict.fromkeys([r[:-1] for r in rois]))
	options_violin = {'linewidth': 0.6, 'inner': None, 'width': 1}
	options_poinpl = {'errorbar': 'se', 'markersize': 2, 'linewidth': 1.5, 'color': 'k'}

	contrasts = [("om4", "om5"), ("om4", "om6"), ("om5", "om6")]
	p  = dict([(roinames[level], [np.nan for contrast in contrasts]) for level in levels])
	d  = dict([(roinames[level], [np.nan for contrast in contrasts]) for level in levels])


	z_dict = {'Subject':[], 'Level':[], 'Hemisphere':[], 'Condition':[], 'z-score': []}
	for roi in rois:
		for sub in subjects:
			if z[sub][roi]['std0'] != []:
				for condition in conditions:
					avg_z = np.mean([run.mean() for run in z[sub][roi][condition]])
					z_dict['Subject']    += [sub]
					z_dict['Level']      += [roinames[roi[:-1]]]
					z_dict['Hemisphere'] += ['Left' if roi[-1]=='L' else 'Right']
					z_dict['Condition']  += [condition]
					z_dict['z-score']    += [avg_z]

	df = pd.DataFrame(z_dict)

	fig, axs = plt.subplots(1, len(levels))

	for level, ax in zip(levels, axs):
		subset = df[df.Level==roinames[level]]
		legend = 'full' if level=='IC' else False
		sns.violinplot(subset, y='z-score', x='Condition', hue='Hemisphere', 
					   palette=colors, split=True, ax=ax, legend=legend, **options_violin)
		sns.pointplot(subset, y='z-score', x='Condition', ax=ax, **options_poinpl)
	
		ax.set_ylim([-1.35, 1.75])
		ax.set_title(roinames[level])
		ax.set_xticks(range(len(conditions)), labels=conditions)
		ax.set_xlabel('')

		if level == 'IC':
			ax.set_yticks(np.arange(-1, 1.6, 0.5), labels=['-1', '-.5', '0', '+.5', '+1', ''])
			ax.set_ylabel(f'z-score')
		else:
			ax.set_yticks(np.arange(-1, 1.6, 0.5), labels='')
			ax.set_ylabel('')

		for n, cont in enumerate(contrasts):
			x = [np.array(subset[subset['Condition']==c].groupby('Subject').mean('Hemisphere')) for c in cont]
			x = [np.array([val for val in x_ix.flatten() if not np.isnan(val)]) for x_ix in x]
			d[roinames[level]][n] = (np.mean(x[0]) - np.mean(x[1])) / ((np.std(x[0])**2 + np.std(x[1])**2)/2)**.5
			p[roinames[level]][n] = scipy.stats.ranksums(x[0], x[1], 'greater').pvalue 

	# Correcting p-values and marking significance using Holm-Bonferroni
	p_corr       = dict([(roinames[level], [np.nan for c in contrasts]) for level in levels])
	significance = dict([(roinames[level], ['' for c in contrasts]) for level in levels])
	group_alphas = sorted(alphas.keys())

	all_keys  = [(roinames[level], n) for level in levels for n in range(len(contrasts))]
	all_pvals = np.array([p[roinames[level]][n] for level in levels for n in range(len(contrasts))])

	all_pvals_corrected = ((-all_pvals).argsort().argsort() + 1) * all_pvals 

	for k in np.argsort(all_pvals):
		level, n = all_keys[k]
		sign_alpha_ix = np.argwhere([all_pvals_corrected[k] < a for a in group_alphas]).min()
		significance[level][n] = alphas[group_alphas[sign_alpha_ix]]
		significance[level][n] = '-' if significance[level][n] == 'ⁿˢ' else significance[level][n]
		p_corr[level][n] = all_pvals_corrected[k]
		group_alphas = group_alphas[sign_alpha_ix:]

	# Annotating significance
	sign_options = {'fontsize': 6, 'horizontalalignment': 'center' ,'verticalalignment': 'center'}

	sigbar_ys = {'IC':  [ 1.15, 0.95, 0.7],
				 'MGB': [-0.6, 1.4, 0.9],
				 'PAC': [-0.3, 1.5, 0.7]}

	for level, ax in zip(levels, axs):
		for n in range(len(contrasts)):
			x0, x1 = conditions.index(contrasts[n][0]), conditions.index(contrasts[n][1])
			y0 = sigbar_ys[roinames[level]][n]
			y1 = sigbar_ys[roinames[level]][n] + .05*np.sign(sigbar_ys[roinames[level]][n])
			ax.plot([x0, x0], [y0, y1], 'k-', linewidth=0.7)
			ax.plot([x0, x1], [y1, y1], 'k-', linewidth=0.7)
			ax.plot([x1, x1], [y1, y0], 'k-', linewidth=0.7)
			y_text = y1 + .05 if y1>0 else y1 - 0.15
			y_text = y_text + 0.07 if significance[roinames[level]][n]=='-' and y_text < 0 else y_text
			ax.text(.5*(x0+x1), y_text, significance[roinames[level]][n], **sign_options)

	axs[0].legend(loc='lower right')

	# Panel captions
	for ax, caption in zip(axs, ['A', 'B', 'C']):
		ax.annotate(caption, xy=(-0.18, 1.05),  xycoords="axes fraction", fontsize=10)

	plt.subplots_adjust(left=0.07, bottom=0.15, right=0.99, top=0.88, wspace=0.08)
	fig.set_size_inches(7.5, 1.8)
	# plt.savefig(f'./figures/figS2-allregs.svg')
	plt.savefig(f'./figures/figS2-allregs.pdf')


def figS3_participants(z, subjects, rois, conditions):

	levels    = list(dict.fromkeys([r[:-1] for r in rois]))
	options   = {'errorbar': 'se', 'markersize': 2, 'linewidth': 1.5}
	contrasts = {'omres': ('om4','std1'), 'sPE': ('om4', 'om6')}

	fig, axs = plt.subplots(len(subjects), len(levels))
	for sub_n, sub in enumerate(subjects):

		avg_z_roi = dict([(roi, dict()) for roi in rois])
		for roi in rois:
			for condition in conditions:
				avg_z_roi[roi][condition] = [run.mean() for run in z[sub][roi][condition]]
				if avg_z_roi[roi][condition] == []:
					avg_z_roi[roi][condition] = np.nan;

		# Compuing p-values
		p = dict([(con, dict()) for con in contrasts])
		for level in levels:
			for con in contrasts:
				x_cond, y_cond = contrasts[con]
				x = avg_z_roi[level+'L'][x_cond]+avg_z_roi[level+'R'][x_cond]
				y = avg_z_roi[level+'L'][y_cond]+avg_z_roi[level+'R'][y_cond]
				p[con][level] = scipy.stats.ranksums(x, y, 'greater').pvalue
		
		# Correcting p-values and marking significance using Holm-Bonferroni
		p_corr       = dict([(c, dict([(l, np.nan) for l in levels])) for c in contrasts])
		significance = dict([(c, dict([(l, '') for l in levels])) for c in contrasts])
		sub_alphas   = sorted(alphas.keys())

		all_keys  = [(c, l) for c in contrasts for l in levels if not np.isnan(p[c][l])]
		all_pvals = np.array([p[c][l] for c in contrasts for l in levels if not np.isnan(p[c][l])])

		all_pvals_corrected = ((-all_pvals).argsort().argsort() + 1)  * all_pvals 

		for n in np.argsort(all_pvals):
			cond, level = all_keys[n]
			sign_alpha_ix = np.argwhere([all_pvals_corrected[n] < a for a in sub_alphas]).min()
			significance[cond][level] = alphas[sub_alphas[sign_alpha_ix]]
			p_corr[cond][level] = all_pvals_corrected[n]
			sub_alphas = sub_alphas[sign_alpha_ix:]

		# Plotting
		for level_n, level in enumerate(levels):
			for k in ['L', 'R']:

				# Plotting errorbars
				data  = avg_z_roi[level+k]
				label = 'Left' if k == 'L' else 'Right'
				sns.pointplot(data, color=colors[k],  ax=axs[sub_n, level_n], label=label, **options)
				
				# Writing down significance of contrasts
				#pvs = f"p_omres = {p['omres'][level]:.2g} ({significance['omres'][level]})\n"
				#pvs += f"p_sPE  = {p['sPE'][level]:.2g} ({significance['sPE'][level]})"
				pvs = f"om{significance['omres'][level]}, pe{significance['sPE'][level]}"
				if np.isnan(p_corr['sPE'][level]):
					 pvs = ''
				plt.text(0.97, 0.75, pvs, transform=axs[sub_n, level_n].transAxes, **textoptions)


		# Aaxes labels and ticks
		for level_n, level in enumerate(levels):
			axs[sub_n, level_n].set_ylim([-1, 1.3])
			if level_n == 0:
				axs[sub_n, level_n].set_yticks([-1, 0, 1.0], labels = ['-1','0','+1'])
			else:
				axs[sub_n, level_n].set_yticks([-1, 0, 1.0], labels=[])

		axs[sub_n, 0].set_ylabel(f'Sub {sub}\n\nz-score')
		
		if sub_n == 0:
			for level_n, level in enumerate(levels):
				axs[sub_n, level_n].set_title(roinames[level])

		if sub_n != len(subjects)-1:
			for level_n, level in enumerate(levels):
				axs[sub_n, level_n].set_xticks(range(len(conditions)), labels=[])
		else:
			for level_n, level in enumerate(levels):
				axs[sub_n, level_n].set_xticks(range(len(conditions)), labels=conditions)


	axs[-1,1].legend(bbox_to_anchor=(-0.05, 1.05))

	plt.subplots_adjust(left=0.1, bottom=0.03, right=0.98, top=0.97, wspace=0.05, hspace=0.2)
	fig.set_size_inches(7.5, 8.0)
	# plt.savefig(f'./figures/figS3-participants.svg')
	plt.savefig(f'./figures/figS3-participants.pdf')


def figS4_residuals(z, noise, subjects, rois, conditions):

	levels  = list(dict.fromkeys([r[:-1] for r in rois]))
	options_violin = {'linewidth': 0.6, 'inner': 'box', 'width': 1, 'inner_kws': {'color': 'k'}}
	options_poinpl = {'errorbar': 'ci', 'markersize': 2, 'linewidth': 1.5, 'color': 'k'}

	contrasts = {'adaptation (std0-std1)':     ('std0', 'std1'), 
				 'omissions (om4-std1)':       ('om4',  'std1'), 
				 'sPE (om4-om6)':              ('om4',  'om6'),
				 'residuals noise':    None}

	p = {cn: dict() for cn in contrasts if 'residual' not in cn}
	d = {cn: dict() for cn in contrasts if 'residual' not in cn}

	fig, axs = plt.subplots(1, len(contrasts))

	for ax, contrast_name in zip(axs, contrasts):

		contrast_dict = {'Level':[], 'Hemisphere':[], 'Subject':[], contrast_name:[]}

		for roi in rois:
			for sub in subjects:
				if z[sub][roi]['std0'] != []:
					if contrasts[contrast_name] is None:
						contrast = np.mean(noise[sub][roi])
					else:
						beta0, beta1 = contrasts[contrast_name]
						run_vals = [v0-v1 for v0, v1 in zip(z[sub][roi][beta0], z[sub][roi][beta1])]
						contrast = np.mean([(val * noise[sub][roi]).mean() for val in run_vals])
				else:
					contrast = np.nan
				contrast_dict['Subject']     += [sub]
				contrast_dict['Level']       += [roinames[roi[:-1]]]
				contrast_dict['Hemisphere']  += [roi[-1]]
				contrast_dict[contrast_name] += [contrast]

		df = pd.DataFrame(contrast_dict)

		sns.violinplot(df, x='Level', y=contrast_name, hue='Hemisphere', ax=ax, split=True,
					  palette=colors, legend=('residual' in contrast_name), **options_violin)

		ax.hlines(0, -0.5, len(levels) - 0.5, color='k', linewidth=0.6)
		ax.set_xlim(-0.5, len(levels) - 0.5)
		ax.set_ylim(-10, 14)
		title, ylabel = contrast_name.split(' ')
		ax.set_title(title)
		ax.set_ylabel(ylabel.replace('(', '').replace(')', '').replace('-', ' - '))
		ax.set_xlabel('')

		if 'residual' not in contrast_name:
			for l1, l2 in zip(levels[:-1], levels[1:]):
				x = df[contrast_name][(df.Level==roinames[l1])][df[contrast_name].notna()]
				y = df[contrast_name][(df.Level==roinames[l2])][df[contrast_name].notna()]
				d[contrast_name][f'{l1}-{l2}'] = (x.mean()-y.mean())/((x.std()**2+y.std()**2)/2)**.5
				p[contrast_name][f'{l1}-{l2}'] = scipy.stats.ranksums(x, y, 'two-sided').pvalue 

	p_corr       = {cn: {comp: []} for cn in contrasts if 'residual' not in cn for comp in d[cn]}
	significance = {cn: {comp: []} for cn in contrasts if 'residual' not in cn for comp in d[cn]}

	# Correcting p-values and marking significance using Holm-Bonferroni
	group_alphas = sorted(alphas.keys())

	all_keys  = [(cn, comp) for cn in d for comp in d[cn]]
	all_pvals = np.array([p[cn][comp] for cn in d for comp in d[cn]])
	all_pvals_corrected = ((-all_pvals).argsort().argsort() + 1) * all_pvals 

	for k in np.argsort(all_pvals):
		cn, comp = all_keys[k]
		sign_alpha_ix = np.argwhere([all_pvals_corrected[k] < a for a in group_alphas]).min()
		significance[cn][comp] = alphas[group_alphas[sign_alpha_ix]]
		if significance[cn][comp] == 'ⁿˢ':
			significance[cn][comp] = '-'
		p_corr[cn][comp] = all_pvals_corrected[k]
		group_alphas = group_alphas[sign_alpha_ix:]


	# Annotations of significance
	sign_options = {'fontsize': 6, 'horizontalalignment': 'center' ,'verticalalignment': 'center'}

	x_eps = 0.15
	sigmark_ys = {'adaptation': [11.5, 7], 'omissions': [11, 12], 'sPE': [11.5, 12]}
	sigmark_xs = [(0, 1), (1, 2)]

	for ax, contrast_name in zip(axs, contrasts):
		if 'residual' not in contrast_name:
			for n, rois in enumerate(d[contrast_name]):
				x_text, y_text = n + 0.5, sigmark_ys[contrast_name.split(' ')[0]][n] + 1.5
				x0, x1 = n, n + 1
				y0, y1 = y_text - 0.5, y_text - 0.3
				ax.plot([x0, x0], [y0, y1], 'k-', linewidth=0.7)
				ax.plot([x0, x1], [y1, y1], 'k-', linewidth=0.7)
				ax.plot([x1, x1], [y1, y0], 'k-', linewidth=0.7)
				ax.text(x_text, y_text, significance[contrast_name][rois], **sign_options)

	# Panel captions
	for ax, caption in zip(axs, ['A', 'B', 'C', 'D']):
		ax.annotate(caption, xy=(-0.18, 1.05),  xycoords="axes fraction", fontsize=10)

	plt.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.88, wspace=0.55)
	fig.set_size_inches(7.5, 2)
	# plt.savefig(f'./figures/figS4-residuals.svg')
	plt.savefig(f'./figures/figS4-residuals.pdf')

	df_d      = pd.DataFrame(d)
	df_p      = pd.DataFrame(p)
	df_p_corr = pd.DataFrame(p_corr)

	# Rename columns to simplified names
	col_rename = {
		'adaptation (std0-std1)': 'adaptation',
		'omissions (om4-std1)': 'omissions',
		'sPE (om4-om6)': 'negative PE'
	}
	df_d = df_d.rename(columns=col_rename)
	df_p = df_p.rename(columns=col_rename)
	df_p_corr = df_p_corr.rename(columns=col_rename)

	# Row label formatter: convert "IC-MGB" to "IC $\neq$ MGB"
	def format_row_label(label):
		return label.replace('-', ' $\\neq$ ')

	# Merge and write single LaTeX file
	merged_latex = _merge_latex_tables_(df_d, df_p, df_p_corr,
									  row_label_formatter=format_row_label)
	with open('./figures/tabS4-residuals.tex', 'w') as f:
		f.write(merged_latex)


# Create figures directory if it doesn't exist
os.makedirs('figures', exist_ok=True)

rois       = ['ICL', 'ICR', 'MGBL', 'MGBR', 'aHGL', 'aHGR']
conditions = ['std0', 'std1', 'std2', 'om4', 'om5', 'om6']


with open(f'./results/zscores_t1w.pickle', 'rb') as f:
	z = pickle.load(f)
	subjects = list(z.keys())

with open(f'./results/funcloc_t1w.pickle', 'rb') as f:
	funcloc = pickle.load(f)

with open(f'./results/noise_t1w.pickle', 'rb') as f:
	noise = pickle.load(f)


# Main figures
fig2_omres(z, subjects, rois, conditions)
fig3_rois(z, subjects, rois)
fig4_emergence(z, funcloc, subjects, rois, conditions)

# Supplementary figures
figS2_allregs(z, subjects, rois, conditions)
figS3_participants(z, subjects, rois, conditions)
figS4_residuals(z, noise, subjects, rois, conditions)
figS5_sweeps(z, funcloc, subjects, rois, conditions)  # Requires sweeps data

