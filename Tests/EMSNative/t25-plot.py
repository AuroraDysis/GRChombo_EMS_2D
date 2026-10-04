#!/usr/bin/env python3
"""Signed numerical expansion envelopes, not horizon-area extrapolations."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
try:
    import scienceplots
    plt.style.use(['science','no-latex'])
except ImportError:
    plt.rcParams.update({'font.family':'serif','font.size':9})
HERE=Path(__file__).resolve().parent
def main():
    fig,axs=plt.subplots(2,2,figsize=(7.2,5.3),layout='constrained')
    cases=[('single','single-sphere','Unboosted single, t=0'),('binary','initial-h0-o+0.0000-c1.00','Binary left, t=0'),
           ('late','late-h0-c1.00','Individual left, t=133'),('late','late-common-c1.00','Midpoint, t=133')]
    for ax,(kind,family,title) in zip(axs.flat,cases):
        rows=[r for r in csv.DictReader((HERE/f't25-map-{kind}-N192.csv').open()) if r['family']==family and r['status']=='RESOLVED']
        r=[float(x['a']) for x in rows];low=[float(x['theta_min']) for x in rows];high=[float(x['theta_max']) for x in rows];mean=[float(x['theta_mean']) for x in rows]
        ax.fill_between(r,low,high,color='#356a9a',alpha=.22,label='angular min–max')
        ax.plot(r,mean,color='#223449',lw=1.2,label='area-weighted mean')
        ax.axhline(0,color='#9b3b2b',lw=.8);ax.set_xscale('log');ax.set_xlim(min(r),max(r))
        ax.set_yscale('symlog',linthresh=.01);ax.set_ylim(min(low)*1.15,max(high)*1.15)
        ax.set_title(title);ax.set_xlabel(r'$r/M_i$');ax.set_ylabel(r'$M_i\theta_+$')
        if kind!='late':ax.axvspan(.00604 if kind=='single' else .00608,.00606 if kind=='single' else .0062,color='#ce963a',alpha=.35)
    axs[0,0].legend(frameon=False,fontsize=7)
    fig.savefig(HERE/'t25-signed-expansion.png',dpi=300);plt.close(fig)
if __name__=='__main__':main()
