#!/usr/bin/env python3
"""Reproduce palette and full-figure color-vision simulations (not certification)."""
from pathlib import Path
import itertools
import json
import os

import numpy as np
from PIL import Image
from daltonlens import simulate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

ROOT=Path(__file__).resolve().parent
PALETTES={
    'Original groups':{'LD SALINE':'#C28E85','LD LEPTIN':'#80A596','NON-LD control':'#8B9FB3'},
    'Revised groups':{'LD SALINE':'#C28E85','LD LEPTIN':'#527D6B','NON-LD control':'#A8BED2'},
    'Effect direction':{'Lower':'#779AAC','Higher':'#C28A80','Not significant':'#CFD3D4'},
}
CONDITIONS={'Normal':None,'Protanopia':simulate.Deficiency.PROTAN,
            'Deuteranopia':simulate.Deficiency.DEUTAN,'Tritanopia':simulate.Deficiency.TRITAN,
            'Grayscale':'gray'}

def rgb(h):return np.array([int(h[i:i+2],16) for i in (1,3,5)],np.uint8)

def linear(a):
    a=np.asarray(a)/255
    return np.where(a<=.04045,a/12.92,((a+.055)/1.055)**2.4)

def lab(a):
    xyz=linear(a)@np.array([[.4124564,.3575761,.1804375],
                           [.2126729,.7151522,.072175],
                           [.0193339,.119192,.9503041]]).T
    xyz=xyz/np.array([.95047,1,1.08883])
    f=np.where(xyz>(6/29)**3,np.cbrt(xyz),xyz/(3*(6/29)**2)+4/29)
    return np.stack([116*f[...,1]-16,500*(f[...,0]-f[...,1]),200*(f[...,1]-f[...,2])],axis=-1)

def transform(a,condition,simulator):
    if condition is None:return a.copy()
    if condition=='gray':
        y=linear(a)@np.array([.2126729,.7151522,.072175])
        encoded=np.where(y<=.0031308,y*12.92,1.055*y**(1/2.4)-.055)
        return np.repeat(np.round(encoded[...,None]*255).astype(np.uint8),3,axis=-1)
    return simulator.simulate_cvd(a,condition,1.0)

def main():
    simulator=simulate.Simulator_Machado2009()
    regular=Path(os.environ.get('NIMBUS_SANS_FONT_DIR','/usr/share/fonts/opentype/urw-base35'))/'NimbusSans-Regular.otf'
    if not regular.is_file():raise RuntimeError('Nimbus Sans is required; no substitution.')
    fm.fontManager.addfont(str(regular))
    fp=fm.FontProperties(fname=str(regular),size=13)
    fig,axes=plt.subplots(3,5,figsize=(17.6,8.8),dpi=120)
    results={}
    for row,(name,palette) in enumerate(PALETTES.items()):
        values=np.array([rgb(h) for h in palette.values()])[None,:,:]
        results[name]={}
        for col,(label,condition) in enumerate(CONDITIONS.items()):
            values_seen=transform(values,condition,simulator)[0]
            labs=lab(values_seen)
            result={'hex':{key:'#'+''.join(f'{int(v):02X}' for v in color)
                           for key,color in zip(palette,values_seen)},
                    'Lstar':{key:round(float(v[0]),2) for key,v in zip(palette,labs)},
                    'deltaE76':{f'{list(palette)[i]} / {list(palette)[j]}':round(float(np.linalg.norm(labs[i]-labs[j])),2)
                                 for i,j in itertools.combinations(range(len(palette)),2)}}
            results[name][label]=result
            ax=axes[row,col]
            for j,(key,color) in enumerate(zip(palette,values_seen)):
                ax.add_patch(plt.Rectangle((0,2-j),1,.54,color=color/255))
                ax.text(.5,2.73-j,key,fontproperties=fp,ha='center',va='center')
            ax.set(xlim=(0,1),ylim=(-.1,3.05));ax.axis('off')
            if row==0:ax.set_title(label,fontproperties=fp,pad=18)
            if col==0:ax.text(-.12,.5,name,transform=ax.transAxes,rotation=90,ha='right',va='center',fontproperties=fp)
    fig.tight_layout(pad=2.5,w_pad=2,h_pad=2)
    fig.savefig(ROOT/'palette_color_vision_check.png')
    plt.close(fig)
    original=np.asarray(Image.open(ROOT.parent/'figure_01_small_effects_synthetic.png').convert('RGB'))
    for label,condition in CONDITIONS.items():
        if condition is None:continue
        seen=transform(original,condition,simulator)
        Image.fromarray(seen).save(ROOT/f'figure_01_{label.lower()}_check.png')
    report={'method':'DaltonLens 0.1.5, Machado et al. 2009, severity 1.0; sRGB decoded and encoded by library',
            'grayscale':'linear-light sRGB/Rec.709 relative luminance, re-encoded to sRGB',
            'metrics':'CIE1976 Lab distances (D65); descriptive only, no universal pass/fail threshold',
            'scope':'Original and revised categorical colors; retained effect hues; full rendered Figure 1 including heatmap and transparency',
            'limitation':'Simulation is an accessibility check, not a guarantee for every observer or grayscale print.',
            'palettes':results}
    (ROOT/'palette_simulation.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Created actual palette and full-Figure-1 protanopia/deuteranopia/tritanopia/grayscale checks.')

if __name__=='__main__':main()
