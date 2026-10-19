#!/usr/bin/env python3
"""Reproduce palette and full-figure color-vision simulations (not certification)."""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import io
import itertools
import json
import os
import tempfile

import numpy as np
from PIL import Image
from daltonlens import simulate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

ROOT=Path(__file__).resolve().parent
REPORT=ROOT/'palette_simulation.json'
FONT_DIR=Path(os.environ.get('NIMBUS_SANS_FONT_DIR','/usr/share/fonts/opentype/urw-base35'))
PALETTES={
    'Figure 1: yellow / blue':{'Higher':'#FDB338','Lower':'#025196'},
    'Figure 2: tan / turquoise':{'Higher':'#FFBE6A','Lower':'#40B0A6'},
    'Figure 3: orange / purple':{'Higher':'#EB6123','Lower':'#512888'},
}
CONDITIONS={'Normal':None,'Protanopia':simulate.Deficiency.PROTAN,
            'Deuteranopia':simulate.Deficiency.DEUTAN,'Tritanopia':simulate.Deficiency.TRITAN,
            'Grayscale':'gray'}

def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def atomic_bytes(path, data):
    with tempfile.NamedTemporaryFile('wb',dir=path.parent,delete=False) as handle:
        temporary=Path(handle.name)
        try:
            handle.write(data)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    temporary.replace(path)

def write_report(report):
    atomic_bytes(REPORT,(json.dumps(report,indent=2,allow_nan=False)+'\n').encode())

def source_snapshot():
    paths=[Path(__file__),FONT_DIR/'NimbusSans-Regular.otf']
    paths.extend(ROOT.parent/f'figure_{number:02d}_gsea_synthetic.png' for number in (1,2,3))
    for path in paths:
        if not path.is_file():raise RuntimeError(f'Missing color-audit input: {path}')
    return {str(path):sha256(path) for path in paths}

def verify_report_current(path=REPORT):
    report=json.loads(path.read_text())
    if report.get('status')!='PASS':raise RuntimeError('Color report is not a completed PASS')
    if report.get('auditor_sha256')!=sha256(Path(__file__)):
        raise RuntimeError('Color-audit script changed; rerun the audit')
    regular=FONT_DIR/'NimbusSans-Regular.otf'
    if not regular.is_file() or report.get('font_sha256')!=sha256(regular):
        raise RuntimeError('Color-audit source font changed or missing')
    inputs=report.get('input_export_sha256',{})
    expected={f'figure_{number:02d}_gsea_synthetic.png' for number in (1,2,3)}
    if set(inputs)!=expected:raise RuntimeError('Incomplete color-audit input evidence')
    for filename,digest in inputs.items():
        source=ROOT.parent/filename
        if not source.is_file() or sha256(source)!=digest:
            raise RuntimeError(f'Color-audit input changed or missing: {filename}')
    outputs=report.get('output_sha256',{})
    expected_outputs={'palette_color_vision_check.png'} | {
        f'figure_{number:02d}_{label.lower()}_check.png' for number in (1,2,3)
        for label,condition in CONDITIONS.items() if condition is not None}
    if set(outputs)!=expected_outputs:raise RuntimeError('Incomplete color-audit output evidence')
    for filename,digest in outputs.items():
        output=ROOT/filename
        if not output.is_file() or sha256(output)!=digest:
            raise RuntimeError(f'Color-audit output changed or missing: {filename}')
    return report

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

def run_audit():
    simulator=simulate.Simulator_Machado2009()
    regular=FONT_DIR/'NimbusSans-Regular.otf'
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
    output=io.BytesIO()
    fig.savefig(output,format='png')
    atomic_bytes(ROOT/'palette_color_vision_check.png',output.getvalue())
    plt.close(fig)
    input_hashes={}
    output_files=[ROOT/'palette_color_vision_check.png']
    for number in (1,2,3):
        path=ROOT.parent/f'figure_{number:02d}_gsea_synthetic.png'
        data=path.read_bytes()
        input_hashes[path.name]=hashlib.sha256(data).hexdigest()
        with Image.open(io.BytesIO(data)) as image:
            if image.size!=(1536,1024):raise RuntimeError(f'Unexpected source figure size: {path}: {image.size}')
            original=np.asarray(image.convert('RGB'))
        for label,condition in CONDITIONS.items():
            if condition is None:continue
            seen=transform(original,condition,simulator)
            path=ROOT/f'figure_{number:02d}_{label.lower()}_check.png'
            output=io.BytesIO()
            Image.fromarray(seen).save(output,format='PNG')
            atomic_bytes(path,output.getvalue())
            output_files.append(path)
    report={'status':'PASS','verified_at_utc':datetime.now(timezone.utc).isoformat(),
            'method':f'DaltonLens {version("daltonlens")}, Machado et al. 2009, severity 1.0; sRGB decoded and encoded by library',
            'grayscale':'linear-light sRGB/Rec.709 relative luminance, re-encoded to sRGB',
            'metrics':'CIE1976 Lab distances (D65); descriptive only, no universal pass/fail threshold',
            'scope':'Three comparison-specific pairs from the RGB codes printed in the supplied PDF; full rendered Figures 1–3 including heatmaps, GSEA and transparency',
            'limitation':'Simulation is an accessibility check, not a guarantee for every observer or grayscale print.',
            'palettes':results,
            'auditor_sha256':sha256(Path(__file__)),
            'font_sha256':sha256(regular),
            'software':{package:version(package) for package in ('daltonlens','matplotlib','numpy','Pillow')},
            'input_export_sha256':input_hashes,
            'output_sha256':{path.name:sha256(path) for path in output_files}}
    return report

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-report',action='store_true',help='Verify evidence matches current sources and simulations.')
    args=parser.parse_args()
    if args.check_report:
        verify_report_current()
        print('PASS: saved color-audit evidence matches current inputs and simulations.')
        return
    write_report({'status':'RUNNING','started_at_utc':datetime.now(timezone.utc).isoformat()})
    try:
        before=source_snapshot()
        report=run_audit()
        if before!=source_snapshot():raise RuntimeError('An input changed during the color audit; rerun')
        write_report(report)
    except BaseException as error:
        write_report({'status':'FAILED','verified_at_utc':datetime.now(timezone.utc).isoformat(),
                      'error':f'{type(error).__name__}: {error}'})
        raise
    print('Created requested-palette and full-Figures-1–3 protanopia/deuteranopia/tritanopia/grayscale checks.')

if __name__=='__main__':main()
