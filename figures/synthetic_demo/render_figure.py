#!/usr/bin/env python3
"""Render the manuscript-layout demonstration from explicitly synthetic tables.

No scientific narrative, methods, or provenance is burned into the artwork;
these belong in FIGURE_LEGEND.md. All plotted observations come from data/.
"""
from pathlib import Path
import base64
import json
import os
import tempfile

os.environ.setdefault('MPLCONFIGDIR', '/tmp/lipohypolep-mpl-v2')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.text import Text
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from fontTools.fontBuilder import FontBuilder
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.pens.cu2quPen import Cu2QuPen
from fontTools.ttLib import TTFont
import fitz
from PIL import Image

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent
DATA = ROOT/'data'
W,H,DPI = 1536,1024,100
NAME = 'figure_01_synthetic'
GROUPS = ('LD SALINE','LD LEPTIN','NON-LD control')
PALETTE = {'LD SALINE':'#C28E85','LD LEPTIN':'#80A596','NON-LD control':'#8B9FB3'}
MARKERS = {'LD SALINE':'o','LD LEPTIN':'s','NON-LD control':'^'}
INK,AXIS,GRID,NS = '#262A2B','#A1A7A7','#C8CDCD','#CFD3D4'
LOWER,HIGHER = '#779AAC','#C28A80'

def prepare_nimbus():
    """Retain Nimbus outlines/metrics while embedding clean TrueType PDFs.

    Matplotlib's CFF OTF embedding produces font-type mismatch warnings. This
    standard cubic-to-quadratic conversion uses the actual URW Nimbus fonts,
    with less than 0.5/1000 em curve error; no font substitution takes place.
    Derivative font files are cached outside the repository.
    """
    source=Path(os.environ.get('NIMBUS_SANS_FONT_DIR','/usr/share/fonts/opentype/urw-base35'))
    cache=Path(tempfile.gettempdir())/'lipohypolep-nimbus-v2'
    cache.mkdir(exist_ok=True)
    # Ensure mathtext and family-based lookups use the same verified TTF files.
    fm.fontManager.ttflist=[entry for entry in fm.fontManager.ttflist if entry.name!='Nimbus Sans']
    result={}
    for style in ('Regular','Bold','Italic'):
        otf=source/f'NimbusSans-{style}.otf'
        if not otf.is_file():
            raise RuntimeError(f'Required Nimbus Sans font missing: {otf}; substitution is forbidden.')
        ttf=cache/f'NimbusSans-{style}.ttf'
        f=TTFont(otf)
        gs=f.getGlyphSet()
        glyphs={}
        for name in f.getGlyphOrder():
            pen=TTGlyphPen(gs)
            gs[name].draw(Cu2QuPen(pen,max_err=0.5,reverse_direction=True))
            glyphs[name]=pen.glyph()
        b=FontBuilder(f['head'].unitsPerEm,isTTF=True)
        b.setupGlyphOrder(f.getGlyphOrder())
        b.setupCharacterMap(f.getBestCmap())
        b.setupGlyf(glyphs)
        b.setupHorizontalMetrics({k:f['hmtx'].metrics[k] for k in f.getGlyphOrder()})
        b.setupHorizontalHeader(ascent=f['hhea'].ascent,descent=f['hhea'].descent,lineGap=f['hhea'].lineGap)
        b.setupNameTable({'familyName':'Nimbus Sans','styleName':style,
                         'uniqueFontIdentifier':f'NimbusSans-{style}-Figure',
                         'fullName':f'Nimbus Sans {style}','psName':f'NimbusSans-{style}'})
        b.setupOS2(sTypoAscender=f['OS/2'].sTypoAscender,sTypoDescender=f['OS/2'].sTypoDescender,
                   usWinAscent=f['OS/2'].usWinAscent,usWinDescent=f['OS/2'].usWinDescent,
                   fsSelection=f['OS/2'].fsSelection,usWeightClass=f['OS/2'].usWeightClass)
        b.setupPost(italicAngle=f['post'].italicAngle)
        b.setupMaxp()
        b.font['head'].macStyle=f['head'].macStyle
        b.save(ttf)
        assert fm.FontProperties(fname=str(ttf)).get_name()=='Nimbus Sans'
        fm.fontManager.addfont(str(ttf))
        result[style]=ttf
    return result

def px(size): return size*72/DPI

def main():
    fonts=prepare_nimbus()
    plt.rcParams.update({
        'font.family':'Nimbus Sans','font.size':px(22),
        'axes.labelsize':px(24),'axes.titlesize':px(24),
        'xtick.labelsize':px(21),'ytick.labelsize':px(21),
        'legend.fontsize':px(20),'axes.labelcolor':INK,
        'text.color':INK,'xtick.color':INK,'ytick.color':INK,
        'axes.edgecolor':AXIS,'axes.linewidth':0.65,
        'axes.spines.top':False,'axes.spines.right':False,
        'axes.titleweight':'normal','xtick.major.width':0.65,'ytick.major.width':0.65,
        'xtick.major.size':3,'ytick.major.size':3,'axes.unicode_minus':True,
        'mathtext.fontset':'custom','mathtext.rm':'Nimbus Sans',
        'mathtext.it':'Nimbus Sans:italic','mathtext.bf':'Nimbus Sans:bold',
        'mathtext.cal':'Nimbus Sans','mathtext.sf':'Nimbus Sans','mathtext.tt':'Nimbus Sans',
        'mathtext.default':'regular','mathtext.fallback':None,
        'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'path',
        'savefig.facecolor':'white','figure.facecolor':'white',
    })
    fig=plt.figure(figsize=(W/DPI,H/DPI),dpi=DPI)
    def axes(left,top,width,height):
        return fig.add_axes([left/W,1-(top+height)/H,width/W,height/H])
    def panel(letter,left,top):
        fig.text(left/W,1-top/H,letter,fontproperties=fm.FontProperties(fname=str(fonts['Bold']),size=px(32)),
                 ha='left',va='top')
    def theme(ax):
        ax.tick_params(direction='out',pad=5)
        ax.grid(False)
        ax.set_axisbelow(True)
    contrasts=[
        ('LD LEPTIN','LD SALINE','contrast_LD_LEPTIN_vs_LD_SALINE.tsv.gz'),
        ('LD SALINE','NON-LD control','contrast_LD_SALINE_vs_NON_LD_control.tsv.gz'),
        ('LD LEPTIN','NON-LD control','contrast_LD_LEPTIN_vs_NON_LD_control.tsv.gz'),
    ]
    tables=[pd.read_csv(DATA/f,sep='\t').set_index('gene_id') for _,_,f in contrasts]
    for table in tables[1:]: assert table.index.equals(tables[0].index)
    xmax=float(np.ceil(max(t.log2FoldChange.abs().max() for t in tables)))
    ymax=float(np.ceil(max((-np.log10(t.padj)).max() for t in tables)/2)*2)
    top_axes=[]
    for i,((first,second,_),table) in enumerate(zip(contrasts,tables)):
        left=110+510*i
        ax=axes(left,88,352,344)
        top_axes.append(ax)
        panel('abc'[i],left-70,30)
        ax.set_title(f'{first} vs {second}',pad=20)
        x=table.log2FoldChange.to_numpy()
        y=-np.log10(table.padj.to_numpy())
        selected=(table.padj.to_numpy()<0.05)&(np.abs(x)>=1)
        ax.scatter(x[~selected],y[~selected],s=8,c=NS,alpha=.55,edgecolors='none',rasterized=True)
        for mask,color in [(selected&(x<0),LOWER),(selected&(x>0),HIGHER)]:
            ax.scatter(x[mask],y[mask],s=10,c=color,alpha=.82,edgecolors='none',rasterized=True)
        for value in [-1,1]: ax.axvline(value,color=GRID,lw=.65,ls=(0,(3,4)),zorder=0)
        ax.axhline(-np.log10(.05),color=GRID,lw=.65,ls=(0,(3,4)),zorder=0)
        ax.set(xlim=(-xmax,xmax),ylim=(-.08,ymax),xticks=[-4,-2,0,2,4])
        ax.set_xlabel(r'log$_2$ fold change',labelpad=8)
        ax.set_ylabel(r'−log$_{10}$(adjusted $P$)',labelpad=8)
        ax.set_yticks(np.arange(0,ymax+.1,2))
        theme(ax)

    # Sample-level PCA, with actual observations and an unframed semantic key.
    ax=axes(110,610,352,302)
    panel('d',40,526)
    scores=pd.read_csv(DATA/'pca_scores.tsv',sep='\t')
    variance=pd.read_csv(DATA/'pca_variance.tsv',sep='\t').set_index('component')
    for group in GROUPS:
        q=scores[scores.group==group]
        ax.scatter(q.PC1,q.PC2,s=62,marker=MARKERS[group],c=PALETTE[group],
                   edgecolors='white',linewidths=.65,zorder=3)
    ax.set_xlabel(f'PC1 ({variance.loc["PC1","explained_variance_ratio"]*100:.1f}%)',labelpad=8)
    ax.set_ylabel(f'PC2 ({variance.loc["PC2","explained_variance_ratio"]*100:.1f}%)',labelpad=8)
    ax.set(xlim=(-55,55),ylim=(-35,46),xticks=[-40,0,40],yticks=[-20,0,20,40])
    handles=[Line2D([],[],marker=MARKERS[g],color='none',markerfacecolor=PALETTE[g],
                     markeredgecolor='white',markersize=7,label=g) for g in GROUPS]
    # Three compact semantic labels, no extra study narration.
    ax.legend(handles=handles,loc='lower left',bbox_to_anchor=(-.03,1.02),
              ncol=1,frameon=False,handlelength=1,handletextpad=.35,labelspacing=.25,
              borderaxespad=0,fontsize=px(19))
    theme(ax)

    # Joint effect display: correlated estimates, deliberately descriptive.
    ax=axes(642,610,352,302)
    panel('e',552,526)
    x=tables[1].log2FoldChange.to_numpy()
    y=tables[2].log2FoldChange.to_numpy()
    ld_selected=(tables[1].padj.to_numpy()<.05)&(np.abs(x)>=1)
    ax.scatter(x[~ld_selected],y[~ld_selected],s=8,c=NS,alpha=.35,edgecolors='none',rasterized=True)
    ax.scatter(x[ld_selected],y[ld_selected],s=10,c=PALETTE['LD LEPTIN'],alpha=.62,edgecolors='none',rasterized=True)
    ax.axhline(0,lw=.65,color=GRID,zorder=0)
    ax.axvline(0,lw=.65,color=GRID,zorder=0)
    ax.plot([-xmax,xmax],[-xmax,xmax],color=AXIS,lw=.8,ls=(0,(4,4)),zorder=0)
    ax.set(xlim=(-xmax,xmax),ylim=(-xmax,xmax),xticks=[-4,-2,0,2,4],yticks=[-4,-2,0,2,4])
    ax.set_xlabel('LD SALINE vs NON-LD control\n'+r'(log$_2$ fold change)',fontsize=px(21),labelpad=10)
    ax.set_ylabel('LD LEPTIN vs NON-LD control\n'+r'(log$_2$ fold change)',fontsize=px(21),labelpad=10)
    theme(ax)

    # Selection is entirely observation-based; truth classes are not used here.
    expression=pd.read_csv(DATA/'normalized_log2_expression.tsv.gz',sep='\t').set_index('gene_id')
    metadata=pd.read_csv(DATA/'sample_metadata.tsv',sep='\t').set_index('sample_id')
    ordered_samples=[s for g in GROUPS for s in metadata.index[metadata.group==g]]
    minima=np.minimum.reduce([t.padj.to_numpy() for t in tables])
    ordered=np.argsort(minima,kind='stable')[:60]
    selected_genes=tables[0].index[ordered]
    matrix=expression.loc[selected_genes,ordered_samples].to_numpy()
    z=(matrix-matrix.mean(axis=1,keepdims=True))/matrix.std(axis=1,keepdims=True,ddof=1)
    rows=leaves_list(linkage(z,method='average',metric='euclidean',optimal_ordering=True))
    z=z[rows]
    selected_genes=selected_genes[rows]
    panel('f',1062,526)
    ax=axes(1150,610,302,302)
    cmap=LinearSegmentedColormap.from_list('soft_expression',[LOWER,'#FAF9F6',HIGHER])
    im=ax.imshow(z,aspect='auto',interpolation='nearest',cmap=cmap,vmin=-2,vmax=2)
    ax.set_xticks([]);ax.set_yticks([])
    ax.set_ylabel('Genes',labelpad=11)
    for spine in ax.spines.values():spine.set_visible(False)
    for value in [5.5,11.5]:ax.axvline(value,color='white',lw=2)
    # Full names span two lines solely for legibility in narrow annotation blocks.
    for i,group in enumerate(GROUPS):
        left=1150+i*302/3
        strip=axes(left,600,302/3-2,6)
        strip.add_patch(plt.Rectangle((0,0),1,1,color=PALETTE[group]))
        strip.set_axis_off()
        label={'LD SALINE':'LD\nSALINE','LD LEPTIN':'LD\nLEPTIN','NON-LD control':'NON-LD\ncontrol'}[group]
        fig.text((left+302/6)/W,1-588/H,label,ha='center',va='bottom',
                 fontsize=px(19),linespacing=1.13)
    cax=axes(1464,668,10,186)
    cb=fig.colorbar(im,cax=cax,ticks=[-2,0,2])
    cb.outline.set_visible(False)
    cb.ax.tick_params(size=2,pad=3,labelsize=px(19))
    cb.set_label('z-score',fontsize=px(21),labelpad=5)
    pd.DataFrame({'gene_id':selected_genes}).to_csv(DATA/'heatmap_genes.tsv',sep='\t',index=False)

    # Lock fonts on all ticks, legend elements, labels, and titles.
    labels=[]
    for t in fig.findobj(Text):
        if not t.get_text():continue
        style='Bold' if (t.get_fontweight() in ['bold',700] or
                         t.get_fontproperties().get_file()==str(fonts['Bold'])) else 'Regular'
        t.set_fontproperties(fm.FontProperties(fname=str(fonts[style]),size=t.get_fontsize()))
        labels.append({'text':t.get_text(),'font':f'Nimbus Sans {style}','size_px':round(t.get_fontsize()*DPI/72,2)})
    fig.canvas.draw()
    for t in fig.findobj(Text):
        if not t.get_text() or not t.get_visible():continue
        b=t.get_window_extent(fig.canvas.get_renderer())
        if b.x0<0 or b.y0<0 or b.x1>W or b.y1>H:
            raise RuntimeError(f'Label outside canvas: {t.get_text()!r}, {b.bounds}')
    common={'Title':'Synthetic hypothalamic RNA-seq figure demonstration',
            'Description':'Entirely synthetic data. See separate figure legend; not experimental results.'}
    fig.savefig(OUT/f'{NAME}.png',dpi=100,metadata=common)
    fig.savefig(OUT/f'{NAME}_300dpi.png',dpi=300,metadata=common)
    fig.savefig(OUT/f'{NAME}.pdf',metadata={'Title':common['Title'],'Subject':common['Description']})
    fig.savefig(OUT/f'{NAME}.svg',metadata=common)
    svg=OUT/f'{NAME}.svg'
    svg.write_text('\n'.join(s.rstrip() for s in svg.read_text().splitlines())+'\n')
    assert '<text' not in svg.read_text()
    with matplotlib.rc_context({'svg.fonttype':'none'}):
        fig.savefig(OUT/f'{NAME}_editable.svg',metadata=common)
    editable=OUT/f'{NAME}_editable.svg'
    code=editable.read_text()
    css=[]
    for style,weight,fontstyle in [('Regular',400,'normal'),('Bold',700,'normal'),('Italic',400,'italic')]:
        encoded=base64.b64encode(fonts[style].read_bytes()).decode()
        css.append(f'@font-face {{ font-family: "Nimbus Sans"; font-weight: {weight}; font-style: {fontstyle}; src: url("data:font/ttf;base64,{encoded}") format("truetype"); }}')
    code=code.replace('<defs>','<defs>\n<style type="text/css"><![CDATA[\n'+'\n'.join(css)+'\n]]></style>',1)
    editable.write_text('\n'.join(s.rstrip() for s in code.splitlines())+'\n')
    plt.close(fig)
    assert Image.open(OUT/f'{NAME}.png').size==(W,H)
    with fitz.open(OUT/f'{NAME}.pdf') as doc:
        pdf_fonts=sorted({f[3] for page in doc for f in page.get_fonts()})
        if not pdf_fonts or not all('NimbusSans' in f for f in pdf_fonts):
            raise RuntimeError(f'Font substitution in PDF: {pdf_fonts}')
    audit={'all_visible_typography':'Nimbus Sans','canvas_pixels':[W,H],
           'pdf_fonts':pdf_fonts,'pdf_text':'embedded, searchable vector text',
           'svg_text':'outlined paths','editable_svg':'text with embedded Nimbus Sans font data',
           'figure_contains_editorial_prose':False,'labels':labels,
           'heatmap_selection':'60 tested genes with smallest minimum adjusted P across three contrasts; average-linkage Euclidean row ordering',
           'volcano_limits':{'absolute_log2fc':xmax,'negative_log10_adjusted_p':ymax}}
    (ROOT/'font_and_plot_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k!='labels'},indent=2))

if __name__=='__main__':main()
