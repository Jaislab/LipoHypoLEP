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
import re
import textwrap
import hashlib
import warnings

os.environ.setdefault('MPLCONFIGDIR', '/tmp/lipohypolep-mpl-v2')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.text import Text
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy import stats
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
# Scale every text role together; keep Nimbus Sans and the artwork dimensions.
FONT_SCALE = 1.20
GROUPS = ('LD SALINE','LD LEPTIN','NON-LD control')
INK,AXIS,GRID,NS = '#262A2B','#A1A7A7','#C8CDCD','#CFD3D4'
LOG2FC_CUTOFF = 0.2
ADJUSTED_P_CUTOFF = 0.05
# RGB values printed beside the requested pairs in the supplied NCEAS guide.
COMPARISON_COLORS = {
    1: {'name':'yellow/blue','higher':'#FDB338','lower':'#025196'},
    2: {'name':'tan/turquoise','higher':'#FFBE6A','lower':'#40B0A6'},
    3: {'name':'orange/purple','higher':'#EB6123','lower':'#512888'},
}

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
        name_table={'familyName':'Nimbus Sans','styleName':style,
                    'uniqueFontIdentifier':f'NimbusSans-{style}-Figure',
                    'fullName':f'Nimbus Sans {style}','psName':f'NimbusSans-{style}'}
        # Keep source attribution and licensing in the converted/embedded fonts.
        for key,identifier in [('copyright',0),('trademark',7),('manufacturer',8),
                               ('designer',9),('description',10),('vendorURL',11),
                               ('designerURL',12),('licenseDescription',13),('licenseInfoURL',14)]:
            value=f['name'].getDebugName(identifier)
            if value:name_table[key]=value
        b.setupNameTable(name_table)
        b.setupOS2(sTypoAscender=f['OS/2'].sTypoAscender,sTypoDescender=f['OS/2'].sTypoDescender,
                   usWinAscent=f['OS/2'].usWinAscent,usWinDescent=f['OS/2'].usWinDescent,
                   fsSelection=f['OS/2'].fsSelection,usWeightClass=f['OS/2'].usWeightClass)
        b.setupPost(italicAngle=f['post'].italicAngle)
        b.setupMaxp()
        b.font['head'].macStyle=f['head'].macStyle
        b.save(ttf)
        if fm.FontProperties(fname=str(ttf)).get_name()!='Nimbus Sans':
            raise RuntimeError(f'Unexpected font family in {ttf}; substitution is forbidden.')
        fm.fontManager.addfont(str(ttf))
        result[style]=ttf
    return result

def px(size): return size*FONT_SCALE*72/DPI

SPECS = (
    (1, 'How does untreated LD differ from non-LD controls?',
     'LD SALINE', 'NON-LD control', 'contrast_LD_SALINE_vs_NON_LD_control.tsv.gz'),
    (2, 'What does leptin treatment change in LD mice?',
     'LD LEPTIN', 'LD SALINE', 'contrast_LD_LEPTIN_vs_LD_SALINE.tsv.gz'),
    (3, 'What differences remain after leptin treatment?',
     'LD LEPTIN', 'NON-LD control', 'contrast_LD_LEPTIN_vs_NON_LD_control.tsv.gz'),
)


def configure_theme():
    # A missing character must not silently disappear or become an empty box.
    warnings.filterwarnings('error',message=r'Glyph .* missing from font',category=UserWarning)
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


def selected_rows(table):
    # Gene-level highlights use the user-specified joint criteria; GSEA is unfiltered.
    selected=table[(table.padj<ADJUSTED_P_CUTOFF)&(table.log2FoldChange.abs()>=LOG2FC_CUTOFF)]
    heatmap=selected.sort_values('padj',kind='stable').head(40)
    positive=selected[selected.log2FoldChange>0].sort_values('padj',kind='stable').head(6)
    negative=selected[selected.log2FoldChange<0].sort_values('padj',kind='stable').head(6)
    forest=pd.concat([positive,negative]).sort_values('log2FoldChange',ascending=False,kind='stable')
    return heatmap,forest


def read_table(path):
    """Preserve the binary64 values written with the tables' 17-digit format."""
    return pd.read_csv(path,sep='\t',float_precision='round_trip')


def validate_inputs(tables, expression, metadata):
    """Fail before drawing if data cannot be aligned without guessing or dropping rows."""
    if len(tables)!=len(SPECS):
        raise ValueError('Exactly three contrast tables are required.')
    def check_index(index, label):
        if index.empty or index.has_duplicates or index.isna().any():
            raise ValueError(f'{label} identifiers must be nonempty and unique.')
        if any(not isinstance(value,str) or not value.strip() for value in index):
            raise ValueError(f'{label} identifiers must be nonempty strings.')
    for number,table in enumerate(tables,1):
        check_index(table.index,f'Figure {number} gene')
        check_index(table.columns,f'Figure {number} column')
        required={'log2FoldChange','lfcSE','padj'}
        if not required.issubset(table.columns):
            raise ValueError(f'Figure {number} is missing {sorted(required-set(table.columns))}.')
        if any(not pd.api.types.is_numeric_dtype(table[column]) or pd.api.types.is_bool_dtype(table[column])
               for column in required):
            raise ValueError(f'Figure {number} estimate and probability columns must be numeric.')
        values=table[list(required)].to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f'Figure {number} contains nonfinite estimates, standard errors, or adjusted P values.')
        if (table.lfcSE<0).any() or ((table.padj<0)|(table.padj>1)).any():
            raise ValueError(f'Figure {number} has invalid standard errors or adjusted P values.')
        if not table.index.equals(tables[0].index):
            raise ValueError('Contrast gene identifiers and order must agree exactly.')
    check_index(expression.index,'Expression gene')
    check_index(metadata.index,'Sample')
    check_index(expression.columns,'Expression sample')
    if not expression.index.equals(tables[0].index):
        raise ValueError('Expression and contrast gene identifiers and order must agree exactly.')
    if set(expression.columns)!=set(metadata.index):
        raise ValueError('Expression samples and metadata sample identifiers must agree exactly.')
    if not np.isfinite(expression.to_numpy(dtype=float)).all():
        raise ValueError('Expression contains nonfinite values.')
    if 'group' not in metadata or set(metadata.group)!=set(GROUPS):
        raise ValueError(f'Metadata must contain all three exact group names: {GROUPS}.')


def volcano_ordinates(padj):
    """Keep zero probabilities drawable without changing source values or selection."""
    values=np.asarray(padj,dtype=float)
    if not np.isfinite(values).all() or ((values<0)|(values>1)).any():
        raise ValueError('Adjusted P values must be finite probabilities.')
    return -np.log10(np.maximum(values,np.nextafter(0.0,1.0)))


def volcano_limits(tables):
    # Minimum ranges retain visible cutoff guides even when every effect is zero.
    xmax=max(1.0,float(np.ceil(max(t.log2FoldChange.abs().max() for t in tables))))
    ymax=max(2.0,float(np.ceil(max(volcano_ordinates(t.padj).max() for t in tables)/2)*2))
    return xmax,ymax


def heatmap_values(expression, genes, samples):
    """Return row z scores; constant rows receive neutral zero rather than NaN."""
    matrix=expression.loc[genes,samples].to_numpy(dtype=float)
    if matrix.shape[0]==0:
        return matrix,[],[]
    centered=matrix-matrix.mean(axis=1,keepdims=True)
    scale=matrix.std(axis=1,keepdims=True,ddof=1)
    z=np.divide(centered,scale,out=np.zeros_like(centered),where=scale>0)
    rows=leaves_list(linkage(z,method='average',metric='euclidean',optimal_ordering=True)) if len(z)>1 else np.array([0])
    constant=list(np.asarray(genes)[scale[:,0]==0])
    return z[rows],list(np.asarray(genes)[rows]),constant


def displayed_gsea(results):
    required={'term','NES','gsea_fdr_q'}
    if not required.issubset(results.columns):
        raise ValueError(f'GSEA display table is missing {sorted(required-set(results.columns))}.')
    if results.term.isna().any() or results.term.duplicated().any() or any(not str(term).strip() for term in results.term):
        raise ValueError('GSEA term labels must be nonempty and unique.')
    values=results[['NES','gsea_fdr_q']].to_numpy(dtype=float)
    if not np.isfinite(values).all() or ((results.gsea_fdr_q<0)|(results.gsea_fdr_q>1)).any():
        raise ValueError('GSEA display scores and q values must be finite and q values must lie in [0,1].')
    positive=results[results.NES>0].sort_values('NES',ascending=False,kind='stable').head(2)
    negative=results[results.NES<0].sort_values('NES',ascending=True,kind='stable').head(2)
    return pd.concat([positive,negative]).sort_values('NES',ascending=False,kind='stable')


def q_label(value):
    """Avoid rounding a significant q to 0.05, or a positive estimate to zero."""
    if 0<value<.05 and float(f'{value:.2f}')>=.05:
        return '<0.05'
    if 0<value<.005:
        return f'{value:.1g}'
    return '0' if value==0 else f'{value:.2f}'


def validate_gsea_layout(fig):
    """Reject crowded terms at the requested font size before saving any artwork."""
    renderer=fig.canvas.get_renderer()
    for region in getattr(fig,'_gsea_layout_regions',[]):
        ax=region['axes']
        labels=[label for label in ax.get_yticklabels() if label.get_visible() and label.get_text()]
        boxes=[label.get_window_extent(renderer) for label in labels]
        problem=None
        for label,box in zip(labels,boxes):
            if box.x0<region['left'] or box.x1>ax.bbox.x0-4 or box.y0<ax.bbox.y0 or box.y1>ax.bbox.y1:
                problem=f'term label exceeds its reserved panel area: {label.get_text()!r}'
                break
        if problem is None:
            for i,box in enumerate(boxes):
                if any(box.padded(2).overlaps(other) for other in boxes[i+1:]):
                    problem=f'term labels overlap around {labels[i].get_text()!r}'
                    break
        if problem is None:
            for text in region['q_labels']:
                box=text.get_window_extent(renderer)
                if box.x0<ax.bbox.x1+4 or box.x1>region['right']:
                    problem=f'q label exceeds its reserved column: {text.get_text()!r}'
                    break
        if problem:
            raise RuntimeError(f"{region['collection']} GSEA layout: {problem}. "
                               'Review the term labels and allocate more panel space or fewer displayed terms; '
                               'font size and full term text have been preserved.')


def export_figure(fig,number,question,fonts,details):
    name=f'figure_{number:02d}_gsea_synthetic'
    labels=[]
    for t in fig.findobj(Text):
        if not t.get_text():continue
        style='Bold' if (t.get_fontweight() in ['bold',700] or
                         t.get_fontproperties().get_file()==str(fonts['Bold'])) else 'Regular'
        t.set_fontproperties(fm.FontProperties(fname=str(fonts[style]),size=t.get_fontsize()))
        labels.append({'text':t.get_text(),'font':f'Nimbus Sans {style}',
                       'size_px':round(t.get_fontsize()*DPI/72,2)})
    fig.canvas.draw()
    validate_gsea_layout(fig)
    for t in fig.findobj(Text):
        if not t.get_text() or not t.get_visible():continue
        b=t.get_window_extent(fig.canvas.get_renderer())
        if b.x0<0 or b.y0<0 or b.x1>W or b.y1>H:
            raise RuntimeError(f'Label outside canvas: {t.get_text()!r}, {b.bounds}')
    metadata={'Title':f'FIGURE {number} | {question} | Synthetic demonstration',
              'Description':'Entirely synthetic data. See separate figure legend; not experimental results.'}
    fig.savefig(OUT/f'{name}.png',dpi=100,metadata=metadata)
    fig.savefig(OUT/f'{name}_300dpi.png',dpi=300,metadata=metadata)
    fig.savefig(OUT/f'{name}.pdf',metadata={'Title':metadata['Title'],'Subject':metadata['Description']})
    fig.savefig(OUT/f'{name}.svg',metadata=metadata)
    svg=OUT/f'{name}.svg'
    svg.write_text('\n'.join(s.rstrip() for s in svg.read_text().splitlines())+'\n')
    if '<text' in svg.read_text():
        raise RuntimeError('Outlined SVG unexpectedly contains editable text.')
    with matplotlib.rc_context({'svg.fonttype':'none'}):
        fig.savefig(OUT/f'{name}_editable.svg',metadata=metadata)
    editable=OUT/f'{name}_editable.svg'
    code=editable.read_text()
    css=[]
    for style,weight,fontstyle in [('Regular',400,'normal'),('Bold',700,'normal'),('Italic',400,'italic')]:
        encoded=base64.b64encode(fonts[style].read_bytes()).decode()
        css.append(f'@font-face {{ font-family: "Nimbus Sans"; font-weight: {weight}; font-style: {fontstyle}; src: url("data:font/ttf;base64,{encoded}") format("truetype"); }}')
    code=code.replace('<defs>','<defs>\n<style type="text/css"><![CDATA[\n'+'\n'.join(css)+'\n]]></style>',1)
    editable.write_text('\n'.join(s.rstrip() for s in code.splitlines())+'\n')
    plt.close(fig)
    with Image.open(OUT/f'{name}.png') as image:
        if image.size!=(W,H):
            raise RuntimeError(f'Unexpected raster dimensions: {image.size}.')
    with fitz.open(OUT/f'{name}.pdf') as doc:
        pdf_fonts=sorted({f[3] for page in doc for f in page.get_fonts()})
        if not pdf_fonts or not all('NimbusSans' in f for f in pdf_fonts):
            raise RuntimeError(f'Font substitution in PDF: {pdf_fonts}')
    return {'figure':number,'file_stem':name,'question':question,
            'all_visible_typography':'Nimbus Sans','canvas_pixels':[W,H],
            'typography_scale':FONT_SCALE,
            'pdf_fonts':pdf_fonts,'pdf_text':'embedded, searchable vector text',
            'svg_text':'outlined paths','editable_svg':'text with embedded Nimbus Sans font data',
            'figure_contains_editorial_prose':False,'pca_panel':False,'labels':labels,**details}


def render_question(spec,table,expression,metadata,fonts,xmax,ymax,forest_limit,df,nes_limit=2.6):
    number,question,first,second,filename=spec
    heatmap,forest=selected_rows(table)
    comparison=f'{first} vs {second}'
    colors=COMPARISON_COLORS[number]
    lower,higher=colors['lower'],colors['higher']
    fig=plt.figure(figsize=(W/DPI,H/DPI),dpi=DPI)
    def axes(left,top,width,height):
        return fig.add_axes([left/W,1-(top+height)/H,width/W,height/H])
    def panel(letter,left,top):
        fig.text(left/W,1-top/H,letter,
                 fontproperties=fm.FontProperties(fname=str(fonts['Bold']),size=px(32)),
                 ha='left',va='top')
    def theme(ax):
        ax.tick_params(direction='out',pad=5)
        ax.grid(False)
        ax.set_axisbelow(True)

    # a: the complete single-contrast result, with common scales across figures.
    ax=axes(110,88,560,310)
    panel('a',40,27)
    ax.set_title(comparison,pad=20)
    x=table.log2FoldChange.to_numpy()
    y=volcano_ordinates(table.padj.to_numpy())
    selected=(table.padj.to_numpy()<ADJUSTED_P_CUTOFF)&(np.abs(x)>=LOG2FC_CUTOFF)
    ax.scatter(x[~selected],y[~selected],s=8,c=NS,alpha=.55,edgecolors='none',rasterized=True)
    for mask,color in [(selected&(x<0),lower),(selected&(x>0),higher)]:
        ax.scatter(x[mask],y[mask],s=11,c=color,alpha=.82,edgecolors='none',rasterized=True)
    ax.axvline(0,color=GRID,lw=.65,ls=(0,(3,4)),zorder=0)
    for cutoff in (-LOG2FC_CUTOFF,LOG2FC_CUTOFF):
        ax.axvline(cutoff,color=AXIS,lw=.8,ls=(0,(4,4)),zorder=0)
    ax.axhline(-np.log10(.05),color=GRID,lw=.65,ls=(0,(3,4)),zorder=0)
    ax.set(xlim=(-xmax,xmax),ylim=(-.08,ymax))
    yticks=MaxNLocator(nbins=5,integer=True,steps=[1,2,2.5,5,10]).tick_values(0,ymax)
    ax.set_yticks(yticks[(yticks>=0)&(yticks<=ymax)])
    ticks=MaxNLocator(nbins=6,steps=[1,2,2.5,5,10],symmetric=True).tick_values(-xmax,xmax)
    ax.set_xticks(ticks[(ticks>=-xmax)&(ticks<=xmax)])
    ax.set_xlabel(r'log$_2$ fold change',labelpad=8)
    ax.set_ylabel(r'−log$_{10}$(adjusted $P$)',labelpad=8)
    theme(ax)

    # b: only the two groups that answer this figure's question.
    samples=[s for group in (first,second) for s in metadata.index[metadata.group==group]]
    z,genes,constant_genes=heatmap_values(expression,heatmap.index,samples)
    panel('b',788,27)
    ax=axes(928,132,450,266)
    cmap=LinearSegmentedColormap.from_list(f'expression_{number}',[lower,'#FAF9F6',higher])
    im=ax.imshow(z,aspect='auto',interpolation='nearest',cmap=cmap,vmin=-2,vmax=2) if len(z) else None
    if not len(z):
        ax.text(.5,.5,'No genes meet cutoffs',transform=ax.transAxes,ha='center',va='center',fontsize=px(21))
    ax.set_xticks([]);ax.set_yticks([])
    ax.set_ylabel('Genes',labelpad=12)
    for spine in ax.spines.values():spine.set_visible(False)
    first_count=int((metadata.group==first).sum())
    if len(z):ax.axvline(first_count-.5,color='white',lw=2)
    offset=0
    for i,group in enumerate((first,second)):
        width=450*int((metadata.group==group).sum())/len(samples)
        left=928+offset
        strip=axes(left,120,max(width-2,1),6)
        strip.add_patch(plt.Rectangle((0,0),1,1,color=higher if i==0 else lower))
        strip.set_axis_off()
        fig.text((left+width/2)/W,1-108/H,group,ha='center',va='bottom',fontsize=px(24))
        offset+=width
    if im is not None:
        cax=axes(1402,185,13,160)
        cb=fig.colorbar(im,cax=cax,ticks=[-2,0,2])
        cb.outline.set_visible(False)
        cb.ax.tick_params(size=2,pad=4,labelsize=px(21))
        cb.set_label('z-score',fontsize=px(24),labelpad=7)
    pd.DataFrame({'gene_id':genes}).to_csv(DATA/f'figure_{number:02d}_heatmap_genes.tsv',sep='\t',index=False)

    # Retain the gene effect estimates in tables; the canvas now prioritizes GSEA.
    crit=stats.t.ppf(.975,df)
    errors=crit*forest.lfcSE.to_numpy()
    fx=forest.log2FoldChange.to_numpy()
    forest=forest.copy()
    forest['ci_lower']=fx-errors
    forest['ci_upper']=fx+errors
    forest['comparison']=comparison
    forest.reset_index().to_csv(DATA/f'figure_{number:02d}_forest_genes.tsv',sep='\t',index=False,float_format='%.17g')

    # c,d: one preranked calculation per collection and biological comparison.
    gsea_details={}
    fig._gsea_layout_regions=[]
    # Reserve space for larger q values, including '<0.05' and scientific notation.
    for collection,letter,left,title in [('kegg','c',518,'KEGG'),('go_bp','d',1300,'GO BP')]:
        path=ROOT/'gsea_data'/f'figure_{number:02d}_{collection}_display.tsv'
        if not path.is_file():raise RuntimeError(f'GSEA results missing: run run_gsea.py first ({path})')
        results=read_table(path)
        shown=displayed_gsea(results)
        panel(letter,40 if collection=='kegg' else 788,494)
        fig.text((390 if collection=='kegg' else 1160)/W,1-524/H,title,
                 ha='center',va='baseline',fontsize=px(25))
        ax=axes(left,568,150,340)
        labels=[]
        for term in shown.term:
            label=re.sub(r'^(?:GO:\d+|mmu\d+)\s+','',str(term)).strip()
            labels.append(textwrap.fill(label,width=42,break_long_words=False,break_on_hyphens=False))
        ax.axvline(0,color=GRID,lw=.75,ls=(0,(3,4)),zorder=0)
        q_labels=[]
        for row,(_,result) in enumerate(shown.iterrows()):
            color=higher if result.NES>0 else lower
            is_significant=result.gsea_fdr_q<.05
            ax.scatter(result.NES,row,s=90,facecolors=color if is_significant else 'white',
                       edgecolors=color,linewidths=1.4,zorder=3)
            q_labels.append(ax.text(1.055,row,q_label(result.gsea_fdr_q),transform=ax.get_yaxis_transform(),
                    ha='left',va='center',fontsize=px(20)))
        ax.set_yticks(np.arange(len(shown)),labels=labels,fontsize=px(21))
        if shown.empty:
            ax.text(.5,.5,'No eligible sets',transform=ax.transAxes,ha='center',va='center',fontsize=px(21))
        limit=max(nes_limit,float(shown.NES.abs().max())*1.1) if len(shown) else nes_limit
        ax.set(xlim=(-limit,limit),ylim=(max(len(shown)-.5,.5),-.5))
        if limit==2.6:ax.set_xticks([-2,0,2])
        else:
            ticks=MaxNLocator(nbins=3,steps=[1,2,2.5,5,10],symmetric=True).tick_values(-limit,limit)
            ax.set_xticks(ticks[(ticks>=-limit)&(ticks<=limit)])
        ax.set_xlabel('NES',labelpad=9)
        ax.spines['left'].set_visible(False)
        theme(ax)
        ax.tick_params(axis='y',length=0,pad=12)
        ax.text(1.15,1.04,'q',transform=ax.transAxes,ha='center',va='bottom',fontsize=px(21))
        fig._gsea_layout_regions.append({'axes':ax,'collection':title,'q_labels':q_labels,
                                         'left':40 if collection=='kegg' else 788,
                                         'right':768 if collection=='kegg' else W})
        gsea_details[collection]={'display_table':str(path.relative_to(OUT)),
                                  'terms':shown.term.tolist(),'NES':shown.NES.tolist(),
                                  'gsea_fdr_q':shown.gsea_fdr_q.tolist(),
                                  'selection':'up to two largest positive and two most negative NES, without significance requirement',
                                  'absolute_nes_limit':limit,
                                  'ranking':'all tested genes, signed t statistic; no DEG cutoff',
                                  'synthetic_random_mouse_assignment':True}
    return export_figure(fig,number,question,fonts,{
        'comparison':comparison,'contrast_table':filename,
        'panels':['volcano','two-group heatmap','KEGG preranked GSEA','GO BP preranked GSEA'],
        'volcano_limits':{'absolute_log2fc':xmax,'negative_log10_adjusted_p':ymax},
        'zero_adjusted_p_plotting_floor':float(np.nextafter(0.0,1.0)),
        'volcano_genes':len(table),'heatmap_samples':samples,'heatmap_genes':genes,
        'heatmap_constant_genes_at_zero_zscore':constant_genes,
        'heatmap_group_order':[first,second],
        'heatmap_selection':'up to 40 genes with smallest adjusted P among genes passing adjusted P < 0.05 and absolute log2FC >= 0.2 in this contrast',
        'forest_genes':forest.index.tolist(),'forest_absolute_limit':forest_limit,
        'forest_confidence_level':.95,'forest_residual_degrees_of_freedom':df,
        'forest_intervals':'marginal, unadjusted; genes selected from the same dataset',
        'forest_selection':'up to six positive and six negative estimates with smallest adjusted P among genes passing the joint adjusted-P and log2FC criteria',
        'absolute_log2fc_cutoff':LOG2FC_CUTOFF,'adjusted_p_cutoff':ADJUSTED_P_CUTOFF,
        'synthetic_effect_scale':'original synthetic generator restored; fitted effects and intervals are neither rescaled nor clipped',
        'comparison_colors':colors,'gsea':gsea_details,
    })


def main():
    audit_path=ROOT/'font_and_plot_audit.json'
    audit_path.unlink(missing_ok=True)  # A failed rerun must not leave a previous success report.
    input_paths=[Path(__file__).resolve(),DATA/'normalized_log2_expression.tsv.gz',
                 DATA/'sample_metadata.tsv',DATA/'simulation_summary.json',
                 ROOT/'gsea_data'/'gsea_provenance.json']+[DATA/spec[4] for spec in SPECS]
    input_hashes={str(path.relative_to(ROOT.parents[1])):hashlib.sha256(path.read_bytes()).hexdigest()
                  for path in input_paths}
    fonts=prepare_nimbus()
    configure_theme()
    tables=[read_table(DATA/spec[4]).set_index('gene_id') for spec in SPECS]
    from run_gsea import validate_completed_run
    validate_completed_run(ROOT/'gsea_data',DATA)
    expression=read_table(DATA/'normalized_log2_expression.tsv.gz').set_index('gene_id')
    metadata=read_table(DATA/'sample_metadata.tsv').set_index('sample_id')
    validate_inputs(tables,expression,metadata)
    summary=json.loads((DATA/'simulation_summary.json').read_text())
    df=summary['analysis']['residual_degrees_of_freedom']
    if not np.isfinite(df) or df<=0:
        raise ValueError('Residual degrees of freedom must be positive and finite.')
    xmax,ymax=volcano_limits(tables)
    crit=stats.t.ppf(.975,df)
    forests=[selected_rows(t)[1] for t in tables]
    forest_limit=max([1.0]+[float(np.ceil((f.log2FoldChange.abs()+crit*f.lfcSE).max()/.2)*.2)
                            for f in forests if len(f)])
    displays=[displayed_gsea(read_table(path))
              for path in sorted((ROOT/'gsea_data').glob('figure_??_*_display.tsv'))]
    nes_limit=max([2.6]+[float(t.NES.abs().max())*1.1 for t in displays if len(t)])
    audits=[render_question(spec,table,expression,metadata,fonts,xmax,ymax,forest_limit,df,nes_limit)
            for spec,table in zip(SPECS,tables)]
    for path in input_paths:
        if hashlib.sha256(path.read_bytes()).hexdigest()!=input_hashes[str(path.relative_to(ROOT.parents[1]))]:
            raise RuntimeError(f'Render input changed while drawing figures: {path}. Rerun rendering.')
    report={'complete':True,'input_sha256':input_hashes,'figures':audits}
    with tempfile.NamedTemporaryFile(mode='w',encoding='utf-8',dir=ROOT,suffix='.json',delete=False) as handle:
        pending=Path(handle.name)
        try:
            handle.write(json.dumps(report,indent=2)+'\n')
            handle.flush()
            os.fsync(handle.fileno())
        except Exception:
            pending.unlink(missing_ok=True)
            raise
    try:
        os.replace(pending,audit_path)
    finally:
        pending.unlink(missing_ok=True)
    print(json.dumps([{'figure':a['figure'],'comparison':a['comparison'],'panels':a['panels'],
                       'pdf_fonts':a['pdf_fonts'],'heatmap_shape':[len(a['heatmap_genes']),len(a['heatmap_samples'])]}
                      for a in audits],indent=2))


if __name__=='__main__':main()
