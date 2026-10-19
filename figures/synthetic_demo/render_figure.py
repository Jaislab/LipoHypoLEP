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
from matplotlib.text import Text
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
GROUPS = ('LD SALINE','LD LEPTIN','NON-LD control')
PALETTE = {'LD SALINE':'#C28E85','LD LEPTIN':'#80A596','NON-LD control':'#8B9FB3'}
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

SPECS = (
    (1, 'How does untreated LD differ from non-LD controls?',
     'LD SALINE', 'NON-LD control', 'contrast_LD_SALINE_vs_NON_LD_control.tsv.gz'),
    (2, 'What does leptin treatment change in LD mice?',
     'LD LEPTIN', 'LD SALINE', 'contrast_LD_LEPTIN_vs_LD_SALINE.tsv.gz'),
    (3, 'What differences remain after leptin treatment?',
     'LD LEPTIN', 'NON-LD control', 'contrast_LD_LEPTIN_vs_NON_LD_control.tsv.gz'),
)


def configure_theme():
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
    selected=table[(table.padj<.05)&(table.log2FoldChange.abs()>=1)]
    heatmap=selected.sort_values('padj',kind='stable').head(40)
    positive=selected[selected.log2FoldChange>0].sort_values('padj',kind='stable').head(6)
    negative=selected[selected.log2FoldChange<0].sort_values('padj',kind='stable').head(6)
    forest=pd.concat([positive,negative]).sort_values('log2FoldChange',ascending=False,kind='stable')
    if len(heatmap)!=40 or len(forest)!=12:
        raise RuntimeError('Synthetic demonstration needs 40 selected heatmap genes and six genes in each direction.')
    return heatmap,forest


def export_figure(fig,number,question,fonts,details):
    name=f'figure_{number:02d}_synthetic'
    labels=[]
    for t in fig.findobj(Text):
        if not t.get_text():continue
        style='Bold' if (t.get_fontweight() in ['bold',700] or
                         t.get_fontproperties().get_file()==str(fonts['Bold'])) else 'Regular'
        t.set_fontproperties(fm.FontProperties(fname=str(fonts[style]),size=t.get_fontsize()))
        labels.append({'text':t.get_text(),'font':f'Nimbus Sans {style}',
                       'size_px':round(t.get_fontsize()*DPI/72,2)})
    fig.canvas.draw()
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
    assert '<text' not in svg.read_text()
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
    assert Image.open(OUT/f'{name}.png').size==(W,H)
    with fitz.open(OUT/f'{name}.pdf') as doc:
        pdf_fonts=sorted({f[3] for page in doc for f in page.get_fonts()})
        if not pdf_fonts or not all('NimbusSans' in f for f in pdf_fonts):
            raise RuntimeError(f'Font substitution in PDF: {pdf_fonts}')
    return {'figure':number,'file_stem':name,'question':question,
            'all_visible_typography':'Nimbus Sans','canvas_pixels':[W,H],
            'pdf_fonts':pdf_fonts,'pdf_text':'embedded, searchable vector text',
            'svg_text':'outlined paths','editable_svg':'text with embedded Nimbus Sans font data',
            'figure_contains_editorial_prose':False,'pca_panel':False,'labels':labels,**details}


def render_question(spec,table,expression,metadata,fonts,xmax,ymax,forest_limit,df):
    number,question,first,second,filename=spec
    heatmap,forest=selected_rows(table)
    comparison=f'{first} vs {second}'
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
    ax=axes(110,88,560,386)
    panel('a',40,27)
    ax.set_title(comparison,pad=20)
    x=table.log2FoldChange.to_numpy()
    y=-np.log10(table.padj.to_numpy())
    selected=(table.padj.to_numpy()<.05)&(np.abs(x)>=1)
    ax.scatter(x[~selected],y[~selected],s=8,c=NS,alpha=.55,edgecolors='none',rasterized=True)
    for mask,color in [(selected&(x<0),LOWER),(selected&(x>0),HIGHER)]:
        ax.scatter(x[mask],y[mask],s=11,c=color,alpha=.82,edgecolors='none',rasterized=True)
    for value in [-1,1]:ax.axvline(value,color=GRID,lw=.65,ls=(0,(3,4)),zorder=0)
    ax.axhline(-np.log10(.05),color=GRID,lw=.65,ls=(0,(3,4)),zorder=0)
    ax.set(xlim=(-xmax,xmax),ylim=(-.08,ymax),xticks=[-4,-2,0,2,4],yticks=np.arange(0,ymax+.1,2))
    ax.set_xlabel(r'log$_2$ fold change',labelpad=8)
    ax.set_ylabel(r'−log$_{10}$(adjusted $P$)',labelpad=8)
    theme(ax)

    # b: only the two groups that answer this figure's question.
    samples=[s for group in (first,second) for s in metadata.index[metadata.group==group]]
    matrix=expression.loc[heatmap.index,samples].to_numpy()
    z=(matrix-matrix.mean(axis=1,keepdims=True))/matrix.std(axis=1,keepdims=True,ddof=1)
    rows=leaves_list(linkage(z,method='average',metric='euclidean',optimal_ordering=True))
    z=z[rows]
    genes=heatmap.index[rows]
    panel('b',788,27)
    ax=axes(928,132,450,772)
    cmap=LinearSegmentedColormap.from_list('soft_expression',[LOWER,'#FAF9F6',HIGHER])
    im=ax.imshow(z,aspect='auto',interpolation='nearest',cmap=cmap,vmin=-2,vmax=2)
    ax.set_xticks([]);ax.set_yticks([])
    ax.set_ylabel('Genes',labelpad=12)
    for spine in ax.spines.values():spine.set_visible(False)
    ax.axvline(5.5,color='white',lw=2)
    for i,group in enumerate((first,second)):
        left=928+i*225
        strip=axes(left,120,223,6)
        strip.add_patch(plt.Rectangle((0,0),1,1,color=PALETTE[group]))
        strip.set_axis_off()
        fig.text((left+112.5)/W,1-108/H,group,ha='center',va='bottom',fontsize=px(24))
    cax=axes(1402,386,13,264)
    cb=fig.colorbar(im,cax=cax,ticks=[-2,0,2])
    cb.outline.set_visible(False)
    cb.ax.tick_params(size=2,pad=4,labelsize=px(21))
    cb.set_label('z-score',fontsize=px(24),labelpad=7)
    pd.DataFrame({'gene_id':genes}).to_csv(DATA/f'figure_{number:02d}_heatmap_genes.tsv',sep='\t',index=False)

    # c: selected effect estimates and marginal, unadjusted model intervals.
    panel('c',40,564)
    ax=axes(238,608,432,296)
    crit=stats.t.ppf(.975,df)
    errors=crit*forest.lfcSE.to_numpy()
    fx=forest.log2FoldChange.to_numpy()
    for row,(effect,error) in enumerate(zip(fx,errors)):
        color=HIGHER if effect>0 else LOWER
        ax.errorbar(effect,row,xerr=error,fmt='o',ms=5.5,mfc=color,mec='white',mew=.5,
                    color=color,ecolor=color,elinewidth=1.0,capsize=2,capthick=.8)
    ax.axvline(0,color=GRID,lw=.8,ls=(0,(3,4)),zorder=0)
    ax.set_yticks(np.arange(len(forest)),labels=forest.index)
    ax.set(xlim=(-forest_limit,forest_limit),ylim=(len(forest)-.5,-.5),xticks=[-4,-2,0,2,4])
    ax.set_xlabel(r'log$_2$ fold change',labelpad=9)
    ax.tick_params(axis='y',length=0,pad=10)
    ax.spines['left'].set_visible(False)
    theme(ax)
    ax.tick_params(axis='y',length=0,pad=10)
    forest=forest.copy()
    forest['ci_lower']=fx-errors
    forest['ci_upper']=fx+errors
    forest['comparison']=comparison
    forest.reset_index().to_csv(DATA/f'figure_{number:02d}_forest_genes.tsv',sep='\t',index=False,float_format='%.17g')
    return export_figure(fig,number,question,fonts,{
        'comparison':comparison,'contrast_table':filename,
        'panels':['volcano','two-group heatmap','selected effects with 95% confidence intervals'],
        'volcano_limits':{'absolute_log2fc':xmax,'negative_log10_adjusted_p':ymax},
        'volcano_genes':len(table),'heatmap_samples':samples,'heatmap_genes':genes.tolist(),
        'heatmap_group_order':[first,second],
        'heatmap_selection':'40 display-threshold genes with smallest adjusted P in this contrast, clustered using only the two displayed groups',
        'forest_genes':forest.index.tolist(),'forest_absolute_limit':forest_limit,
        'forest_confidence_level':.95,'forest_residual_degrees_of_freedom':df,
        'forest_intervals':'marginal, unadjusted; genes selected from the same dataset',
    })


def main():
    fonts=prepare_nimbus()
    configure_theme()
    tables=[pd.read_csv(DATA/spec[4],sep='\t').set_index('gene_id') for spec in SPECS]
    for table in tables[1:]:assert table.index.equals(tables[0].index)
    expression=pd.read_csv(DATA/'normalized_log2_expression.tsv.gz',sep='\t').set_index('gene_id')
    metadata=pd.read_csv(DATA/'sample_metadata.tsv',sep='\t').set_index('sample_id')
    summary=json.loads((DATA/'simulation_summary.json').read_text())
    df=summary['analysis']['residual_degrees_of_freedom']
    xmax=float(np.ceil(max(t.log2FoldChange.abs().max() for t in tables)))
    ymax=float(np.ceil(max((-np.log10(t.padj)).max() for t in tables)/2)*2)
    crit=stats.t.ppf(.975,df)
    forest_limit=float(np.ceil(max((f.log2FoldChange.abs()+crit*f.lfcSE).max()
                                  for f in [selected_rows(t)[1] for t in tables])))
    audits=[render_question(spec,table,expression,metadata,fonts,xmax,ymax,forest_limit,df)
            for spec,table in zip(SPECS,tables)]
    (ROOT/'font_and_plot_audit.json').write_text(json.dumps({'figures':audits},indent=2)+'\n')
    print(json.dumps([{'figure':a['figure'],'comparison':a['comparison'],'panels':a['panels'],
                       'pdf_fonts':a['pdf_fonts'],'heatmap_shape':[len(a['heatmap_genes']),len(a['heatmap_samples'])]}
                      for a in audits],indent=2))


if __name__=='__main__':main()
