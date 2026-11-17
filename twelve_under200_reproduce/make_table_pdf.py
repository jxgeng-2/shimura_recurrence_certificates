"""One consistent table for twelve optimal Shimura elliptic factors.

sage -python make_table_pdf.py
Only consumes exact exported coefficients; it does not recognize numbers.
"""
from pathlib import Path
import argparse
import html
import json
import re
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak

PAGE_SIZE=landscape(A4)
TABLE_WIDTH=PAGE_SIZE[0]-92  # Page margins plus the frame's six-point padding.


def register_serif_fonts():
    candidates=[(Path('C:/Windows/Fonts'),['times.ttf','timesi.ttf','timesbd.ttf','timesbi.ttf']),
                (Path('/usr/share/fonts/truetype/dejavu'),['DejaVuSerif.ttf','DejaVuSerif-Italic.ttf',
                  'DejaVuSerif-Bold.ttf','DejaVuSerif-BoldItalic.ttf'])]
    try:
        import matplotlib
        candidates.append((Path(matplotlib.get_data_path())/'fonts'/'ttf',candidates[1][1]))
    except ImportError:pass
    base,files=next((base,files) for base,files in candidates if all((base/f).exists() for f in files))
    for name,file in zip(['Serif','SerifI','SerifB','SerifBI'],files):
        pdfmetrics.registerFont(TTFont(name,str(base/file)))
    pdfmetrics.registerFontFamily('Serif',normal='Serif',italic='SerifI',bold='SerifB',boldItalic='SerifBI')


def polynomial_markup(coeffs):
    pieces=[]
    for power in range(len(coeffs)-1,-1,-1):
        c=coeffs[power]
        if not c:continue
        sign=('-' if c<0 else '') if not pieces else (' -&nbsp;' if c<0 else ' +&nbsp;')
        v=abs(c)
        term=str(v) if power==0 or v!=1 else ''
        if power:
            term+='<i>x</i>'
            if power!=1:term+='<super>%s</super>'%power
        pieces.append(sign+term)
    return ''.join(pieces) or '0'


def value_markup(row):
    if row['status']=='identity':return 'point at infinity'
    if 'rational_numerator' in row:
        n,d=row['rational_numerator'],row['rational_denominator']
        text='<i>x</i> = '+str(n)+(' / '+str(d) if d!=1 else '')
    else:text=polynomial_markup(row['polynomial_coefficients'])
    if row.get('annotation'):text+=' &nbsp;('+html.escape(row['annotation'])+')'
    return text


def build_table(curve,font_size):
    base=ParagraphStyle('Cell',fontName='Serif',fontSize=font_size,leading=font_size*1.23,
                        alignment=TA_LEFT,spaceBefore=0,spaceAfter=0)
    centered=ParagraphStyle('Centered',parent=base,alignment=TA_CENTER)
    header=ParagraphStyle('Head',parent=base,fontSize=font_size+1,leading=font_size*1.3,alignment=TA_CENTER)
    P=lambda text,style=base:Paragraph(text,style)
    rows=[[P('<i>E</i>',header),P('<i>D</i>',header),P('<i>h</i>',header),P('<i>τ</i>',header),
           P('minimal polynomial of <i>x</i>(<i>P</i>(<i>D</i>, τ))',header)]]
    commands=[('GRID',(0,0),(-1,-1),.55,colors.black),('VALIGN',(0,0),(-1,-1),'TOP'),
              ('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),
              ('TOPPADDING',(0,0),(-1,-1),3.5),('BOTTOMPADDING',(0,0),(-1,-1),3.5),
              ('TOPPADDING',(0,0),(-1,0),5),('BOTTOMPADDING',(0,0),(-1,0),5)]
    groups={}
    for i,row in enumerate(curve['rows'],1):
        first=row['D'] not in groups
        groups.setdefault(row['D'],[]).append(i)
        label=curve['E']
        value=value_markup(row)
        # A long integer may move to the next line, but must never split digits.
        for digits in re.findall(r'\d+',value):
            assert pdfmetrics.stringWidth(digits,'Serif',font_size)<TABLE_WIDTH-147-16, (label,row['D'],'integer too wide')
        rows.append([P('%s<i>%s</i><sub>%s</sub>'%(curve['N'],label[-2].upper(),label[-1]),centered) if i==1 else '',
                     P(str(row['D']),centered) if first else '',P(str(row['h']),centered) if first else '',
                     P('τ<sub>%s</sub>'%row['tau'],centered),P(value)])
    if len(rows)>2:commands.append(('SPAN',(0,1),(0,len(rows)-1)))
    for indices in groups.values():
        for col in (1,2):commands.append(('SPAN',(col,indices[0]),(col,indices[-1])))
    table=Table(rows,colWidths=[48,37,24,38,TABLE_WIDTH-147],repeatRows=1,hAlign='LEFT')
    table.setStyle(TableStyle(commands))
    return table


def paginated_tables(curve,cap,continued_cap=440):
    """Keep each discriminant together; balance long tables across pages."""
    for font in (9.8,):
        table=build_table(curve,font);_,height=table.wrap(TABLE_WIDTH,10000)
        if height<=cap:return [(table,font,height,len(curve['rows']))]
    groups=[]
    for row in curve['rows']:
        if not groups or groups[-1][0]['D']!=row['D']:groups.append([])
        groups[-1].append(row)
    font=9.8;memo={}
    def segment(a,b):
        if (a,b) not in memo:
            subset=[r for group in groups[a:b] for r in group]
            table=build_table(dict(curve,rows=subset),font)
            _,height=table.wrap(TABLE_WIDTH,10000)
            memo[a,b]=(table,font,height,len(subset))
        return memo[a,b]
    # Dynamic programming avoids exponential cut enumeration for longer ranges.
    # First minimize page count, then balance whitespace without splitting D.
    count=len(groups)
    best={0:(0,0.0,[])}
    for end in range(1,count+1):
        choices=[]
        for start in range(end):
            if start not in best:continue
            part=segment(start,end)
            page_cap=cap if start==0 else continued_cap
            if part[2]>page_cap:continue
            pages,cost,parts=best[start]
            choices.append((pages+1,cost+(page_cap-part[2])**2,parts+[part]))
        if choices:best[end]=min(choices,key=lambda x:(x[0],x[1]))
    if count in best:return best[count][2]
    raise ValueError('A discriminant group is too tall to fit: '+curve['E'])


def make_pdf(data_path,output_path):
    register_serif_fonts()
    data=json.loads(Path(data_path).read_text(encoding='utf-8'))
    output_path=Path(output_path);output_path.parent.mkdir(parents=True,exist_ok=True)
    doc=SimpleDocTemplate(str(output_path),pagesize=PAGE_SIZE,leftMargin=40,rightMargin=40,
        topMargin=35,bottomMargin=35,title='CM divisor images on twelve optimal elliptic curves',
        author='Sage and Magma computations')
    title=ParagraphStyle('Title',fontName='SerifB',fontSize=16,leading=20,spaceAfter=7)
    body=ParagraphStyle('Body',fontName='Serif',fontSize=8.2,leading=11.1,spaceAfter=6)
    curvehead=ParagraphStyle('Curve',fontName='SerifB',fontSize=12.5,leading=16,spaceAfter=5)
    note=ParagraphStyle('Note',parent=body,fontSize=7.6,leading=10,spaceBefore=6)
    story=[];layout=[]
    for index,curve in enumerate(data['curves']):
        if index:story.append(PageBreak())
        if index==0:
            story.append(Paragraph('CM divisors on twelve optimal elliptic curves',title))
            scope='All admissible negative fundamental <i>D</i> with |<i>D</i>| &lt; 200; no class-number cutoff. '
            scope+='All ideal classes, local orientations, and three nonidentity Atkin-Lehner operators. '
            scope+='Targets use Sage\'s <i>X</i><sub>0</sub>(<i>N</i>)-optimal convention.'
            scope+=' Exact curve equations and inclusions into <i>H</i><sub>D</sub> are checked; analytic identification is numerical, not interval-certified.'
            if not data['complete']:scope='<b>LAYOUT PREVIEW - incomplete curve set.</b> '+scope
            story.append(Paragraph(scope,body))
        story.append(Paragraph('X(%s,1) &nbsp; / &nbsp; %s'%(curve['N'],curve['E']),curvehead))
        meta=('Base discriminant %s; %s and %s terms at %s bits. '
              '%s divisors: %s finite points, %s identities. Minimum finite coordinate match: %s bits.'%
              (curve['base_discriminant'],*curve['truncations'],curve['prec_bits'],curve['individual_divisors'],
               curve['finite_divisors'],curve['identity_divisors'],curve['min_xy_match_bits']))
        if 'kernel' in curve.get('evaluation_method',''):
            meta+=' Verified kernel-AL preconditioning is used for selected evaluations.'
        if 'affine' in curve.get('evaluation_method',''):
            meta+=' w33 images use P(w33) = P(w3) + (0,0), with direct-coordinate checks retained.'
        story.append(Paragraph(meta,body));story.append(Spacer(1,3))
        # Large optimal-model polynomials can require a continuation page.
        cap=365 if index==0 else (405 if 'affine' in curve.get('evaluation_method','') else 425)
        parts=paginated_tables(curve,cap)
        for part_index,(table,font,height,count) in enumerate(parts):
            if part_index:
                story.append(PageBreak())
                story.append(Paragraph('X(%s,1) &nbsp; / &nbsp; %s (continued)'%
                                       (curve['N'],curve['E']),curvehead))
            story.append(table)
            layout.append(dict(level=curve['N'],target=curve['E'],part=part_index+1,
                               font_size=font,table_height=round(height,2),rows=count))
        story.append(Paragraph('Native elliptic quotient: %s. Its isogeny to the displayed optimal target has degree %s.'%
                               (curve['native_quotient'],curve['isogeny_degree']),note))
    def footer(canvas,document):
        canvas.saveState();canvas.setFont('Serif',7.5)
        canvas.drawString(40,20,'CM divisors | Twelve Sage-optimal targets | fundamental |D| < 200')
        canvas.drawRightString(PAGE_SIZE[0]-40,20,str(document.page));canvas.restoreState()
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    output_path.with_suffix('.layout.json').write_text(json.dumps(layout,indent=2),encoding='utf-8')
    print(json.dumps(dict(output=str(output_path),complete=data['complete'],layout=layout),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    default_dir=Path(__file__).resolve().parent/'results'
    parser.add_argument('--data',default=str(default_dir/'cm_divisors_twelve_under200_data.json'))
    parser.add_argument('--output',default=str(default_dir/'cm_divisors_twelve_under200.pdf'))
    args=parser.parse_args();make_pdf(args.data,args.output)
