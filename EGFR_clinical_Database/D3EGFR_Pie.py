#coded by ylshi
#2020-12-29
#pip install pyecharts

from pyecharts.charts import Pie
import pyecharts.options as opts
import sys

drug=sys.argv[1]
CR=sys.argv[2]
PR=sys.argv[3]
SD=sys.argv[4]
PD=sys.argv[5]

attr2 = ['CR', 'PR', 'SD', 'PD']
# Good response = CR + PR + SD
# Bad response = PD
v2 = [CR,PR,SD,PD]
data2 = [(attr2[i],v2[i]) for i in range(len(attr2))]

#pie = (Pie(init_opts=opts.InitOpts(width="300px",
#                            height="350px",
pie = (Pie(init_opts=opts.InitOpts(page_title='EGFR-mutated clinical database'))
.add("clinical response",data2,radius = ["0%","70%"],
#label_opts=opts.LabelOpts(position="inside",font_size=25,formatter="{b}:{d}%")))
label_opts=opts.LabelOpts(position='insideleft', font_size=28,formatter="{b}:{d}%")))
pie.set_global_opts(title_opts=opts.TitleOpts(title=drug,pos_left='center',title_textstyle_opts=(opts.TextStyleOpts(font_size=32))),legend_opts=opts.LegendOpts(is_show=False))
#pie.render_notebook()
pie.render("D3EGFR_clinical_statistics.html")
