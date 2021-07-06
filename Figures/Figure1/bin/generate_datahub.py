# Generate data hub for WashU Epigenome browser

import json
import glob
import os
import numpy as np
from matplotlib import cm
from matplotlib.colors import rgb2hex

datahub = []

colors={"egg":"#207AB5",
"fertilized_egg":"#207AB5",
"1-cell":"#207AB5",
"16-cell":"#3B77AC",
"32-cell":"#3B77AC",
"64-cell":"#4975A6",
"128-cell":"#5C709B",
"256-cell":"#656E94",
"512-cell":"#786886",
"1k-cell":"#81647D",
"high":"",
"oblong":"#833B88",
"sphere":"#896076",
"dome":"#985769",
"dome-30per_epiboly":"#A64D5B",
"30per_epiboly":"#A64D5B",
"germ_ring":"#AC4754",
"shield":"#C62437",
"75per_epiboly":"#C62437",
"1-4_somites":"#C62437",
"5-8_somites":"#C62437",
"14-19_somites":"#C62437",
"prim-5":"#C62437",
"prim-25":"#C62437",
"high-pec":"#C62437",
"long-pec":"#C62437",
"protruding_mouth":"#C62437",
"90_Days-2_Years":"#C62437"}

color_paired = [
"#A6CEE3",
"#1F78B4",
"#33A02C",
"#B2DF8A",
"#E31A1C",
"#FB9A99",
"#FF7F00",
"#6A3D9A",
"#CAB2D6",
"#bdbdbd"]

assays = ["RNA_pos", "RNA_neg",
            "nanti_CAGE_pos","nanti_CAGE_neg",
            "new_CAGE_pos","new_CAGE_neg",
            "CAGE_pos", "CAGE_neg",
            "Methyl","H3K4me1", "H3K4me3",
            "4C_six2a","ATAC", "H3K27ac", "H3K27me3"]
for assay in assays:
    with open('shifted/' + assay + '/log.json') as f:
        data = json.load(f)
    bws = []
    data.reverse()

    # generate individual bigwig tracks
    for i, bw in enumerate(data):
        stage = bw['stage']
        track = {
            "type": 'bigwig',
            "name": stage,
            "url": 'https://export.uppmax.uu.se/uppstore2017255/danio-code_trackhub/' + bw['new_file'],
            "options": {
                "color": colors[stage]
            },
            "metadata": {
                "stage": stage,
                "assay": bw['assay'],
                "file": bw['new_file'],
            }
        }
        bws.append(track)
    # add individual bigwig files to matplot track
    datahub.append({
        "type": "matplot",
        "name": bw['assay'],
        "tracks": bws,
        "options": {
            "height": 120
        },
        "metadata": {
            "assay": bw['assay']
        }
    })

    for i, bw in enumerate(data):
        datahub.append({
            "type": 'bigwig',
            "name": bw['stage'],
            "url": 'https://export.uppmax.uu.se/uppstore2017255/danio-code_trackhub/' + bw['og_file'],
            "options": {
                "color": colors[stage]
            },
            "metadata": {
                "stage": bw['stage'],
                "assay": bw['assay'],
                "displayMode": "heatmap" if assay == "H3K27ac" else "bar"
            }
        })

for assay in ["chromHMM","pdre"]:
    categories = {}
    for i,cat in enumerate(["1_TssA1","2_TssA2","3_TssFlank1","4_TssFlank2","5_EnhA1","6_EnhFlank","7_EnhWk1","8_Pois","9_ReprPC","10_Quies"]):
        categories[cat]={"name":cat,"color":color_paired[i]}
    for file in sorted(glob.glob('shifted/'+assay+'/*.gz')):
        stage = os.path.basename(file).split(".")[1]

        datahub.append({
            "type":"categorical",
            "name": stage+"_"+assay,
            "url": 'https://export.uppmax.uu.se/uppstore2017255/danio-code_trackhub/'+file,
            "options":{
                "height": 12,
                "category":categories
            },
            "metadata": {
                "stage": stage,
                "assay": assay
            }
        })

for assay in ["tissue_specific"]:
    categories = {
        "state_1":{"name":"state_1","color":"#A3539D"},
        "state_3":{"name":"state_3","color":"#70BC4C"},
        "state_4":{"name":"state_4","color":"#EC3C2C"},
        "state_5":{"name":"state_5","color":"#EB2237"},
        "state_10":{"name":"state_10","color":"#6FBB4C"},
        "state_12":{"name":"state_12","color":"#BF509C"},
        "state_14":{"name":"state_14","color":"#F1652F"},
        "state_16":{"name":"state_16","color":"#D84897"},
        "state_18":{"name":"state_18","color":"#F59231"},
        "state_20":{"name":"state_20","color":"#4156A1"},
        "state_21":{"name":"state_21","color":"#436FB4"},
        "state_23":{"name":"state_23","color":"#6FCCDC"},
        "state_24":{"name":"state_24","color":"#FDC331"},
        "state_25":{"name":"state_25","color":"#F7EC3C"},
        "state_26":{"name":"state_26","color":"#DADE39"},
        "state_29":{"name":"state_29","color":"#EC252D"}
    }

    for file in sorted(glob.glob('shifted/'+assay+'/*.gz')):
        state = os.path.basename(file).split(".")[0]

        datahub.append({
            "type":"categorical",
            "name": state+"_"+assay,
            "url": 'https://export.uppmax.uu.se/uppstore2017255/danio-code_trackhub/'+file,
            "options":{
                "height": 12,
                "category":categories
            },
            "metadata": {
                "state": state,
                "assay": assay
            }
        })

for assay in ['mm_H3K27me3']:
    for file in sorted(glob.glob('shifted/'+assay+'/*.bw')):
        stage = "embryo"
        track = {
            "type":"bigwig",
            "name": "embryo",
            "url": 'https://export.uppmax.uu.se/uppstore2017255/danio-code_trackhub/'+file,
            "options":{
                "color":"#1F7AB4"
            },
            "metadata": {
                "stage": stage,
                "assay": assay
            }
        }

        datahub.append({
        "type": "matplot",
        "name": assay,
        "tracks": [track],
        "options": {
            "height":60
        },
        "metadata": {
            "assay": assay
        }
    })
        datahub.append({

        })

outfile = open('/proj/uppstore2017255/webexport/danio-code_trackhub/datahub.json', 'w')
datahub = json.dumps(datahub, indent=4)
outfile.write(datahub)

outfile.close()
